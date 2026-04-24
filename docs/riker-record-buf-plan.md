# RikerRecordBuf — design notes

Working notes for a follow-on PR. Captures what we learned while profiling
`riker multi`, why a custom record type looks attractive, the shape of the
type, and the open questions to resolve when we come back to this.

## Background: the bottleneck

After landing the pooled `RecordBuf` recycling (PR on `tf_multi_recordbuf_pool_proto`)
and the crossbeam MPMC work queue (same branch), `riker multi --tools wgs`
on a 12× 1KG BAM (~313M reads, 18GB) runs in **~3:01 wall** with **225s user
time**. Sampling profile (`samply`, `bench-prof` profile with LTO off):

| Thread            | CPU     | Top cost                                                    |
|-------------------|---------|-------------------------------------------------------------|
| Reader            | 186.9s  | BGZF decompress 92s · BAM decode 76s · Data churn 18s       |
| Pool workers (×2) | ~34s ea | `wait_timeout` gone after crossbeam; now ~10s real work each |

The reader thread is fully pegged at 98% of wall. Of its 187s:

- **49%  (92s)** — BGZF decompress (libdeflate 72s + CRC 19s). Irreducible
  unless we parallelise BGZF itself.
- **41%  (76s)** — BAM record decode inside `read_record_buf`:
  - **28%  (53s)** — `read_data` / `read_value` / `read_string` — parses
    every BAM aux tag into a `Data` `HashMap`
  - 7%   (13s) — `read_sequence`
  - **6%  (18s)** — `Data::clear` + `Data::insert` — resets and refills the
    aux HashMap on every record (because we reuse one `RecordBuf`)
  - ~5%   (10s) — `read_value`/`read_name`/`decode` misc.
- ~4% (7s) — rpmalloc deallocate, mostly tied to aux tag churn.

**No wgs collector touches aux tags.** The error collector reads NM. Everything
else is pure decode-and-throw-away cost.

## Goal

Eliminate the ~71s of per-record aux-tag work on the reader thread for
collectors that don't need it. Target wall time: **2:00–2:15**. The BGZF floor
sits around 92s; the rest of the reader work is parsable but not deeply
reducible without a decoder rewrite.

Secondary: keep the `name`/`cigar`/`sequence`/`quality` decode on the reader
thread (workers need these), but consider making *them* lazy if profiling
shows the wins are there.

## Proposed design: `RikerRecordBuf`

A record type that's:

1. **Lazy like `bam::Record`** — stores the raw BAM bytes plus computed
   field offsets. Decoded fields are cached on first access.
2. **Reusable like `RecordBuf`** — `read_record_buf(&mut RikerRecordBuf)`
   overwrites the byte buffer in place; no per-record allocation.
3. **Owned and `Send`** — goes behind an `Arc` and fans out across threads.
4. **Drop-in at the collector seam** — collectors call the same accessor
   shapes they call on `RecordBuf` today.

### Internals sketch

```rust
pub struct RikerRecordBuf {
    /// Raw BAM record bytes (without the 4-byte block_size prefix).
    /// Reused across records — `read_into(&mut self, reader)` overwrites
    /// this buffer and recomputes `offsets`.
    buf: Vec<u8>,

    /// Byte offsets into `buf` for each variable-length field, plus
    /// the fixed-layout scalar fields copied out for fast access.
    offsets: FieldOffsets,

    /// Lazily-decoded cigar (walk_depth iterates this often).
    cigar_cache: OnceCell<Cigar>,

    /// Lazily-unpacked sequence (BAM stores it 4-bit-packed; wgs and
    /// error walk bytes, so unpack once and cache).
    seq_cache: OnceCell<Vec<u8>>,

    /// Lazily-decoded aux data — skipped entirely for collectors that
    /// never call `data()`.
    data_cache: OnceCell<Data>,
}

struct FieldOffsets {
    // Variable-length field slices
    name: Range<u32>,
    cigar: Range<u32>,
    sequence: Range<u32>,
    quality: Range<u32>,
    data: Range<u32>,

    // Fixed-layout scalars (copied out once during read)
    flags: Flags,
    ref_id: Option<usize>,
    pos: Option<Position>,
    mapq: Option<MappingQuality>,
    mate_ref_id: Option<usize>,
    mate_pos: Option<Position>,
    tlen: i32,
}
```

### Access patterns

| Method                         | Cost           | Notes                                          |
|--------------------------------|----------------|------------------------------------------------|
| `flags()`                      | copy a u16     | from `offsets`                                 |
| `reference_sequence_id()`      | copy           | from `offsets`                                 |
| `alignment_start()`            | copy           | from `offsets`                                 |
| `mapping_quality()`            | copy           | from `offsets`                                 |
| `mate_reference_sequence_id()` | copy           | from `offsets`                                 |
| `mate_alignment_start()`       | copy           | from `offsets`                                 |
| `template_length()`            | copy           | from `offsets`                                 |
| `name()`                       | slice          | `&buf[offsets.name]`, zero-copy                |
| `quality_scores()`             | slice          | `&buf[offsets.quality]`, zero-copy             |
| `cigar()`                      | decode once    | `OnceCell<Cigar>`                              |
| `sequence()`                   | unpack once    | `OnceCell<Vec<u8>>` — 4-bit → 1-byte expansion |
| `data()`                       | decode once    | `OnceCell<Data>` — **this is the big win**     |

### Reset semantics

Reusing the buffer across records requires clearing the `OnceCell`s. `OnceCell`
doesn't have a public `reset()`, but since we own the struct we can `std::mem::take`
each cache and drop it, which releases the decoded field back to the pool/heap.
Alternatively, wrap the caches in our own `Cell<Option<T>>`-like primitive with
a `reset()` method and an invalidation counter — the ergonomic choice matters
less than getting the invariant right: **after `read_into`, every cache is
empty.**

## Where the wins come from

| Current cost              | Fate under `RikerRecordBuf`                        |
|---------------------------|----------------------------------------------------|
| 53s aux-tag parse         | Deferred; zero if no collector calls `data()`       |
| 18s `Data::clear`/`insert`| No HashMap in the hot path                          |
| ~7s rpmalloc              | Mostly gone (Data churn was a big chunk)            |
| 13s sequence unpack       | Still paid when wgs/error asks; same cost but once  |
| 92s BGZF                  | Unchanged — fundamental                             |
| 76s BAM decode overall    | Drops roughly to `92 + (per-field cache fills)`     |

Back-of-envelope: reader falls from 187s → ~110–120s CPU → wall ~2:00–2:15.

## Migration path (incremental)

1. **`RikerRecordBuf` alongside `RecordBuf`** — introduce the type with the
   exact method surface `wgs`' `walk_depth`, `build_non_n_bitvec`, and the
   `MateBuffer` need. Do *not* try to make it `impl sam::alignment::Record`
   immediately; just match `RecordBuf`'s inherent method signatures where
   collectors currently call them.

2. **`AlignmentReader::read_riker_record_buf`** — new in-place read method
   that populates `RikerRecordBuf` directly from BAM: memcpy the record body
   into `buf`, compute offsets, populate the fixed-layout scalars, skip
   aux-tag parsing entirely.

3. **Swap the `multi` reader thread** — change `Batch`'s inner type from
   `Vec<RecordBuf>` to `Vec<RikerRecordBuf>`, and thread the type through
   `accept_multiple`. Collectors now see `&RikerRecordBuf`.

4. **Collector trait change** — `Collector::accept` signature becomes
   `fn accept(&mut self, record: &RikerRecordBuf, header: &Header) -> Result<()>`.
   Every collector edited; mostly mechanical. `MateBuffer<T>` keeps its `T`
   generic, so its internal cache story doesn't change, but callers that
   cached `RecordBuf` (error) switch to caching `RikerRecordBuf`.

5. **Standalone paths** — commands that don't go through `multi` (`wgs`,
   `error`, etc.) also migrate to `RikerRecordBuf` by calling the same
   reader method.

6. **Benchmark, then decide whether to go further** — if reader still has
   hot spots (sequence unpack, cigar decode), those become the next targets.

## Open questions

- **Trait vs. concrete?** Start concrete (`RikerRecordBuf` directly). Adding
  a `RikerRecord` trait later is easy if we ever need to accept noodles'
  types polymorphically.
- **CRAM fallback.** Noodles exposes no lazy CRAM record, so we'd eagerly
  build a `RikerRecordBuf` from a `RecordBuf` for CRAM inputs. That's a
  performance regression for CRAM vs. current behaviour (we'd decode aux
  tags eagerly then discard them unless needed). Acceptable: CRAM is the
  minority case and the eager path is what we have today.
- **Where does `data()` caching live when a collector *does* want it?** For
  `error`'s NM-tag lookup, a `HashMap<Tag, Value>` is overkill — we parse
  the whole aux block to read one tag. A narrower helper like
  `find_aux_tag(record, Tag::NM) -> Option<Value>` that scans the raw aux
  bytes without building a HashMap might be even faster for single-tag
  lookups. Worth measuring once the basic path works.
- **Mate-buffer storage.** The error collector's `MateBuffer<RecordBuf>`
  clones the record into the cache. Under `RikerRecordBuf` the `buf: Vec<u8>`
  means the clone still allocates for the byte buffer itself. We could
  either accept this (paid only for buffered reads, not all reads) or
  extend the buffer pool story to reuse those buffers too — but that adds
  complexity. Start with the simple clone.
- **Handling of `l_read_name`, `l_seq`, `n_cigar_op`.** These are fixed
  offsets in the BAM record header; straightforward to pull out during
  `read_into` and cache in `FieldOffsets`.

## Risks

- **Invariant fragility** — the `OnceCell` reset on buffer reuse is a
  correctness land-mine. A forgotten reset means a cached field returned
  from a previous record would be served for the next. Build a debug-build
  invariant check (e.g. a generation counter on the buffer that the caches
  compare against, panicking if they see a stale generation).
- **Maintenance overhead** — custom record type is another thing to
  understand when contributing. We should document the BAM field layout
  and the reset semantics thoroughly in the module docs.
- **Binary compatibility across noodles updates** — our BAM field-offset
  parsing effectively hard-codes noodles' current encoding. A noodles
  bump that changes an internal layout detail could silently break us.
  We should keep a property-test round-tripping `RikerRecordBuf` against
  `RecordBuf` decoded from the same bytes.

## Starting point

Cheapest experiment that validates the hypothesis:

1. Create `RikerRecordBuf` with **only `data` made lazy** — keep `cigar`,
   `sequence`, `quality` eagerly decoded (so step 1 is strictly smaller
   than the final design).
2. Wire it into `multi`'s reader thread for wgs only.
3. Time it.

If reader CPU drops by ~50s as predicted, expand to the full lazy design.
If it doesn't, something's wrong with the model and we should stop and
re-profile before investing more.

## References

- Profile: `/tmp/hg03953-multi-prof2.json.gz` (captured against `1aacac6`).
- Current reader: `src/commands/multi.rs::reader_thread_loop`.
- Current batch wrapper: `src/commands/multi.rs::RecyclableBatch`.
- Collector trait: `src/collector.rs`.
- Noodles BAM decoder: `noodles-bam-0.88.0::record::codec::decoder`.

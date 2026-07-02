//! Integration tests for the `rna` command.
//!
//! All gene-model and BAM inputs are built programmatically (no checked-in files).

mod helpers;

use std::path::Path;

use helpers::{SamBuilder, SortOrder, read_metrics_tsv};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record_buf::data::field::Value as DataValue;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use riker_lib::commands::command::Command;
use riker_lib::commands::common::{InputOptions, OptionalReferenceOptions, OutputOptions};
use riker_lib::commands::rna::{
    BIOTYPE_PLOT_SUFFIX, BIOTYPE_SUFFIX, COVERAGE_PLOT_SUFFIX, GENE_EXPRESSION_PLOT_SUFFIX,
    ISIZE_SUFFIX, METRICS_SUFFIX, Rna, RnaBiotypeMetric, RnaInsertSizeMetric, RnaOptions,
    RnaSeqMetric, StrandSpec,
};
use tempfile::{NamedTempFile, TempDir};

// ─── Helpers ───────────────────────────────────────────────────────────────────

/// Write `content` to a temp file and return the handle (kept alive by the caller).
fn temp_with(content: &str) -> NamedTempFile {
    let f = NamedTempFile::new().unwrap();
    std::fs::write(f.path(), content).unwrap();
    f
}

/// Run `rna` over `bam` with `gene_model`, returning the output directory and the single
/// summary-metrics row.
fn run_rna(
    bam: &Path,
    gene_model: &Path,
    mutate: impl FnOnce(&mut RnaOptions),
) -> (TempDir, RnaSeqMetric) {
    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    let mut options = RnaOptions { gene_model: gene_model.to_path_buf(), ..Default::default() };
    mutate(&mut options);
    let cmd = Rna {
        input: InputOptions { input: bam.to_path_buf() },
        output: OutputOptions { output: prefix.clone() },
        reference: OptionalReferenceOptions { reference: None },
        options,
    };
    cmd.execute(None).unwrap();
    let mut metrics: Vec<RnaSeqMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{METRICS_SUFFIX}"))).unwrap();
    assert_eq!(metrics.len(), 1, "expected exactly one summary row");
    (dir, metrics.pop().unwrap())
}

/// A two-exon coding gene `MYGENE` on chr1, positive strand, in each of the three formats.
/// 1-based inclusive: exon1 [1000,2000], exon2 [3000,4000], CDS [1500,3500].
/// So coding = 1500..2000 and 3000..3500; UTR = 1000..1499 and 3501..4000; intron = 2001..2999.
fn refflat_two_exon() -> NamedTempFile {
    // refFlat is 0-based half-open: txStart=999, txEnd=4000, cdsStart=1499, cdsEnd=3500,
    // exonStarts=999,2999 exonEnds=2000,4000.
    temp_with("MYGENE\tNM_1\tchr1\t+\t999\t4000\t1499\t3500\t2\t999,2999,\t2000,4000,\n")
}

fn gtf_two_exon() -> NamedTempFile {
    let attrs = "gene_id \"G1\"; transcript_id \"NM_1\"; gene_name \"MYGENE\"; gene_type \"protein_coding\";";
    temp_with(&format!(
        "chr1\tx\texon\t1000\t2000\t.\t+\t.\t{attrs}\n\
         chr1\tx\tCDS\t1500\t2000\t.\t+\t0\t{attrs}\n\
         chr1\tx\texon\t3000\t4000\t.\t+\t.\t{attrs}\n\
         chr1\tx\tCDS\t3000\t3500\t.\t+\t0\t{attrs}\n"
    ))
}

fn gff3_two_exon() -> NamedTempFile {
    temp_with(
        "##gff-version 3\n\
         chr1\tx\tgene\t1000\t4000\t.\t+\t.\tID=G1;gene_name=MYGENE;gene_type=protein_coding\n\
         chr1\tx\tmRNA\t1000\t4000\t.\t+\t.\tID=NM_1;Parent=G1;transcript_id=NM_1\n\
         chr1\tx\texon\t1000\t2000\t.\t+\t.\tID=e1;Parent=NM_1\n\
         chr1\tx\tCDS\t1500\t2000\t.\t+\t0\tID=c1;Parent=NM_1\n\
         chr1\tx\texon\t3000\t4000\t.\t+\t.\tID=e2;Parent=NM_1\n\
         chr1\tx\tCDS\t3000\t3500\t.\t+\t0\tID=c2;Parent=NM_1\n",
    )
}

/// A SamBuilder with one `chr1` contig long enough for these tests. Coordinate-sorted, as `rna`
/// requires.
fn builder() -> SamBuilder {
    SamBuilder::with_contigs(&[("chr1".to_string(), 100_000)]).sort_order(SortOrder::Coordinate)
}

/// Render CIGAR ops to a string (e.g. `[(Match,50),(Skip,2000)]` → `"50M2000N"`).
fn cigar_string(ops: &[(Kind, usize)]) -> String {
    let mut s = String::new();
    for &(kind, len) in ops {
        let c = match kind {
            Kind::Skip => 'N',
            Kind::Deletion => 'D',
            Kind::Insertion => 'I',
            Kind::SoftClip => 'S',
            _ => 'M',
        };
        s.push_str(&len.to_string());
        s.push(c);
    }
    s
}

/// Build one mate of an RNA read pair, with the mate's CIGAR in the `MC` tag.
#[allow(clippy::too_many_arguments)]
fn rna_mate(
    name: &str,
    pos: usize,
    ops: &[(Kind, usize)],
    reverse: bool,
    first: bool,
    mate_pos: usize,
    mate_ops: &[(Kind, usize)],
    mate_reverse: bool,
    tlen: i32,
) -> RecordBuf {
    let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
    flags |= if first { Flags::FIRST_SEGMENT } else { Flags::LAST_SEGMENT };
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    if mate_reverse {
        flags |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    let read_len: usize = ops
        .iter()
        .filter(|(k, _)| matches!(k, Kind::Match | Kind::Insertion | Kind::SoftClip))
        .map(|(_, l)| *l)
        .sum();
    let cigar: Cigar = ops.iter().map(|&(k, l)| Op::new(k, l)).collect();

    let mut data = Data::default();
    data.insert((*b"MC").into(), DataValue::String(cigar_string(mate_ops).into()));

    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(pos).unwrap())
        .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(mate_pos).unwrap())
        .set_template_length(tlen)
        .set_sequence(Sequence::from(vec![b'A'; read_len]))
        .set_quality_scores(QualityScores::from(vec![30u8; read_len]))
        .set_data(data)
        .build()
}

/// Build one unpaired read with an explicit CIGAR (e.g. `50M999N50M` for a spliced read).
fn spliced_unpaired(name: &str, pos: usize, ops: &[(Kind, usize)]) -> RecordBuf {
    let read_len: usize = ops
        .iter()
        .filter(|(k, _)| matches!(k, Kind::Match | Kind::Insertion | Kind::SoftClip))
        .map(|(_, l)| *l)
        .sum();
    let cigar: Cigar = ops.iter().map(|&(k, l)| Op::new(k, l)).collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(pos).unwrap())
        .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; read_len]))
        .set_quality_scores(QualityScores::from(vec![30u8; read_len]))
        .build()
}

/// Build one read of a 50M FR pair with **no** `MC` tag, to exercise the TLEN-based insert-size
/// fallback. `tlen` is the signed genomic template length (+ on the leftmost read).
#[allow(clippy::fn_params_excessive_bools)]
fn fr_read_no_mc(
    name: &str,
    pos: usize,
    reverse: bool,
    first: bool,
    mate_pos: usize,
    mate_reverse: bool,
    tlen: i32,
) -> RecordBuf {
    let mut flags = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
    flags |= if first { Flags::FIRST_SEGMENT } else { Flags::LAST_SEGMENT };
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    if mate_reverse {
        flags |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(pos).unwrap())
        .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
        .set_cigar(cigar)
        .set_mate_reference_sequence_id(0)
        .set_mate_alignment_start(Position::new(mate_pos).unwrap())
        .set_template_length(tlen)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30u8; 50]))
        .build()
}

// ─── Base classification ────────────────────────────────────────────────────────

/// Place one 100bp read squarely in each functional region and check the base tallies.
fn assert_classification(gene_model: &NamedTempFile) {
    let mut b = builder();
    b.add_unpaired("coding", 0, 1500, 60, 100, false, false, false, None); // 1500..1599 coding
    b.add_unpaired("utr", 0, 1000, 60, 100, false, false, false, None); // 1000..1099 UTR
    b.add_unpaired("intron", 0, 2400, 60, 100, false, false, false, None); // 2400..2499 intron
    b.add_unpaired("inter", 0, 6000, 60, 100, false, false, false, None); // 6000..6099 intergenic
    let bam = b.to_temp_bam().unwrap();

    let (_dir, m) = run_rna(bam.path(), gene_model.path(), |_| {});
    assert_eq!(m.coding_bases, 100, "coding");
    assert_eq!(m.utr_bases, 100, "utr");
    assert_eq!(m.intronic_bases, 100, "intronic");
    assert_eq!(m.intergenic_bases, 100, "intergenic");
    assert_eq!(m.aligned_bases, 400);
    assert_eq!(m.bases, 400);
}

#[test]
fn classifies_bases_from_refflat() {
    assert_classification(&refflat_two_exon());
}

#[test]
fn rejects_non_coordinate_sorted_input() {
    // A builder left in the default (unsorted) state has no @HD SO:coordinate.
    let mut b = SamBuilder::with_contigs(&[("chr1".to_string(), 100_000)]);
    b.add_unpaired("r", 0, 1500, 60, 100, false, false, false, None);
    let bam = b.to_temp_bam().unwrap();
    let gm = gtf_two_exon();

    let dir = TempDir::new().unwrap();
    let cmd = Rna {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: dir.path().join("out") },
        reference: OptionalReferenceOptions { reference: None },
        options: RnaOptions { gene_model: gm.path().to_path_buf(), ..Default::default() },
    };
    let err = cmd.execute(None).expect_err("unsorted input must be rejected");
    assert!(err.to_string().contains("coordinate-sorted"), "unexpected error: {err}");
}

#[test]
fn classifies_bases_from_gtf() {
    assert_classification(&gtf_two_exon());
}

#[test]
fn classifies_bases_from_gff3() {
    assert_classification(&gff3_two_exon());
}

#[test]
fn all_three_formats_agree() {
    // Same reads, three formats → identical base composition (validates parsers + auto-detect).
    let make = |gm: &NamedTempFile| {
        let mut b = builder();
        b.add_unpaired("c", 0, 1600, 60, 100, false, false, false, None);
        b.add_unpaired("u", 0, 1000, 60, 100, false, false, false, None);
        let bam = b.to_temp_bam().unwrap();
        let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
        (m.coding_bases, m.utr_bases, m.intronic_bases, m.intergenic_bases)
    };
    let rf = make(&refflat_two_exon());
    let gtf = make(&gtf_two_exon());
    let gff = make(&gff3_two_exon());
    assert_eq!(rf, gtf, "refFlat vs GTF");
    assert_eq!(rf, gff, "refFlat vs GFF3");
}

// ─── Read accounting & genomic origin ──────────────────────────────────────────────

#[test]
fn read_origin_partitions_mapped_reads() {
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_unpaired("coding", 0, 1500, 60, 100, false, false, false, None); // exonic
    b.add_unpaired("intron", 0, 2400, 60, 100, false, false, false, None); // intronic
    b.add_unpaired("inter", 0, 6000, 60, 100, false, false, false, None); // intergenic
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.total_reads, 3);
    assert_eq!(m.mapped_reads, 3);
    assert_eq!(m.unmapped_reads, 0);
    assert_eq!(m.exonic_reads, 1);
    assert_eq!(m.intronic_reads, 1);
    assert_eq!(m.intergenic_reads, 1);
    assert_eq!(m.ambiguous_reads, 0);
    assert_eq!(m.assigned_reads, 1);
    // The four read-origin categories partition the mapped reads.
    assert_eq!(m.exonic_reads + m.intronic_reads + m.intergenic_reads, m.mapped_reads);
}

#[test]
fn unmapped_reads_are_counted_in_the_chain() {
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_unpaired("coding", 0, 1500, 60, 100, false, false, false, None);
    b.add_unmapped("u1", 100, false, None);
    b.add_unmapped("u2", 100, false, None);
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.total_reads, 3, "1 mapped + 2 unmapped primary reads");
    assert_eq!(m.mapped_reads, 1);
    assert_eq!(m.unmapped_reads, 2);
}

#[test]
fn genes_detected_respects_the_read_threshold() {
    let gm = gtf_two_exon();
    let mut b = builder();
    for i in 0..5 {
        b.add_unpaired(&format!("r{i}"), 0, 1500, 60, 100, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    // Default threshold is 5 → the single gene clears it.
    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.assigned_reads, 5);
    assert_eq!(m.genes_detected, 1);

    // Raising the threshold above the observed count drops it back to zero.
    let (_d2, m2) = run_rna(bam.path(), gm.path(), |o| o.genes_detected_min_reads = 6);
    assert_eq!(m2.genes_detected, 0);
}

#[test]
fn writes_biotype_and_gene_expression_plots() {
    let gm = gtf_two_exon();
    let mut b = builder();
    for i in 0..5 {
        b.add_unpaired(&format!("r{i}"), 0, 1500, 60, 100, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();
    let (dir, _m) = run_rna(bam.path(), gm.path(), |_| {});

    // Both charts are always produced.
    for suffix in [BIOTYPE_PLOT_SUFFIX, GENE_EXPRESSION_PLOT_SUFFIX] {
        let path = dir.path().join(format!("out{suffix}"));
        assert!(path.exists(), "expected output file {}", path.display());
    }
}

#[test]
fn reads_exonic_to_two_genes_are_ambiguous() {
    // Two protein-coding genes with exons overlapping over [1500,2000].
    let a =
        "gene_id \"A\"; transcript_id \"TA\"; gene_name \"GENEA\"; gene_type \"protein_coding\";";
    let c =
        "gene_id \"B\"; transcript_id \"TB\"; gene_name \"GENEB\"; gene_type \"protein_coding\";";
    let gm = temp_with(&format!(
        "chr1\tx\texon\t1000\t2000\t.\t+\t.\t{a}\n\
         chr1\tx\texon\t1500\t2500\t.\t+\t.\t{c}\n"
    ));

    let mut b = builder();
    b.add_unpaired("both", 0, 1600, 60, 100, false, false, false, None); // 1600..1699 in both exons
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.exonic_reads, 1);
    assert_eq!(m.ambiguous_reads, 1);
    assert_eq!(m.assigned_reads, 0);
    assert_eq!(m.genes_detected, 0);
}

// ─── Splice junctions ──────────────────────────────────────────────────────────────

#[test]
fn classifies_splice_junctions_known_partial_novel() {
    use Kind::{Match as M, Skip as N};
    // gtf_two_exon exons [1000,2000] and [3000,4000] → one annotated intron (2001,2999).
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_record(spliced_unpaired("known", 1951, &[(M, 50), (N, 999), (M, 50)])); // (2001,2999)
    b.add_record(spliced_unpaired("partial", 1951, &[(M, 50), (N, 1049), (M, 50)])); // donor known
    b.add_record(spliced_unpaired("novel", 1500, &[(M, 50), (N, 200), (M, 50)])); // both novel
    b.add_record(spliced_unpaired("shortgap", 1500, &[(M, 50), (N, 30), (M, 50)])); // < min-intron
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.spliced_reads, 3, "the sub-min-intron gap is not a splice");
    assert_eq!(m.known_splice_obs, 1);
    assert_eq!(m.partial_splice_obs, 1);
    assert_eq!(m.novel_splice_obs, 1);
    assert_eq!(m.known_juncs, 1);
    assert_eq!(m.partial_juncs, 1);
    assert_eq!(m.novel_juncs, 1);
    assert!((m.frac_known_juncs - 1.0 / 3.0).abs() < 1e-4); // serialized to 5 dp
}

#[test]
fn repeated_junction_counts_observations_but_one_distinct() {
    use Kind::{Match as M, Skip as N};
    let gm = gtf_two_exon();
    let mut b = builder();
    // Three reads spanning the same annotated junction.
    for i in 0..3 {
        b.add_record(spliced_unpaired(&format!("r{i}"), 1951, &[(M, 50), (N, 999), (M, 50)]));
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.known_splice_obs, 3, "observation-level counts multiplicity");
    assert_eq!(m.known_juncs, 1, "distinct-junction count collapses them");
    assert_eq!(m.spliced_reads, 3);
}

// ─── Transcript integrity (TIN) ──────────────────────────────────────────────────────

/// A single-exon gene over `[start, end]` (1-based inclusive) as GTF.
fn single_exon_gtf(start: usize, end: usize) -> NamedTempFile {
    let attrs =
        "gene_id \"G1\"; transcript_id \"T1\"; gene_name \"G\"; gene_type \"protein_coding\";";
    temp_with(&format!("chr1\tx\texon\t{start}\t{end}\t.\t+\t.\t{attrs}\n"))
}

#[test]
fn tin_is_high_for_uniform_coverage() {
    let gm = single_exon_gtf(1000, 2000); // 1001 bp
    let mut b = builder();
    for i in 0..12 {
        // Each read spans the whole exon → perfectly uniform coverage of depth 12.
        b.add_unpaired(&format!("r{i}"), 0, 1000, 60, 1001, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.tin_transcripts, 1);
    assert!(m.median_tin > 99.0, "uniform coverage → TIN≈100, got {}", m.median_tin);
}

#[test]
fn tin_is_low_for_deep_but_spiky_coverage() {
    let gm = single_exon_gtf(1000, 2000); // 1001 bp
    let mut b = builder();
    // 150 reads pile on the first 100 bp: mean coverage ~15 (clears the mean gate) but the profile
    // is highly uneven — exactly the deeply-but-unevenly-covered case the mean gate should keep.
    for i in 0..150 {
        b.add_unpaired(&format!("r{i}"), 0, 1000, 60, 100, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.tin_transcripts, 1, "deep enough to clear the mean-coverage gate");
    assert!(m.median_tin < 20.0, "spiky coverage → low TIN, got {}", m.median_tin);
}

#[test]
fn tin_excludes_low_mean_coverage_transcripts() {
    let gm = single_exon_gtf(1000, 2000);
    let mut b = builder();
    for i in 0..5 {
        // Mean coverage 5 < default --tin-min-coverage 10 → excluded.
        b.add_unpaired(&format!("r{i}"), 0, 1000, 60, 1001, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.tin_transcripts, 0);
    assert!(m.median_tin.abs() < 1e-9, "no contributing transcripts → TIN 0");
}

#[test]
fn tin_excludes_short_transcripts() {
    let gm = single_exon_gtf(1000, 1300); // 301 bp < --minimum-length (500)
    let mut b = builder();
    for i in 0..20 {
        // Deeply covered (mean ~20) but too short to estimate integrity reliably.
        b.add_unpaired(&format!("r{i}"), 0, 1000, 60, 301, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.tin_transcripts, 0, "transcript shorter than --minimum-length is excluded");
}

// ─── Biotype file ──────────────────────────────────────────────────────────────────

#[test]
fn biotype_file_groups_assigned_reads_by_biotype() {
    let pc =
        "gene_id \"G1\"; transcript_id \"T1\"; gene_name \"MYGENE\"; gene_type \"protein_coding\";";
    let lc = "gene_id \"G2\"; transcript_id \"T2\"; gene_name \"LINC\"; gene_type \"lincRNA\";";
    let gm = temp_with(&format!(
        "chr1\tx\texon\t1000\t2000\t.\t+\t.\t{pc}\n\
         chr1\tx\tCDS\t1500\t2000\t.\t+\t0\t{pc}\n\
         chr1\tx\texon\t5000\t6000\t.\t+\t.\t{lc}\n"
    ));

    let mut b = builder();
    for i in 0..3 {
        b.add_unpaired(&format!("pc{i}"), 0, 1500, 60, 100, false, false, false, None);
    }
    for i in 0..2 {
        b.add_unpaired(&format!("lc{i}"), 0, 5500, 60, 100, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();

    let (dir, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.assigned_reads, 5);
    let rows: Vec<RnaBiotypeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{BIOTYPE_SUFFIX}"))).unwrap();
    assert_eq!(rows.len(), 2);
    // Sorted by reads desc: protein_coding (3) before lincRNA (2).
    assert_eq!(rows[0].biotype, "protein_coding");
    assert_eq!(rows[0].reads, 3);
    assert!((rows[0].frac_reads - 0.6).abs() < 1e-4);
    assert_eq!(rows[1].biotype, "lincRNA");
    assert_eq!(rows[1].reads, 2);
}

#[test]
fn refflat_biotype_inferred_from_cds_presence() {
    // MYGENE has a CDS (coding); NCGENE has cdsStart == cdsEnd (non-coding).
    let gm = temp_with(
        "MYGENE\tNM_1\tchr1\t+\t999\t4000\t1499\t3500\t2\t999,2999,\t2000,4000,\n\
         NCGENE\tNR_1\tchr1\t+\t5000\t6000\t6000\t6000\t1\t5000,\t6000,\n",
    );
    let mut b = builder();
    for i in 0..3 {
        b.add_unpaired(&format!("c{i}"), 0, 1500, 60, 100, false, false, false, None); // coding
    }
    for i in 0..2 {
        b.add_unpaired(&format!("n{i}"), 0, 5500, 60, 100, false, false, false, None); // non-coding
    }
    let bam = b.to_temp_bam().unwrap();

    let (dir, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.assigned_reads, 5);
    let rows: Vec<RnaBiotypeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{BIOTYPE_SUFFIX}"))).unwrap();
    let pc = rows.iter().find(|r| r.biotype == "protein_coding").expect("protein_coding row");
    assert_eq!(pc.reads, 3);
    let nc = rows.iter().find(|r| r.biotype == "noncoding").expect("noncoding row");
    assert_eq!(nc.reads, 2);
}

// ─── Ribosomal ───────────────────────────────────────────────────────────────────

#[test]
fn ribosomal_derived_from_gtf_biotype() {
    // A single-exon rRNA gene over [1000,2000]; a read fully inside is counted ribosomal.
    let attrs = "gene_id \"R1\"; transcript_id \"RT1\"; gene_name \"RRNA\"; gene_type \"rRNA\";";
    let gm = temp_with(&format!("chr1\tx\texon\t1000\t2000\t.\t+\t.\t{attrs}\n"));

    let mut b = builder();
    b.add_unpaired("rrna", 0, 1200, 60, 100, false, false, false, None); // inside the rRNA gene
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.ribosomal_bases, Some(100), "all 100 bases ribosomal");
    // Ribosomal early-return means the read is NOT also counted as UTR/coding.
    assert_eq!(m.utr_bases, 0);
    assert_eq!(m.coding_bases, 0);
    assert!(m.frac_ribosomal_bases.is_some_and(|f| (f - 1.0).abs() < 1e-9));
}

#[test]
fn ribosomal_from_explicit_interval_list_with_refflat() {
    // refFlat carries no biotype, so ribosomal must come from the explicit interval list.
    let gm = refflat_two_exon();
    // Picard IntervalList: header @SQ then 1-based inclusive interval over the gene.
    let ribo = temp_with("@SQ\tSN:chr1\tLN:100000\nchr1\t1000\t2000\t+\trRNA\n");

    let mut b = builder();
    b.add_unpaired("ribo", 0, 1000, 60, 100, false, false, false, None);
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |o| {
        o.ribosomal_intervals = Some(ribo.path().to_path_buf());
    });
    assert_eq!(m.ribosomal_bases, Some(100));
}

#[test]
fn ribosomal_is_blank_without_a_source() {
    // refFlat (no biotype) and no explicit intervals → ribosomal metrics left blank.
    let gm = refflat_two_exon();
    let mut b = builder();
    b.add_unpaired("x", 0, 1600, 60, 100, false, false, false, None);
    let bam = b.to_temp_bam().unwrap();
    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.ribosomal_bases, None);
    assert_eq!(m.frac_ribosomal_bases, None);
}

// ─── Transcript-space insert size ───────────────────────────────────────────────

#[test]
fn insert_size_collapses_introns() {
    // Gene with exons [1000,2000] and [4000,5000] (intron 2001..3999).
    let attrs =
        "gene_id \"G\"; transcript_id \"T\"; gene_name \"G\"; gene_type \"protein_coding\";";
    let gm = temp_with(&format!(
        "chr1\tx\texon\t1000\t2000\t.\t+\t.\t{attrs}\n\
         chr1\tx\texon\t4000\t5000\t.\t+\t.\t{attrs}\n"
    ));

    // R1 forward 50M at 1500 (5' = 1500); R2 reverse 50M at 4451 (5' = end = 4500).
    let r1 = rna_mate(
        "p",
        1500,
        &[(Kind::Match, 50)],
        false,
        true,
        4451,
        &[(Kind::Match, 50)],
        true,
        3050,
    );
    let r2 = rna_mate(
        "p",
        4451,
        &[(Kind::Match, 50)],
        true,
        false,
        1500,
        &[(Kind::Match, 50)],
        false,
        -3050,
    );
    let mut b = builder();
    b.add_record(r1);
    b.add_record(r2);
    let bam = b.to_temp_bam().unwrap();

    let dir = TempDir::new().unwrap();
    let prefix = dir.path().join("out");
    Rna {
        input: InputOptions { input: bam.path().to_path_buf() },
        output: OutputOptions { output: prefix },
        reference: OptionalReferenceOptions { reference: None },
        options: RnaOptions { gene_model: gm.path().to_path_buf(), ..Default::default() },
    }
    .execute(None)
    .unwrap();

    let isize: Vec<RnaInsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{ISIZE_SUFFIX}"))).unwrap();
    let fr = isize.iter().find(|m| m.pair_orientation == "FR").unwrap();
    assert_eq!(fr.read_pairs, 1, "one FR pair counted");
    // 5'-to-5' exonic distance = (1500..2000) + (4000..4500) = 501 + 501 = 1002.
    assert!(
        (fr.median_insert_size - 1002.0).abs() < 1e-9,
        "transcript-space insert size should collapse the intron, got {}",
        fr.median_insert_size
    );
}

fn two_exon_isize_model() -> NamedTempFile {
    // Gene with exons [1000,2000] and [4000,5000] (intron 2001..3999).
    let attrs =
        "gene_id \"G\"; transcript_id \"T\"; gene_name \"G\"; gene_type \"protein_coding\";";
    temp_with(&format!(
        "chr1\tx\texon\t1000\t2000\t.\t+\t.\t{attrs}\n\
         chr1\tx\texon\t4000\t5000\t.\t+\t.\t{attrs}\n"
    ))
}

#[test]
fn insert_size_estimated_from_tlen_when_mc_absent() {
    let gm = two_exon_isize_model();
    // Forward read [1500,1549] in exon 1; reverse mate [4451,4500] in exon 2; genomic span
    // 1500..4500 → |TLEN| 3001. With no MC tag the FR fallback derives the mate end from TLEN and
    // must recover the same transcript-space insert size as the MC path.
    let mut b = builder();
    b.add_record(fr_read_no_mc("p", 1500, false, true, 4451, true, 3001));
    b.add_record(fr_read_no_mc("p", 4451, true, false, 1500, false, -3001));
    let bam = b.to_temp_bam().unwrap();

    let (dir, _m) = run_rna(bam.path(), gm.path(), |_| {});
    let isize: Vec<RnaInsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{ISIZE_SUFFIX}"))).unwrap();
    let fr = isize.iter().find(|m| m.pair_orientation == "FR").unwrap();
    assert_eq!(fr.read_pairs, 1, "FR pair recovered from TLEN");
    // Identical to the MC path: (1500..2000) + (4000..4500) = 501 + 501 = 1002.
    assert!(
        (fr.median_insert_size - 1002.0).abs() < 1e-9,
        "TLEN fallback should match the MC result, got {}",
        fr.median_insert_size
    );
}

#[test]
fn tlen_fallback_rejects_a_mate_outside_the_exons() {
    let gm = two_exon_isize_model();
    // Reverse mate [2500,2549] sits in the intron. Without MC we can't see its splicing, so the
    // both-ends-in-an-exon guard must reject the pair rather than guess.
    let mut b = builder();
    b.add_record(fr_read_no_mc("p", 1500, false, true, 2500, true, 1050));
    b.add_record(fr_read_no_mc("p", 2500, true, false, 1500, false, -1050));
    let bam = b.to_temp_bam().unwrap();

    let (dir, _m) = run_rna(bam.path(), gm.path(), |_| {});
    let isize: Vec<RnaInsertSizeMetric> =
        read_metrics_tsv(&dir.path().join(format!("out{ISIZE_SUFFIX}"))).unwrap();
    let fr = isize.iter().find(|m| m.pair_orientation == "FR").unwrap();
    assert_eq!(fr.read_pairs, 0, "an intronic mate must not be counted");
}

// ─── Strand auto-detection ───────────────────────────────────────────────────────

#[test]
fn auto_detects_forward_strandedness() {
    // FR pairs enclosed in the + gene with R1 forward → read-1 on transcription strand.
    let gm = gtf_two_exon();
    let mut b = builder();
    for i in 0..10 {
        // R1 forward 50M in exon1 (1500); R2 reverse 50M in exon2 (3000), enclosed in [1000,4000].
        let name = format!("p{i}");
        let r1 = rna_mate(
            &name,
            1500,
            &[(Kind::Match, 50)],
            false,
            true,
            3000,
            &[(Kind::Match, 50)],
            true,
            1550,
        );
        let r2 = rna_mate(
            &name,
            3000,
            &[(Kind::Match, 50)],
            true,
            false,
            1500,
            &[(Kind::Match, 50)],
            false,
            -1550,
        );
        b.add_record(r1);
        b.add_record(r2);
    }
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.r1_tx_strand_reads, 10, "expected one template count per FR pair");
    assert_eq!(m.r2_tx_strand_reads, 0);
    assert_eq!(m.detected_strand, "forward");
    // Under forward strandedness, R1 reads agreeing with the + gene are correct.
    assert!(m.correct_strand_reads > 0);
}

#[test]
fn explicit_strand_override_is_respected() {
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_unpaired("c", 0, 1600, 60, 100, false, false, false, None);
    let bam = b.to_temp_bam().unwrap();
    let (_d, m) = run_rna(bam.path(), gm.path(), |o| o.strand = StrandSpec::Reverse);
    assert_eq!(m.detected_strand, "reverse");
}

// ─── Coverage / bias (bug-fix behaviours) ────────────────────────────────────────

#[test]
fn uniform_coverage_gives_flat_bias_and_zero_cv() {
    // Single-exon non-coding gene [1,600] (≥ minimum_length and ≥ 2×end_bias_bases).
    let gm = temp_with("UNI\tT\tchr1\t+\t0\t600\t600\t600\t1\t0,\t600,\n");
    let mut b = builder();
    // Five reads each fully spanning the transcript → perfectly uniform depth.
    for i in 0..5 {
        b.add_unpaired(&format!("r{i}"), 0, 1, 60, 600, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();
    let (dir, m) = run_rna(bam.path(), gm.path(), |_| {});

    assert!(
        m.median_cv_coverage.abs() < 1e-9,
        "uniform coverage → CV 0, got {}",
        m.median_cv_coverage
    );
    assert!((m.median_5prime_bias - 1.0).abs() < 1e-9);
    assert!((m.median_3prime_bias - 1.0).abs() < 1e-9);
    assert!((m.median_5prime_to_3prime_bias - 1.0).abs() < 1e-9);

    // The coverage chart is produced (the per-position table is no longer emitted).
    assert!(dir.path().join(format!("out{COVERAGE_PLOT_SUFFIX}")).exists());
}

#[test]
fn short_transcripts_are_excluded_from_bias() {
    // A 300bp transcript is shorter than 2×end_bias_bases (200) is fine, but shorter than the
    // 500bp default minimum_length, so it is excluded → no coverage metrics.
    let gm = temp_with("SHORT\tT\tchr1\t+\t0\t300\t300\t300\t1\t0,\t300,\n");
    let mut b = builder();
    for i in 0..5 {
        b.add_unpaired(&format!("r{i}"), 0, 1, 60, 300, false, false, false, None);
    }
    let bam = b.to_temp_bam().unwrap();
    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert!(m.median_cv_coverage.abs() < 1e-9);
    assert!(m.median_5prime_to_3prime_bias.abs() < 1e-9, "short transcript excluded → 0");
}

// ─── Filtering ────────────────────────────────────────────────────────────────────

#[test]
fn duplicates_are_included_by_default_and_excludable() {
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_unpaired("keep", 0, 1600, 60, 100, false, false, false, None);
    b.add_unpaired("dup", 0, 1600, 60, 100, false, true, false, None); // duplicate-flagged
    let bam = b.to_temp_bam().unwrap();

    // Default: duplicates counted → 200 coding bases, duplicate_rate 0.5.
    let (_d, included) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(included.coding_bases, 200);
    assert!((included.duplicate_rate - 0.5).abs() < 1e-9);

    // With --exclude-duplicates the dup read is dropped from the metrics.
    let (_d, excluded) = run_rna(bam.path(), gm.path(), |o| o.exclude_duplicates = true);
    assert_eq!(excluded.coding_bases, 100);
}

#[test]
fn ignore_chrom_counts_only_as_ignored_reads() {
    let gm = gtf_two_exon();
    let mut b =
        SamBuilder::with_contigs(&[("chr1".to_string(), 100_000), ("chrM".to_string(), 16_000)])
            .sort_order(SortOrder::Coordinate);
    b.add_unpaired("gene", 0, 1600, 60, 100, false, false, false, None);
    b.add_unpaired("mito", 1, 1000, 60, 100, false, false, false, None); // on chrM
    let bam = b.to_temp_bam().unwrap();

    let (_d, m) = run_rna(bam.path(), gm.path(), |o| {
        o.ignore_chrom = Some(vec!["chrM".to_string()]);
    });
    assert_eq!(m.ignored_reads, 1);
    assert_eq!(m.coding_bases, 100, "only the gene read contributes aligned bases");
    assert_eq!(m.bases, 200, "ignored read still counts in total bases");
}

// ─── QC-fail / unmapped ────────────────────────────────────────────────────────────

#[test]
fn qc_fail_reads_are_dropped_entirely() {
    let gm = gtf_two_exon();
    let mut b = builder();
    b.add_unpaired("ok", 0, 1600, 60, 100, false, false, false, None);
    b.add_unpaired("qc", 0, 1600, 60, 100, false, false, true, None); // QC fail
    let bam = b.to_temp_bam().unwrap();
    let (_d, m) = run_rna(bam.path(), gm.path(), |_| {});
    assert_eq!(m.bases, 100, "QC-fail read excluded from bases");
    assert_eq!(m.coding_bases, 100);
}

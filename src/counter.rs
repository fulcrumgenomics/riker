use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{AddAssign, Index};

use num_traits::{PrimInt, ToPrimitive};

/// A frequency counter backed by `HashMap<K, u64>`.
///
/// Methods are gated by key type via separate `impl` blocks with progressively
/// tighter trait bounds:
///
/// - `Eq + Hash + Copy` — basic counting, iteration, merge
/// - `+ Ord` — sorted output, mode, min, max
/// - `+ ToPrimitive` — mean, stddev, trimmed stats
/// - `PrimInt + ToPrimitive + Hash` — median, MAD (integer keys only)
#[derive(Clone, Debug)]
pub struct Counter<K> {
    map: HashMap<K, u64>,
    total_count: u64,
}

// ─── impl block 1: K: Eq + Hash + Copy — basic counting ─────────────────────

impl<K: Eq + Hash + Copy> Counter<K> {
    /// Create an empty counter.
    #[must_use]
    pub fn new() -> Self {
        Self { map: HashMap::new(), total_count: 0 }
    }

    /// Increment the count for `key` by 1.
    pub fn count(&mut self, key: K) {
        *self.map.entry(key).or_insert(0) += 1;
        self.total_count += 1;
    }

    /// Increment the count for `key` by `n`.
    pub fn count_n(&mut self, key: K, n: u64) {
        *self.map.entry(key).or_insert(0) += n;
        self.total_count += n;
    }

    /// Return the count for `key`, or 0 if absent.
    #[must_use]
    pub fn count_of(&self, key: &K) -> u64 {
        self.map.get(key).copied().unwrap_or(0)
    }

    /// Return the sum of all counts.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.total_count
    }

    /// Return the number of distinct keys.
    #[must_use]
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Return true if no keys have been counted.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Merge all counts from `other` into `self`.
    pub fn merge(&mut self, other: &Counter<K>) {
        for (&key, &count) in &other.map {
            *self.map.entry(key).or_insert(0) += count;
        }
        self.total_count += other.total_count;
    }

    /// Iterate over `(&K, &u64)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &u64)> {
        self.map.iter()
    }

    /// Iterate over keys.
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.map.keys()
    }

    /// Iterate over counts.
    pub fn values(&self) -> impl Iterator<Item = &u64> {
        self.map.values()
    }
}

impl<K: Eq + Hash + Copy> Default for Counter<K> {
    fn default() -> Self {
        Self::new()
    }
}

// ─── impl block 2: K: Eq + Hash + Copy + Ord — ordered operations ───────────

impl<K: Eq + Hash + Copy + Ord> Counter<K> {
    /// Return entries sorted by key.
    #[must_use]
    pub fn sorted(&self) -> Vec<(K, u64)> {
        let mut entries: Vec<(K, u64)> = self.map.iter().map(|(&k, &v)| (k, v)).collect();
        entries.sort_unstable_by_key(|&(k, _)| k);
        entries
    }

    /// Return the most frequent key (mode).  On tie, the largest key wins.
    #[must_use]
    pub fn mode(&self) -> Option<K> {
        self.map.iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| a.0.cmp(b.0))).map(|(&k, _)| k)
    }

    /// Return the smallest key with a nonzero count.
    #[must_use]
    pub fn min(&self) -> Option<K> {
        self.map.keys().min().copied()
    }

    /// Return the largest key with a nonzero count.
    #[must_use]
    pub fn max(&self) -> Option<K> {
        self.map.keys().max().copied()
    }
}

// ─── impl block 3: K: Eq + Hash + Copy + Ord + ToPrimitive — numeric stats ──

impl<K: Eq + Hash + Copy + Ord + ToPrimitive> Counter<K> {
    /// Compute the population mean.
    #[must_use]
    pub fn mean(&self) -> f64 {
        if self.total_count == 0 {
            return 0.0;
        }
        let sum: f64 = self.map.iter().map(|(k, &c)| k.to_f64().unwrap_or(0.0) * c as f64).sum();
        sum / self.total_count as f64
    }

    /// Compute the population standard deviation.
    #[must_use]
    pub fn stddev(&self) -> f64 {
        self.mean_and_stddev().1
    }

    /// Compute mean and population standard deviation in a single pass.
    #[must_use]
    pub fn mean_and_stddev(&self) -> (f64, f64) {
        self.trimmed_mean_and_stddev(f64::MAX)
    }

    /// Compute mean and stddev including only keys where `key.to_f64() <= max_val`.
    #[must_use]
    pub fn trimmed_mean_and_stddev(&self, max_val: f64) -> (f64, f64) {
        let mut n = 0u64;
        let mut sum = 0.0_f64;
        let mut sum_sq = 0.0_f64;

        for &(k, count) in &self.sorted() {
            let v = k.to_f64().unwrap_or(0.0);
            if v > max_val {
                break;
            }
            n += count;
            sum += v * count as f64;
            sum_sq += v * v * count as f64;
        }

        if n == 0 {
            return (0.0, 0.0);
        }
        let mean = sum / n as f64;
        let variance = (sum_sq / n as f64) - mean * mean;
        (mean, variance.max(0.0).sqrt())
    }
}

// ─── impl block 4: PrimInt + ToPrimitive + Hash — integer-only stats ─────────

impl<K: PrimInt + ToPrimitive + Hash> Counter<K> {
    /// Compute the weighted median (interpolated for even totals).
    /// Returns 0.0 for an empty counter.
    #[must_use]
    pub fn median(&self) -> f64 {
        let total = self.total_count;
        if total == 0 {
            return 0.0;
        }

        let sorted = self.sorted();

        // 0-indexed positions of the two "middle" elements.
        let lower_pos = (total - 1) / 2;
        let upper_pos = total / 2;

        let mut running = 0u64;
        let mut lower_val: Option<f64> = None;

        for &(val, count) in &sorted {
            let prev = running;
            running += count;
            let fval = val.to_f64().unwrap_or(0.0);

            if lower_val.is_none() && prev <= lower_pos && lower_pos < running {
                lower_val = Some(fval);
            }
            if prev <= upper_pos && upper_pos < running {
                return f64::midpoint(lower_val.unwrap_or(fval), fval);
            }
        }

        lower_val.unwrap_or(0.0)
    }

    /// Compute the median absolute deviation (MAD).
    /// Returns 0.0 for an empty counter.
    #[must_use]
    pub fn mad(&self) -> f64 {
        self.median_and_mad().1
    }

    /// Compute the median and MAD in a single pass over the sorted entries.
    #[must_use]
    pub fn median_and_mad(&self) -> (f64, f64) {
        let median = self.median();
        // Store deviations×2 as u32 keys to avoid floating-point map keys.
        let mut dev_hist = Counter::<u32>::new();
        for (&val, &count) in &self.map {
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "deviation×2 is non-negative and fits in u32 for realistic data"
            )]
            let dev2 = ((val.to_f64().unwrap_or(0.0) - median).abs() * 2.0).round() as u32;
            dev_hist.count_n(dev2, count);
        }
        (median, dev_hist.median() / 2.0)
    }
}

// ─── Index trait ─────────────────────────────────────────────────────────────

/// A static zero value for the `Index` implementation.
static ZERO: u64 = 0;

impl<K: Eq + Hash + Copy> Index<&K> for Counter<K> {
    type Output = u64;

    fn index(&self, key: &K) -> &u64 {
        self.map.get(key).unwrap_or(&ZERO)
    }
}

// ─── AddAssign trait ─────────────────────────────────────────────────────────

impl<K: Eq + Hash + Copy> AddAssign<&Counter<K>> for Counter<K> {
    fn add_assign(&mut self, rhs: &Counter<K>) {
        self.merge(rhs);
    }
}

// ─── IntoIterator traits ─────────────────────────────────────────────────────

impl<K: Eq + Hash + Copy> IntoIterator for Counter<K> {
    type Item = (K, u64);
    type IntoIter = std::collections::hash_map::IntoIter<K, u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.map.into_iter()
    }
}

impl<'a, K: Eq + Hash + Copy> IntoIterator for &'a Counter<K> {
    type Item = (&'a K, &'a u64);
    type IntoIter = std::collections::hash_map::Iter<'a, K, u64>;

    fn into_iter(self) -> Self::IntoIter {
        self.map.iter()
    }
}

// ─── FromIterator trait ──────────────────────────────────────────────────────

impl<K: Eq + Hash + Copy> FromIterator<(K, u64)> for Counter<K> {
    fn from_iter<I: IntoIterator<Item = (K, u64)>>(iter: I) -> Self {
        let mut counter = Counter::new();
        for (key, count) in iter {
            counter.count_n(key, count);
        }
        counter
    }
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn counter(pairs: &[(u32, u64)]) -> Counter<u32> {
        pairs.iter().copied().collect()
    }

    // ── Basic counting ──────────────────────────────────────────────────────

    #[test]
    fn test_new_is_empty() {
        let c = Counter::<u32>::new();
        assert!(c.is_empty());
        assert_eq!(c.len(), 0);
        assert_eq!(c.total(), 0);
    }

    #[test]
    fn test_count_single() {
        let mut c = Counter::new();
        c.count(42u32);
        assert_eq!(c.count_of(&42), 1);
        assert_eq!(c.total(), 1);
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn test_count_multiple() {
        let mut c = Counter::new();
        c.count(10u32);
        c.count(10);
        c.count(20);
        assert_eq!(c.count_of(&10), 2);
        assert_eq!(c.count_of(&20), 1);
        assert_eq!(c.total(), 3);
        assert_eq!(c.len(), 2);
    }

    #[test]
    fn test_count_n() {
        let mut c = Counter::new();
        c.count_n(100u32, 50);
        assert_eq!(c.count_of(&100), 50);
        assert_eq!(c.total(), 50);
    }

    #[test]
    fn test_count_of_absent() {
        let c = Counter::<u32>::new();
        assert_eq!(c.count_of(&999), 0);
    }

    #[test]
    fn test_default() {
        let c = Counter::<u32>::default();
        assert!(c.is_empty());
    }

    // ── Merge ───────────────────────────────────────────────────────────────

    #[test]
    fn test_merge() {
        let mut a = counter(&[(10, 2), (20, 3)]);
        let b = counter(&[(20, 7), (30, 1)]);
        a.merge(&b);
        assert_eq!(a.count_of(&10), 2);
        assert_eq!(a.count_of(&20), 10);
        assert_eq!(a.count_of(&30), 1);
        assert_eq!(a.total(), 13);
    }

    #[test]
    fn test_add_assign() {
        let mut a = counter(&[(10, 1)]);
        let b = counter(&[(10, 2), (20, 3)]);
        a += &b;
        assert_eq!(a.count_of(&10), 3);
        assert_eq!(a.count_of(&20), 3);
        assert_eq!(a.total(), 6);
    }

    // ── Iteration ───────────────────────────────────────────────────────────

    #[test]
    fn test_iter() {
        let c = counter(&[(10, 1), (20, 2)]);
        let mut items: Vec<_> = c.iter().map(|(&k, &v)| (k, v)).collect();
        items.sort_unstable();
        assert_eq!(items, vec![(10, 1), (20, 2)]);
    }

    #[test]
    fn test_keys() {
        let c = counter(&[(5, 1), (10, 2)]);
        let mut keys: Vec<_> = c.keys().copied().collect();
        keys.sort_unstable();
        assert_eq!(keys, vec![5, 10]);
    }

    #[test]
    fn test_values() {
        let c = counter(&[(5, 1), (10, 2)]);
        let mut vals: Vec<_> = c.values().copied().collect();
        vals.sort_unstable();
        assert_eq!(vals, vec![1, 2]);
    }

    #[test]
    fn test_into_iter_owned() {
        let c = counter(&[(5, 3)]);
        let items: Vec<_> = c.into_iter().collect();
        assert_eq!(items.len(), 1);
        assert_eq!(items[0], (5, 3));
    }

    #[test]
    fn test_into_iter_borrowed() {
        let c = counter(&[(5, 3)]);
        let items: Vec<_> = (&c).into_iter().collect();
        assert_eq!(items.len(), 1);
    }

    #[test]
    fn test_from_iterator() {
        let c: Counter<u32> = vec![(1u32, 10u64), (2, 20)].into_iter().collect();
        assert_eq!(c.count_of(&1), 10);
        assert_eq!(c.count_of(&2), 20);
        assert_eq!(c.total(), 30);
    }

    // ── Index trait ─────────────────────────────────────────────────────────

    #[test]
    fn test_index_present() {
        let c = counter(&[(42, 7)]);
        assert_eq!(c[&42], 7);
    }

    #[test]
    fn test_index_absent() {
        let c = Counter::<u32>::new();
        assert_eq!(c[&999], 0);
    }

    // ── Sorted, mode, min, max ──────────────────────────────────────────────

    #[test]
    fn test_sorted() {
        let c = counter(&[(30, 1), (10, 2), (20, 3)]);
        assert_eq!(c.sorted(), vec![(10, 2), (20, 3), (30, 1)]);
    }

    #[test]
    fn test_mode_single_peak() {
        let c = counter(&[(90, 1), (100, 3), (110, 2)]);
        assert_eq!(c.mode(), Some(100));
    }

    #[test]
    fn test_mode_tie_returns_largest() {
        let c = counter(&[(90, 2), (100, 2)]);
        assert_eq!(c.mode(), Some(100));
    }

    #[test]
    fn test_mode_empty() {
        let c = Counter::<u32>::new();
        assert_eq!(c.mode(), None);
    }

    #[test]
    fn test_min_max() {
        let c = counter(&[(50, 1), (100, 2), (200, 3)]);
        assert_eq!(c.min(), Some(50));
        assert_eq!(c.max(), Some(200));
    }

    #[test]
    fn test_min_max_empty() {
        let c = Counter::<u32>::new();
        assert_eq!(c.min(), None);
        assert_eq!(c.max(), None);
    }

    // ── Mean and stddev ─────────────────────────────────────────────────────

    #[test]
    fn test_mean_single() {
        let c = counter(&[(100, 5)]);
        assert_eq!(c.mean(), 100.0);
    }

    #[test]
    fn test_mean_weighted() {
        let c = counter(&[(100, 3), (200, 1)]);
        // (100*3 + 200*1) / 4 = 500/4 = 125
        assert_eq!(c.mean(), 125.0);
    }

    #[test]
    fn test_mean_empty() {
        let c = Counter::<u32>::new();
        assert_eq!(c.mean(), 0.0);
    }

    #[test]
    fn test_stddev_constant() {
        let c = counter(&[(100, 5)]);
        assert_eq!(c.stddev(), 0.0);
    }

    #[test]
    fn test_mean_and_stddev_symmetric() {
        let c = counter(&[(99, 1), (100, 1), (101, 1)]);
        let (mean, sd) = c.mean_and_stddev();
        assert_eq!(mean, 100.0);
        let expected_sd = (2.0_f64 / 3.0).sqrt();
        assert!((sd - expected_sd).abs() < 1e-9, "sd={sd} expected={expected_sd}");
    }

    #[test]
    fn test_mean_and_stddev_empty() {
        let c = Counter::<u32>::new();
        assert_eq!(c.mean_and_stddev(), (0.0, 0.0));
    }

    // ── Trimmed mean and stddev ─────────────────────────────────────────────

    #[test]
    fn test_trimmed_mean_constant() {
        let c = counter(&[(100, 5)]);
        let (mean, sd) = c.trimmed_mean_and_stddev(200.0);
        assert_eq!(mean, 100.0);
        assert_eq!(sd, 0.0);
    }

    #[test]
    fn test_trimmed_excludes_outliers() {
        let c = counter(&[(100, 1), (200, 1), (1000, 1)]);
        let (mean, _sd) = c.trimmed_mean_and_stddev(500.0);
        assert_eq!(mean, 150.0);
    }

    #[test]
    fn test_trimmed_mean_empty() {
        let c = Counter::<u32>::new();
        assert_eq!(c.trimmed_mean_and_stddev(100.0), (0.0, 0.0));
    }

    #[test]
    fn test_trimmed_all_excluded() {
        let c = counter(&[(1000, 5)]);
        assert_eq!(c.trimmed_mean_and_stddev(100.0), (0.0, 0.0));
    }

    // ── Median ──────────────────────────────────────────────────────────────

    #[test]
    fn test_median_single_value() {
        assert_eq!(counter(&[(100, 1)]).median(), 100.0);
    }

    #[test]
    fn test_median_odd_uniform() {
        let c = counter(&[(90, 1), (95, 1), (100, 1), (105, 1), (110, 1)]);
        assert_eq!(c.median(), 100.0);
    }

    #[test]
    fn test_median_even_same_bin() {
        assert_eq!(counter(&[(100, 4)]).median(), 100.0);
    }

    #[test]
    fn test_median_even_adjacent_bins() {
        assert_eq!(counter(&[(100, 1), (101, 1)]).median(), 100.5);
    }

    #[test]
    fn test_median_with_counts() {
        let c = counter(&[(100, 2), (101, 1), (102, 1), (103, 1)]);
        assert_eq!(c.median(), 101.0);
    }

    #[test]
    fn test_median_empty() {
        assert_eq!(Counter::<u32>::new().median(), 0.0);
    }

    #[test]
    fn test_median_heavily_weighted() {
        let c = counter(&[(10, 1), (20, 1000)]);
        assert_eq!(c.median(), 20.0);
    }

    #[test]
    fn test_median_large_counts() {
        let half = u64::MAX / 2;
        let c = counter(&[(100, half), (200, half)]);
        assert_eq!(c.median(), 150.0);
    }

    // ── MAD ─────────────────────────────────────────────────────────────────

    #[test]
    fn test_mad_constant() {
        assert_eq!(counter(&[(100, 5)]).mad(), 0.0);
    }

    #[test]
    fn test_mad_symmetric() {
        let c = counter(&[(98, 1), (99, 1), (100, 1), (101, 1), (102, 1)]);
        assert_eq!(c.mad(), 1.0);
    }

    #[test]
    fn test_mad_empty() {
        assert_eq!(Counter::<u32>::new().mad(), 0.0);
    }

    #[test]
    fn test_mad_single_element() {
        assert_eq!(counter(&[(50, 1)]).mad(), 0.0);
    }

    // ── Generic type tests ──────────────────────────────────────────────────

    #[test]
    fn test_u8_counter() {
        let c: Counter<u8> = [(1u8, 2u64), (2, 3)].into_iter().collect();
        assert_eq!(c.median(), 2.0);
        assert_eq!(c.mean(), 1.6);
    }

    #[test]
    fn test_u16_counter() {
        let c: Counter<u16> = [(100u16, 1u64), (200, 1)].into_iter().collect();
        assert_eq!(c.median(), 150.0);
    }

    #[test]
    fn test_u64_counter() {
        let c: Counter<u64> = [(1000u64, 3u64)].into_iter().collect();
        assert_eq!(c.median(), 1000.0);
        assert_eq!(c.mean(), 1000.0);
    }

    #[test]
    fn test_i32_counter() {
        let c: Counter<i32> = [(-10i32, 1u64), (0, 1), (10, 1)].into_iter().collect();
        assert_eq!(c.median(), 0.0);
        assert_eq!(c.mean(), 0.0);
    }

    // ── String keys (Ord, but not numeric) ──────────────────────────────────

    #[test]
    fn test_string_counter_sorted() {
        let mut c = Counter::<&str>::new();
        c.count("banana");
        c.count("apple");
        c.count("banana");
        assert_eq!(c.sorted(), vec![("apple", 1), ("banana", 2)]);
        assert_eq!(c.mode(), Some("banana"));
    }
}

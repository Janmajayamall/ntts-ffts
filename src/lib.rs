use rustfft::{
    num_complex::{Complex, ComplexFloat},
    num_traits::Pow,
    FftPlanner,
};

/// School book negacylic mult in ring Z[X]/(X^N+1)
fn school_book_negacylic_multiplication(p1: &[f64], p2: &[f64]) -> Vec<f64> {
    let n = p1.len();
    let mut res = vec![0f64; n];
    for i in 0..n {
        for j in 0..i + 1 {
            res[i] += p1[j] * p2[i - j];
        }

        for j in (i + 1)..n {
            res[i] -= p1[j] * p2[n - (j - i)];
        }
    }
    res
}

/// School book negacylic mul in ring Z_{q}[X]/X^N+1
fn school_book_negacylic_multiplication_mod_z_q(p1: &[u64], p2: &[u64], q: u64) -> Vec<u64> {
    let n = p1.len();
    let mut res = vec![0; n];
    for i in 0..n {
        for j in 0..i + 1 {
            res[i] = (res[i] + ((p1[j] as u128) * (p2[i - j] as u128) % (q as u128)) as u64) % q;
        }

        for j in (i + 1)..n {
            res[i] = (res[i] + q
                - ((p1[j] as u128) * (p2[n - (j - i)] as u128) % (q as u128)) as u64)
                % q;
        }
    }
    res
}

/// School book negacylic mul in ring Z_{q=u64}[X]/X^N+1
fn school_book_negacylic_multiplication_native_mod(p1: &[u64], p2: &[u64]) -> Vec<u64> {
    let n = p1.len();
    let mut res = vec![0u64; n];
    for i in 0..n {
        for j in 0..i + 1 {
            res[i] = res[i].wrapping_add(p1[j].wrapping_mul(p2[i - j]));
        }

        for j in (i + 1)..n {
            res[i] = res[i].wrapping_sub(p1[j].wrapping_mul(p2[n - (j - i)]));
        }
    }
    res
}

fn normalize_fft_output(p: &mut [Complex<f64>]) {
    let n = p.len() as f64;
    p.iter_mut().for_each(|v| *v = *v / n);
}

/// Multiplies (p1 x p2) mod X^N+1 using FFNT (https://eprint.iacr.org/2021/480.pdf)
fn negacylic_mul_with_twist(p1: &[f64], p2: &[f64]) -> Vec<f64> {
    /// Maps $p \in Z[X]/X^N+1$ to element $\in Z[X]/X^N-1$. One cannot map $p$ directly
    /// from big ring to small ring. Instead observe that X^N+1 factors into X^{N/2}-i and X^{N/2}+i.
    /// Ring $Z[X]/X^{N/2}-i$ maps to Z[X]/X^{N/2}-1 by replacing X with $\omega_{4(N/2)}X$
    /// (refer to accompaniying notes to understand why).
    ///
    /// Thus to map $p$ to element in smaller ring, first map it to X^N+1 by consdering upper half coefficients as imaginery counterparts of lower half coefficients. Then element-wise multiply resultling coefficient vector with vector of powers of $\omega_{4(N/2)}$.
    fn _map_from_big_ring_to_small_ring(p: &[f64]) -> Vec<Complex<f64>> {
        let n = p.len();
        let middle = p
            .iter()
            .take(p.len() / 2)
            .zip(p.iter().skip(p.len() / 2))
            .map(|(r, i)| Complex { re: *r, im: *i })
            .collect::<Vec<Complex<f64>>>();

        let powers_of_unity = _powers_of_root_of_unity(n);

        middle
            .iter()
            .zip(powers_of_unity.iter())
            .map(|(v, omega_power)| v * omega_power)
            .collect::<Vec<Complex<f64>>>()
    }

    /// Maps $p \in Z[X]/X^N-1$ to element $\in Z[X]/X^N+1$.
    /// This calculates the inverse of `_map_from_big_ring_to_small_ring`. It first multiplies $p$ with inverse of powers of $\omega_{4(N/2)}$. Then to map
    /// resulting element in $Z[X]/X^N-i$ to element in $Z[X]/X^N+1$, it takes reals in coefficient vector and places them in lower half of new coefficient vector.
    /// Then it considers imaginary counterparts as reals and appends them to the new coefficient vector.
    fn _map_from_small_ring_to_big_ring(p: &[Complex<f64>]) -> Vec<f64> {
        let n = p.len() * 2;

        let power_of_unity_inverse = _powers_of_root_of_unity(n)
            .iter()
            .map(|c| c.inv())
            .collect::<Vec<Complex<f64>>>();

        let middle = p
            .iter()
            .zip(power_of_unity_inverse.iter())
            .map(|(v, omega_power)| v * omega_power)
            .collect::<Vec<Complex<f64>>>();

        let mut res = vec![];
        middle.iter().for_each(|v| {
            res.push(v.re().round());
        });
        middle.iter().for_each(|v| {
            res.push(v.im().round());
        });
        res
    }

    /// Calculates e^{(2i \pi)/2k}
    fn _root_of_unity(k: usize) -> Complex<f64> {
        Complex::from_polar(1.0, (std::f64::consts::PI / (k as f64)))
    }

    /// Calculates powers (\omega_{4(k/2)})^i for i \in [0,k/2)
    fn _powers_of_root_of_unity(k: usize) -> Vec<Complex<f64>> {
        let root = _root_of_unity(k);
        (0..(k / 2))
            .into_iter()
            .map(|exp| root.pow(&(exp as f64)))
            .collect::<Vec<Complex<f64>>>()
    }

    let mut mapped_p1 = _map_from_big_ring_to_small_ring(&p1);
    let mut mapped_p2 = _map_from_big_ring_to_small_ring(&p2);

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(mapped_p1.len());
    fft.process(&mut mapped_p1);
    fft.process(&mut mapped_p2);

    // Hadamard product
    let mut product = mapped_p1
        .iter()
        .zip(mapped_p2.iter())
        .map(|(v1, v2)| v1 * v2)
        .collect::<Vec<Complex<f64>>>();
    let fft = planner.plan_fft_inverse(mapped_p1.len());
    fft.process(&mut product);
    normalize_fft_output(&mut product);

    _map_from_small_ring_to_big_ring(&product)
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use itertools::{izip, Itertools};
    use rand::{thread_rng, Rng};
    use rustfft::num_traits::ToPrimitive;

    use super::*;

    pub(crate) struct Stats<T> {
        pub(crate) samples: Vec<T>,
    }

    impl<T> Default for Stats<T> {
        fn default() -> Self {
            Stats { samples: vec![] }
        }
    }

    impl<T: Copy + ToPrimitive + Debug> Stats<T>
    where
        // T: for<'a> Sum<&'a T>,
        T: for<'a> std::iter::Sum<&'a T> + std::iter::Sum<T>,
    {
        pub(crate) fn new() -> Self {
            Self { samples: vec![] }
        }

        pub(crate) fn mean(&self) -> f64 {
            self.samples.iter().sum::<T>().to_f64().unwrap() / (self.samples.len() as f64)
        }

        pub(crate) fn variance(&self) -> f64 {
            let mean = self.mean();

            // diff
            let diff_sq = self
                .samples
                .iter()
                .map(|v| {
                    let t = v.to_f64().unwrap() - mean;
                    t * t
                })
                .into_iter()
                .sum::<f64>();

            diff_sq / (self.samples.len() as f64 - 1.0)
        }

        pub(crate) fn std_dev(&self) -> f64 {
            self.variance().sqrt()
        }

        pub(crate) fn add_many_samples(&mut self, values: &[T]) {
            self.samples.extend(values.iter());
        }

        pub(crate) fn add_sample(&mut self, value: T) {
            self.samples.push(value)
        }

        pub(crate) fn merge_in(&mut self, other: &Self) {
            self.samples.extend(other.samples.iter());
        }
    }

    fn map_f64_to_z_q(values: &[f64], q: u64) -> Vec<u64> {
        values
            .iter()
            .map(|v| {
                // v.to_bits() % q
                let vint = v.rem_euclid(q as f64).to_i64().unwrap() % (q as i64);
                if vint.is_negative() {
                    q - (vint.abs().to_u64().unwrap())
                } else {
                    vint.abs().to_u64().unwrap()
                }
            })
            .collect()
    }

    fn map_to_i64_from_z_q(values: &[u64], q: u64) -> Vec<i64> {
        values
            .iter()
            .map(|v| {
                if *v > (q >> 1) {
                    (q - v) as i64
                } else {
                    *v as i64
                }
            })
            .collect()
    }

    #[test]
    fn test_nega_cylic_and_cylic_mul() {
        let mut rng = thread_rng();
        let logn = 11;
        let n = 1 << logn;
        let loga = 50;
        let logb = 17;
        let q = 1u64 << loga;

        let mut p1 = vec![];
        let mut p2 = vec![];

        let mut stats = Stats::new();

        const K: usize = 1;

        for _ in 0..K {
            for _ in 0..n {
                p1.push(rng.gen_range(0u64..1u64 << loga));
                p2.push(rng.gen_range(0u64..1u64 << logb));
            }

            let p1_f64 = p1.iter().map(|v| v.to_f64().unwrap()).collect::<Vec<f64>>();
            let p2_f64 = p2.iter().map(|v| v.to_f64().unwrap()).collect::<Vec<f64>>();

            let result = map_f64_to_z_q(&negacylic_mul_with_twist(&p1_f64, &p2_f64), q);
            let expected_result = school_book_negacylic_multiplication_mod_z_q(&p1, &p2, q);
            let diff = izip!(result.iter(), expected_result.iter())
                .map(|(a, b)| (a + q - b) % q)
                .collect_vec();

            stats.add_many_samples(&map_to_i64_from_z_q(&diff, q));
        }

        println!("Error Mean: {}", stats.mean().abs().log2());
        println!("Error log 2 std_dev: {}", stats.std_dev().abs().log2());
        // Eq 35 2021/480
        println!(
            "Expected error log 2 std_dev <= {}",
            (4 * logn + 2 * loga + 2 * logb - 2 * 53 - 3) / 2
        );
    }
}

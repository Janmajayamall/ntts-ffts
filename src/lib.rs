use rand::{thread_rng, Rng};
use rustfft::{
    num_complex::{Complex, ComplexFloat},
    num_traits::Pow,
    FftPlanner,
};

fn school_book_cylic_multiplication(p1: &[f32], p2: &[f32]) -> Vec<f32> {
    let n = p1.len();
    let mut res = vec![0f32; n];
    for i in 0..n {
        for j in 0..i + 1 {
            res[i] += p1[j] * p2[i - j];
        }
        for j in (i + 1)..n {
            res[i] += p1[j] * p2[n - (j - i)];
        }
    }
    res
}

fn school_book_negacylic_multiplication(p1: &[f32], p2: &[f32]) -> Vec<f32> {
    let n = p1.len();
    let mut res = vec![0f32; n];
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

fn normalize_fft_output(p: &mut [Complex<f32>]) {
    let n = p.len() as f32;
    p.iter_mut().for_each(|v| *v = *v / n);
}
/// Multiplies (p1 x p2) mod X^N+1
fn negacylic_mul_with_twist(p1: &[f32], p2: &[f32]) -> Vec<f32> {
    /// Maps $p \in Z[X]/X^N+1$ to element $\in Z[X]/X^N-1$. One cannot map $p$ directly
    /// from big ring to small ring. Instead observe that X^N+1 factors into X^{N/2}-i and X^{N/2}+i.
    /// Ring $Z[X]/X^{N/2}-i$ maps to Z[X]/X^{N/2}-1 by replacing X with $\omega_{4(N/2)}X$
    /// (refer to accompaniying notes to understand why).
    ///
    /// Thus to map $p$ to element in smaller ring, first map it to X^N+1 by consdering upper half coefficients as imaginery counterparts of lower half coefficients. Then element-wise multiply resultling coefficient vector with vector of powers of $\omega_{4(N/2)}$.
    fn _map_from_big_ring_to_small_ring(p: &[f32]) -> Vec<Complex<f32>> {
        let n = p.len();
        let middle = p
            .iter()
            .take(p.len() / 2)
            .zip(p.iter().skip(p.len() / 2))
            .map(|(r, i)| Complex { re: *r, im: *i })
            .collect::<Vec<Complex<f32>>>();

        let powers_of_unity = _powers_of_root_of_unity(n);

        middle
            .iter()
            .zip(powers_of_unity.iter())
            .map(|(v, omega_power)| v * omega_power)
            .collect::<Vec<Complex<f32>>>()
    }

    /// Maps $p \in Z[X]/X^N-1$ to element $\in Z[X]/X^N+1$.
    /// This calculates the inverse of `_map_from_big_ring_to_small_ring`. It first multiplies $p$ with inverse of powers of $\omega_{4(N/2)}$. Then to map
    /// resulting element in $Z[X]/X^N-i$ to element in $Z[X]/X^N+1$, it takes reals in coefficient vector and places them in lower half of new coefficient vector.
    /// Then it considers imaginary counterparts as reals and appends them to the new coefficient vector.
    fn _map_from_small_ring_to_big_ring(p: &[Complex<f32>]) -> Vec<f32> {
        let n = p.len() * 2;

        let power_of_unity_inverse = _powers_of_root_of_unity(n)
            .iter()
            .map(|c| c.inv())
            .collect::<Vec<Complex<f32>>>();

        let middle = p
            .iter()
            .zip(power_of_unity_inverse.iter())
            .map(|(v, omega_power)| v * omega_power)
            .collect::<Vec<Complex<f32>>>();

        let mut res = vec![];
        middle.iter().for_each(|v| {
            res.push(v.re().round());
        });
        middle.iter().for_each(|v| {
            res.push(v.im().round());
        });
        res
    }

    /// Calculates e^{(2i\pi)/4(k/2)}
    fn _root_of_unity(k: usize) -> Complex<f32> {
        Complex::from_polar(1.0, (std::f32::consts::PI / (k as f32)))
    }

    /// Calculates powers (\omega_{4(k/2)})^i for i \in [0,k/2)
    fn _powers_of_root_of_unity(k: usize) -> Vec<Complex<f32>> {
        let root = _root_of_unity(k);
        (0..(k / 2))
            .into_iter()
            .map(|exp| root.pow(&(exp as f32)))
            .collect::<Vec<Complex<f32>>>()
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
        .collect::<Vec<Complex<f32>>>();
    let fft = planner.plan_fft_inverse(mapped_p1.len());
    fft.process(&mut product);
    normalize_fft_output(&mut product);

    _map_from_small_ring_to_big_ring(&product)
}

/// Multiplies (p1 x p2) / (X^N + 1)
fn negacylic_mul_duplicate_mapping(p1: &[f32], p2: &[f32]) -> Vec<f32> {
    /// Maps $p \in Z[X]/X^N+1$ to $p-p(X^N) \in Z[X]/X^{2N}-1$. Since the upper half coefficients
    /// in Big ring are -ve of lower half coefficients, mapping equals negative lower half coefficients
    /// and append at the end.
    fn _preimage_map_big_ring(p: &[f32]) -> Vec<Complex<f32>> {
        let mut buffer1 = p
            .iter()
            .map(|v| Complex { re: *v, im: 0f32 })
            .collect::<Vec<Complex<f32>>>();

        p.iter()
            .for_each(|v| buffer1.push(Complex { re: -*v, im: 0f32 }));

        buffer1
    }

    /// Reduces product of $f(1-f(x^N))g(1-g(x^N)) \in Z[X]/X^{2N}-1$ where $f,g \in Z[X]/X^N+1$ to element $\in Z[X]/X^{N}+1$. Observing
    /// $f(1-f(x^N))g(1-g(x^N)) = 2(fg-fg(X^N))$ where $fg \in Z[X]/X^{N}+1$, fg can be obtained from product by taking lower half of coefficient vector
    /// and dividing the coefficients by 2.
    fn _reduce_big_ring_product_to_small_ring(p: &[Complex<f32>]) -> Vec<f32> {
        let res = p
            .iter()
            .take(p.len() / 2)
            .map(|v| (0.5 * v.re).round())
            .collect::<Vec<f32>>();
        res
    }

    let mut buffer1 = _preimage_map_big_ring(p1);
    let mut buffer2 = _preimage_map_big_ring(p2);

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(buffer1.len());
    fft.process(&mut buffer1);
    fft.process(&mut buffer2);

    // Hadamard product
    let mut product = buffer1
        .iter()
        .zip(buffer2.iter())
        .map(|(v1, v2)| v1 * v2)
        .collect::<Vec<Complex<f32>>>();

    let fft = planner.plan_fft_inverse(buffer1.len());
    fft.process(&mut product);
    normalize_fft_output(&mut product);

    _reduce_big_ring_product_to_small_ring(&product)
}

/// Multiplies (p1 x p2) / (X^N - 1)
fn cyclic_mul(p1: &[f32], p2: &[f32]) -> Vec<f32> {
    let mut planner = FftPlanner::new();
    let mut buffer1 = p1
        .iter()
        .map(|v| Complex { re: *v, im: 0.0 })
        .collect::<Vec<Complex<f32>>>();
    let mut buffer2 = p2
        .iter()
        .map(|v| Complex { re: *v, im: 0.0 })
        .collect::<Vec<Complex<f32>>>();

    let fft = planner.plan_fft_forward(p1.len());
    fft.process(&mut buffer1);
    fft.process(&mut buffer2);

    // Hadamard product
    let mut product = buffer1
        .iter()
        .zip(buffer2.iter())
        .map(|(v1, v2)| v1 * v2)
        .collect::<Vec<Complex<f32>>>();

    let fft = planner.plan_fft_inverse(p1.len());
    fft.process(&mut product);
    normalize_fft_output(&mut product);

    product
        .iter()
        .map(|complex_v| complex_v.re())
        .collect::<Vec<f32>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nega_cylic_and_cylic_mul() {
        let mut rng = thread_rng();
        let n = 8;

        let mut p1 = vec![];
        let mut p2 = vec![];

        for _ in 0..n {
            p1.push(rng.gen_range(0f32..8.0).round());
            p2.push(rng.gen_range(0f32..8.0).round());
        }

        // let p1 = vec![2.0, 5.0, 3.0, 1.0, 1.0, 4.0, 1.0, 4.0];
        // let p2 = vec![5.0, 4.0, 5.0, 1.0, 3.0, 5.0, 1.0, 2.0];

        let result1 = negacylic_mul_duplicate_mapping(&p1, &p2);
        let result2 = negacylic_mul_with_twist(&p1, &p2);
        let expected_result = school_book_negacylic_multiplication(&p1, &p2);
        assert_eq!(result1, expected_result);
        assert_eq!(result2, expected_result);

        let result = cyclic_mul(&p1, &p2);
        let expected_result = school_book_cylic_multiplication(&p1, &p2);
        assert_eq!(result, expected_result);
    }
}

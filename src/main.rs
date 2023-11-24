use rand::{thread_rng, Rng};
// Perform a forward FFT of size 1234
use rustfft::{
    num_complex::{Complex, ComplexFloat},
    FftPlanner,
};

fn school_book_negacylic_multiplication(p1: &[i32], p2: &[i32]) -> Vec<i32> {
    let n = p1.len();
    let mut res = vec![0i32; n];
    for i in 0..n {
        dbg!(i);
        for j in 0..i {
            res[i] += p1[j] * p2[i - j];
        }

        for j in (i + 1)..n {
            res[i] -= p1[j] * p2[n - j];
        }
    }
    res
}

fn preimage_map_big_ring(p: &[i32]) -> Vec<Complex<i32>> {
    let mut buffer1 = p
        .iter()
        .map(|v| Complex { re: *v, im: 0i32 })
        .collect::<Vec<Complex<i32>>>();

    p.iter()
        .for_each(|v| buffer1.push(Complex { re: -*v, im: 0i32 }));

    buffer1
}

fn reduce_big_ring_product_to_small_ring(p: &[Complex<i32>]) -> Vec<i32> {
    let res = p
        .iter()
        .take(p.len() / 2)
        .map(|v| v.re / 2)
        .collect::<Vec<i32>>();
    res
}

fn negacylic_mul_duplicate_mapping(p1: &[i32], p2: &[i32]) -> Vec<i32> {
    let mut planner = FftPlanner::new();
    let mut buffer1 = preimage_map_big_ring(p1);
    let mut buffer2 = preimage_map_big_ring(p2);

    let fft = planner.plan_fft_forward(p1.len());
    fft.process(&mut buffer1);
    fft.process(&mut buffer2);

    // Hadamard product
    let mut product = buffer1
        .iter()
        .zip(buffer2.iter())
        .map(|(v1, v2)| v1 * v2)
        .collect::<Vec<Complex<i32>>>();

    let fft = planner.plan_fft_inverse(p1.len());
    fft.process(&mut product);

    reduce_big_ring_product_to_small_ring(&product)
}

fn main() {
    let mut rng = thread_rng();
    let n = 8;

    let mut p1 = vec![];
    let mut p2 = vec![];

    for _ in 0..n {
        p1.push(rng.gen_range(0i32..8));
        p2.push(rng.gen_range(0i32..8));
    }

    let result = negacylic_mul_duplicate_mapping(&p1, &p2);
    let expected_result = school_book_negacylic_multiplication(&p1, &p2);
    dbg!(&result);
    dbg!(&expected_result);

    assert_eq!(result, expected_result);
}

use fast_modulo::powmod_u64 as modular_exponentiation;
use rayon::prelude::*;
use std::simd::Simd;
use std::simd::num::SimdUint;

pub struct PythagoreanTriples {
    a_values: Vec<u64>,
    b_values: Vec<u64>,
    c_values: Vec<u64>,
}

#[derive(Default)]
struct TemporaryBuffer {
    ascending: Vec<usize>,
    indexes: Vec<usize>,
    values: Vec<u64>,
}

type SimdU64 = Simd::<u64, { crate::SIMD_LANES }>;

impl PythagoreanTriples {
    pub fn new(num_primes: usize) -> Self {
        let mut a_values = Vec::with_capacity(num_primes);
        let mut b_values = Vec::with_capacity(num_primes);
        let mut c_values = Vec::with_capacity(num_primes);

        let mut primes = primal::Primes::all().filter(|p| p % 4 == 1).take(num_primes);
        let mut chunk = Vec::with_capacity(1_000_000);
        let mut tuples = Vec::with_capacity(1_000_000);

        loop {
            chunk.clear();
            chunk.extend(primes.by_ref().map(|p| p as u64).take(1_000_000));
            if chunk.is_empty() { break; }

            tuples.par_extend(chunk.par_iter().map(|&p| (Self::compute(p), p)));

            for ((a, b), c) in tuples.drain(..) {
                a_values.push(a);
                b_values.push(b);
                c_values.push(c);
            }
        }

        Self { a_values, b_values, c_values }
    }

    pub fn len(&self) -> usize {
        self.a_values.len()
    }

    pub fn resize(&mut self, new_len: usize, value: u64) {
        self.a_values.resize(new_len, value);
        self.b_values.resize(new_len, value);
        self.c_values.resize(new_len, value);
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            a_values: Vec::with_capacity(capacity),
            b_values: Vec::with_capacity(capacity),
            c_values: Vec::with_capacity(capacity),
        }
    }

    pub fn push(&mut self, (a, b, c): (u64, u64, u64)) {
        self.a_values.push(a);
        self.b_values.push(b);
        self.c_values.push(c);
    }

    // Apply the Brahmagupta–Fibonacci identity to combine two Pythagorean triples
    // into two new primitive Pythagorean triples for the product of hypotenuses.
    pub fn product(&self, (x, y, z): (u64, u64, u64), output: &mut Self) {
        let num_triples = self.len();
        let existing_len = output.len();
        output.resize(existing_len + num_triples * 2, 0);

        let x_vector = SimdU64::splat(x);
        let y_vector = SimdU64::splat(y);
        let z_vector = SimdU64::splat(z);

        let remainder = num_triples % crate::SIMD_LANES;
        let simd_end = num_triples - remainder;

        for chunk_start in (0..simd_end).step_by(crate::SIMD_LANES) {
            let chunk_end = chunk_start + crate::SIMD_LANES;
            let a_vector = SimdU64::from_slice(&self.a_values[chunk_start..chunk_end]);
            let b_vector = SimdU64::from_slice(&self.b_values[chunk_start..chunk_end]);
            let c_vector = SimdU64::from_slice(&self.c_values[chunk_start..chunk_end]);

            let ax_vector = a_vector * x_vector;
            let ay_vector = a_vector * y_vector;
            let bx_vector = b_vector * x_vector;
            let by_vector = b_vector * y_vector;
            let cz_vector = c_vector * z_vector;

            let first_slot = existing_len + chunk_start * 2;
            let second_slot = first_slot + crate::SIMD_LANES;
            let second_slot_end = second_slot + crate::SIMD_LANES;

            ax_vector.abs_diff(by_vector).copy_to_slice(&mut output.a_values[first_slot..second_slot]);
            (ay_vector + bx_vector).copy_to_slice(&mut output.b_values[first_slot..second_slot]);
            cz_vector.copy_to_slice(&mut output.c_values[first_slot..second_slot]);

            (ax_vector + by_vector).copy_to_slice(&mut output.a_values[second_slot..second_slot_end]);
            ay_vector.abs_diff(bx_vector).copy_to_slice(&mut output.b_values[second_slot..second_slot_end]);
            cz_vector.copy_to_slice(&mut output.c_values[second_slot..second_slot_end]);
        }

        for i in simd_end..num_triples {
            let a = self.a_values[i];
            let b = self.b_values[i];
            let c = self.c_values[i];

            let ax = a * x;
            let ay = a * y;
            let bx = b * x;
            let by = b * y;
            let cz = c * z;

            let first_slot = existing_len + i * 2;
            let second_slot = first_slot + 1;

            output.a_values[first_slot] = ax.abs_diff(by);
            output.b_values[first_slot] = ay + bx;
            output.c_values[first_slot] = cz;

            output.a_values[second_slot] = ax + by;
            output.b_values[second_slot] = ay.abs_diff(bx);
            output.c_values[second_slot] = cz;
        }
    }

    pub fn remove_trivial(&mut self, buffer: &mut TemporaryBuffer) {
        buffer.indexes.clear();

        self.b_values.iter().enumerate().for_each(|(i, &b)| if b != 0 { buffer.indexes.push(i); });
        buffer.values.resize(buffer.indexes.len(), 0);

        for values in [&mut self.a_values, &mut self.b_values, &mut self.c_values] {
            buffer.indexes.iter().enumerate().for_each(|(i, &j)| buffer.values[i] = values[j]);
            values.clear();
            values.extend_from_slice(&buffer.values);
        }
    }

    pub fn sort_and_dedup<F: Fn(&Self, usize) -> O, O: Ord>(&mut self, buffer: &mut TemporaryBuffer, key: F) {
        let num_triples = self.len();

        buffer.ascending.extend(buffer.ascending.len()..num_triples);
        buffer.indexes.resize(num_triples, 0);
        buffer.indexes.copy_from_slice(&buffer.ascending[..num_triples]);

        buffer.indexes.sort_by_key(|&i| key(self, i));
        buffer.indexes.dedup_by_key(|&mut i| key(self, i));
        buffer.values.resize(buffer.indexes.len(), 0);

        for values in [&mut self.a_values, &mut self.b_values, &mut self.c_values] {
            buffer.indexes.iter().enumerate().for_each(|(i, &j)| buffer.values[i] = values[j]);
            values.clear();
            values.extend_from_slice(&buffer.values);
        }
    }

    pub fn sort_and_dedup_by_c_and_a(&mut self, buffer: &mut TemporaryBuffer) {
        self.sort_and_dedup(buffer, |triples, i| (triples.c_values[i], triples.a_values[i]));
    }

    // Use Cornacchia's algorithm to solve a^2 + b^2 = p then apply Euclid's
    // parameterization to find the primitive Pythagorean triple for the prime.
    fn compute(pythagorean_prime: u64) -> (u64, u64) {
        let root = Self::modular_sqrt_of_one_less_than(pythagorean_prime);
        let (m, n) = Self::modified_euclidean_algorithm(pythagorean_prime, root);

        (m * m - n * n, 2 * m * n)
    }

    // Find a quadratic non-residue modulo p. Half of numbers in the field are
    // quadratic non-residues so this is fairly efficient. Use Euler's criterion
    // to check if candidate^2 is congruent to -1 mod p without having to factor.
    fn modular_sqrt_of_one_less_than(pythagorean_prime: u64) -> u64 {
        let k = (pythagorean_prime - 1) / 4;

        for candidate in 2.. {
            let euler_criterion = modular_exponentiation(candidate, 2 * k, pythagorean_prime);

            if euler_criterion == pythagorean_prime - 1 {
                return modular_exponentiation(candidate, k, pythagorean_prime);
            }
        }

        unreachable!()
    }

    // Apply Cornacchia's algorithm which will always terminate for a Pythagorean
    // prime since it can be expressed as the sum of two squares (Fermat's theorem).
    fn modified_euclidean_algorithm(pythagorean_prime: u64, root: u64) -> (u64, u64) {
        let mut a = pythagorean_prime;
        let mut b = root;

        loop {
            let remainder = a % b;

            let sum_of_squares = remainder * remainder + b * b;
            if sum_of_squares == pythagorean_prime { return (b, remainder); }

            a = b;
            b = remainder;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn it_can_compute_the_primitive_pythagorean_triple_for_a_pythagorean_prime() {
        assert_eq!(PythagoreanTriples::compute(5), (3, 4));
        assert_eq!(PythagoreanTriples::compute(13), (5, 12));
        assert_eq!(PythagoreanTriples::compute(17), (15, 8));

        let pythagorean_primes = primal::Primes::all().filter(|p| p % 4 == 1);

        for c in pythagorean_primes.take(100) {
            let (a, b) = PythagoreanTriples::compute(c as u64);
            assert_eq!(a * a + b * b, (c * c) as u64);
        }
    }

    #[test]
    fn it_can_compute_the_first_n_primitive_pythagorean_triples() {
        let triples = PythagoreanTriples::new(100);
        assert_eq!(triples.len(), 100);

        assert_eq!(&triples.a_values[0..5], &[3, 5, 15, 21, 35]);
        assert_eq!(&triples.b_values[0..5], &[4, 12, 8, 20, 12]);
        assert_eq!(&triples.c_values[0..5], &[5, 13, 17, 29, 37]);

        for i in 0..triples.len() {
            let a = triples.a_values[i];
            let b = triples.b_values[i];
            let c = triples.c_values[i];
            assert_eq!(a * a + b * b, c * c);
        }
    }

    #[test]
    fn it_can_calculate_the_product_of_primitive_triples() {
        let mut triples = PythagoreanTriples::new(100);
        let mut output = PythagoreanTriples::with_capacity(203);

        // Existing triples in output should be preserved.
        output.push((3, 4, 5));
        output.push((5, 12, 13));
        output.push((15, 8, 17));

        triples.product((3, 4, 5), &mut output);
        assert_eq!(output.len(), 203);

        assert_eq!(&output.a_values[0..8], &[3, 5, 15, 7, 33, 13, 17, 57]);
        assert_eq!(&output.b_values[0..8], &[4, 12, 8, 24, 56, 84, 144, 176]);
        assert_eq!(&output.c_values[0..8], &[5, 13, 17, 25, 65, 85, 145, 185]);

        for i in 0..output.len() {
            let a = output.a_values[i];
            let b = output.b_values[i];
            let c = output.c_values[i];
            assert_eq!(a * a + b * b, c * c);
        }
    }

    #[test]
    fn it_returns_a_trivial_triple_with_b_set_to_zero_if_the_products_are_the_same() {
        let mut triples = PythagoreanTriples::new(1);
        let mut output = PythagoreanTriples::with_capacity(2);

        triples.product((3, 4, 5), &mut output);
        assert_eq!(output.len(), 2);

        assert_eq!(&output.a_values, &[7, 25]);
        assert_eq!(&output.b_values, &[24, 0]);
    }

    #[test]
    fn it_can_remove_trivial_triples() {
        let mut triples = PythagoreanTriples::with_capacity(5);
        let mut buffer = TemporaryBuffer::default();

        triples.a_values.extend_from_slice(&[3, 5, 5, 13, 0]);
        triples.b_values.extend_from_slice(&[4, 0, 12, 0, 0]);
        triples.c_values.extend_from_slice(&[5, 5, 13, 13, 0]);

        triples.remove_trivial(&mut buffer);
        assert_eq!(&triples.a_values, &[3, 5]);
        assert_eq!(&triples.b_values, &[4, 12]);
        assert_eq!(&triples.c_values, &[5, 13]);
    }

    #[test]
    fn it_can_sort_and_dedup_triples() {
        let mut triples = PythagoreanTriples::with_capacity(5);
        let mut buffer = TemporaryBuffer::default();

        triples.a_values.extend_from_slice(&[3, 5, 3, 5, 3]);
        triples.b_values.extend_from_slice(&[4, 12, 4, 12, 4]);
        triples.c_values.extend_from_slice(&[5, 13, 5, 13, 5]);

        triples.sort_and_dedup_by_c_and_a(&mut buffer);
        assert_eq!(&triples.a_values, &[3, 5]);
        assert_eq!(&triples.b_values, &[4, 12]);
        assert_eq!(&triples.c_values, &[5, 13]);
    }
}

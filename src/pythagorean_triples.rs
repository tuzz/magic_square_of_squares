use fast_modulo::powmod_u64 as modular_exponentiation;
use rayon::prelude::*;

pub struct PythagoreanTriples {
    a_values: Vec<u64>,
    b_values: Vec<u64>,
    c_values: Vec<u64>,
}

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

        for i in 0..100 {
            let a = triples.a_values[i];
            let b = triples.b_values[i];
            let c = triples.c_values[i];
            assert_eq!(a * a + b * b, c * c);
        }
    }
}

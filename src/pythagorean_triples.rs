use fast_modulo::powmod_u64 as modular_exponentiation;

pub struct PythagoreanTriples;

impl PythagoreanTriples {
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
}

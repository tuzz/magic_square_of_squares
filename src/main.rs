#![feature(portable_simd)]

mod composite_number;
mod magic_hourglass;
mod pythagorean_triples;

use composite_number::*;
use magic_hourglass::*;
use pythagorean_triples::*;

const PRINT_FACTORS: bool = false;
const SIMD_LANES: usize = 64;

fn main() {
    let pythagorean_triples = PythagoreanTriples::new(1_000_000);
    let max_factors = u64::MAX.ilog(5) as usize;

    CompositeNumber::new(2..=max_factors, 0..100_000_000_000, pythagorean_triples).for_each(|primitive_start, a_values, b_values, c| {
        detect_magic_hourglass(primitive_start, a_values, b_values, c, |square1, square2, square3, magic_sum| {
             println!("EUREKA! {square1} + {square2} + {square3} = {magic_sum}")
        });
    });
}

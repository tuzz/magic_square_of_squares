#![feature(portable_simd)]

mod composite_number;
mod pythagorean_triples;

use composite_number::*;
use pythagorean_triples::*;

const PRINT_FACTORS: bool = false;
const SIMD_LANES: usize = 64;

fn main() {
    let pythagorean_triples = PythagoreanTriples::new(1_000_000);

    CompositeNumber::new(2..=20, 0..1_000_000_000, pythagorean_triples).for_each(|primitive_start, a_values, b_values, c| {

    });
}

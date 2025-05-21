#![feature(portable_simd)]

mod composite_number;
mod magic_hourglass;
mod magic_patterns;
mod pythagorean_triples;

use composite_number::*;
use magic_hourglass::*;
use magic_patterns::*;
use pythagorean_triples::*;

const NUM_TRIPLES: usize = 500_000_000;
const DETECT_HOURGLASS: bool = true;
const SEARCH_INTERVAL: u64 = 100_000_000_000;
const SIMD_LANES: usize = 64;
const PRINT_FACTORS: bool = false;
const HIDE_KNOWN_SOLUTION: bool = true;

fn main() {
    let pythagorean_triples = PythagoreanTriples::new(NUM_TRIPLES);
    let max_factors = u64::MAX.ilog(5) as usize;

    CompositeNumber::new(2..=max_factors, 0..SEARCH_INTERVAL, pythagorean_triples).for_each(|primitive_start, a_values, b_values, c| {
        if crate::DETECT_HOURGLASS {
            detect_magic_hourglass(primitive_start, a_values, b_values, c);
        } else {
            check_magic_patterns(a_values, b_values, c);
        }
    });
}

#![feature(portable_simd)]

mod composite_number;
mod patterns_16;
mod patterns_234;
mod pythagorean_triples;

use composite_number::*;
use patterns_16::*;
use patterns_234::*;
use pythagorean_triples::*;

const NUM_TRIPLES: usize = 500_000_000;
const SEARCH_MODE: SearchMode = SearchMode::Patterns16;
const SEARCH_INTERVAL: u64 = 100_000_000_000;
const SIMD_LANES: usize = 64;
const PRINT_FACTORS: bool = false;
const HIDE_KNOWN_SOLUTION: bool = true;

#[allow(dead_code)]
enum SearchMode { Patterns16, Patterns234 }

fn main() {
    let pythagorean_triples = PythagoreanTriples::new(NUM_TRIPLES);
    let max_factors = u64::MAX.ilog(5) as usize;

    CompositeNumber::new(2..=max_factors, 0..SEARCH_INTERVAL, pythagorean_triples).for_each(|primitive_start, a_values, b_values, c| {
        match SEARCH_MODE {
            SearchMode::Patterns16 => check_patterns_1_and_6(primitive_start, a_values, b_values, c),
            SearchMode::Patterns234 => check_patterns_2_3_and_4(a_values, b_values, c),
        }
    });
}

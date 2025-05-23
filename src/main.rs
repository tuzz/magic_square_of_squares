#![feature(portable_simd)]

mod composite_number;
mod pattern_1;
mod pattern_234;
mod pattern_6;
mod pythagorean_triples;

use composite_number::*;
use pattern_1::*;
use pattern_234::*;
use pattern_6::*;
use pythagorean_triples::*;

const NUM_TRIPLES: usize = 500_000_000;
const SEARCH_MODE: SearchMode = SearchMode::Pattern1;
const SEARCH_INTERVAL: u64 = 100_000_000_000;
const SIMD_LANES: usize = 64;
const PRINT_FACTORS: bool = false;
const HIDE_KNOWN_SOLUTION: bool = true;

#[allow(dead_code)]
enum SearchMode { Pattern1, Pattern234, Pattern6 }

fn main() {
    let pythagorean_triples = PythagoreanTriples::new(NUM_TRIPLES);
    let max_factors = u64::MAX.ilog(5) as usize;

    CompositeNumber::new(2..=max_factors, 0..SEARCH_INTERVAL, pythagorean_triples).for_each(|primitive_start, a_values, b_values, c| {
        match SEARCH_MODE {
            SearchMode::Pattern1 => check_pattern_1(primitive_start, a_values, b_values, c),
            SearchMode::Pattern234 => check_patterns_234(a_values, b_values, c),
            SearchMode::Pattern6 => check_pattern_6(primitive_start, a_values, b_values, c),
        }
    });
}

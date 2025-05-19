#![feature(portable_simd)]

mod composite_number;
mod pythagorean_triples;

use composite_number::*;
use pythagorean_triples::*;

const PRINT_FACTORS: bool = false;
const SIMD_LANES: usize = 64;

fn main() {
    println!("Hello, world!");
}

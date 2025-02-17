#![feature(portable_simd)]

mod checkpoints;
mod hashing;
mod shared_vec;

use rayon::prelude::*;
use checkpoints::*;
use hashing::*;
use shared_vec::*;
use std::simd::Simd;
use std::sync::{Arc, Mutex};

const SIMD_LANES: usize = 64;

fn main() {
    let (mut squares_by_class, mut sums_by_class, mut centers_to_check, mut next_checkpoint, next_number) = read_checkpoint_or_default();

    for number in next_number.. {
        let square = number as u64 * number as u64;
        if square % 24 != 1 { continue; }

        while let Some(&center) = centers_to_check.front() {
            let center_square = center as u64 * center as u64;
            let center_sum = center_square + center_square;

            if center_sum > square { break; }
            centers_to_check.pop_front();

            let center_class = (center_square % 72 / 24) as usize;
            let complement_class = (6 - center_class) % 3;

            let sums = &mut sums_by_class[complement_class];
            let Some(SharedVec(numbers)) = sums.remove(&hash(center_sum)) else { continue };

            let numbers = Mutex::into_inner(Arc::into_inner(numbers).unwrap()).unwrap();
            // println!("{}, {:?}", center, numbers);
        }

        let center_sum = square + square;
        let center_sum_class = (center_sum % 72 / 24) as usize;

        sums_by_class[center_sum_class].insert(hash(center_sum), SharedVec::default());
        centers_to_check.push_back(number);

        let square_class = (square % 72 / 24) as usize;
        let square_vector = Simd::splat(square);

        sums_by_class.par_iter_mut().enumerate().for_each(|(i, sums)| {
            let residue_class = (6 - square_class + i) % 3;
            let squares = &squares_by_class[residue_class];

            let chunks = squares.par_chunks_exact(SIMD_LANES);
            let remainder = chunks.remainder();

            chunks.for_each(|chunk| {
                let sum_vector = square_vector + Simd::from_slice(chunk);
                for hash in parallel_hash(sum_vector).as_array() {
                    if let Some(vec) = sums.get(hash) {
                        vec.0.lock().unwrap().push(number);
                    }
                }
            });

            for &square2 in remainder {
                let sum = square + square2;
                if let Some(vec) = sums.get_mut(&hash(sum)) {
                    vec.0.lock().unwrap().push(number);
                }
            }
        });

        squares_by_class[square_class].push(square);

        if square >= next_checkpoint {
            write_checkpoint(&squares_by_class, &sums_by_class, &centers_to_check, next_checkpoint, number);
            next_checkpoint += CHECKPOINT_FREQUENCY;
        }
    }
}

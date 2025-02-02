use ahash::AHashMap;

fn main() {
    let target_sum = read_target_sum_from_cli();
    let triples = pythagorean_triples(target_sum);
    let _kernel = magic_square_kernel(&triples);

    for squares in triples.chunks_exact(3) {
        println!("{} = {} + {} + {}", target_sum, squares[0], squares[1], squares[2]);
    }
}

fn read_target_sum_from_cli() -> u64 {
    let target_sum = match std::env::args().nth(1) {
        Some(string) => string.parse::<u64>().unwrap(),
        None => { eprintln!("Usage: ./magic_squares <target_sum>");
            std::process::exit(1);
        }
    };

    if target_sum % 72 != 3 {
        eprintln!("The target sum must be congruent to 3 modulo 72.");
        std::process::exit(1);
    }

    target_sum
}

fn pythagorean_triples(target_sum: u64) -> Vec<u64> {
    let mut pythagorean_triples = vec![];
    let mut squares_by_residue_class = [vec![], vec![], vec![]];

    for number in 1_u64.. {
        // Skip perfect squares that are not congruent to 1 modulo 24.
        let square1 = number * number;
        if square1 % 24 != 1 { continue; }

        // Once the square exceeds the target sum, no more triples exist.
        if square1 >= target_sum { return pythagorean_triples; }

        // Map the residues modulo 72 (1, 25, 49) to indexes (0, 1, 2).
        let residue_class1 = (square1 % 72 / 24) as usize;

        for (residue_class2, squares2) in squares_by_residue_class.iter().enumerate() {
            let residue_class3 = (6 - residue_class1 - residue_class2) % 3;

            // Skip enumerating symmetries in B + C when the arrays are equal.
            if residue_class2 == residue_class3 {
                for (i, square2) in squares2.iter().enumerate() {
                    let partial_sum = square1 + square2;
                    if partial_sum >= target_sum { break; }

                    for square3 in &squares2[0..i] {
                        let magic_sum = partial_sum + square3;
                        if magic_sum > target_sum { break; }
                        if magic_sum == target_sum { pythagorean_triples.extend([square1, *square2, *square3]); }
                    }
                }
            // Otherwise, enumerate all pairwise combinations of the arrays.
            } else {
                let squares3 = &squares_by_residue_class[residue_class3];

                for square2 in squares2 {
                    let partial_sum = square1 + square2;
                    if partial_sum >= target_sum { break; }

                    for square3 in squares3 {
                        // Skip symmetrical triples where B + C are swapped.
                        if square3 > square2 { break; }

                        let magic_sum = partial_sum + square3;
                        if magic_sum > target_sum { break; }
                        if magic_sum == target_sum { pythagorean_triples.extend([square1, *square2, *square3]); }
                    }
                }
            }
        }

        // Skip duplicate perfect squares by pushing after the for loop.
        squares_by_residue_class[residue_class1].push(square1);
    }

    unreachable!();
}

fn magic_square_kernel(triples: &[u64]) -> Option<()> {
    // No magic square exists if fewer than 8 ways to make the magic sum.
    if triples.len() < 8 * 3 { return None; }

    // Count the number of times each square appears in a magic sum.
    let mut occurrences = AHashMap::<u64, u32>::new();
    for &square in triples { *occurrences.entry(square).or_insert(0) += 1; }

    let mut center_candidates = 0;
    let mut corner_candidates = 0;
    let mut edge_candidates = 0;

    for &count in occurrences.values() {
        if count >= 4 { center_candidates += 1; }
        if count >= 3 { corner_candidates += 1; }
        if count >= 2 { edge_candidates += 1; }
    }

    // No magic square exists if there aren't enough of each candidate.
    if center_candidates < 1 { return None; }
    if corner_candidates < 5 { return None; } // Includes the center.
    if edge_candidates < 9 { return None; } // Includes the center and corners.

    None // TODO
}

// These are the residues modulo 72 that are congruent to 1 modulo 24.
const RESIDUES: [u8; 3] = [1, 25, 49];

fn main() {
    let mut squares_by_residue = [const { vec![] }; 72];

    for number in 1_u64.. {
        let square1 = number * number;

        // Skip perfect squares that are not congruent to 1 modulo 24.
        if square1 % 24 != 1 { continue; }
        let residue1 = (square1 % 72) as u8;

        for &residue2 in &RESIDUES {
            // Skip magic sums that are not congruent to 3 modulo 72.
            let squares2 = &squares_by_residue[residue2 as usize];
            let residue3 = (147 - residue1 - residue2) % 72;

            // Skip enumerating symmetries in B + C when the arrays are equal.
            if residue2 == residue3 {
                for (i, square2) in squares2.iter().enumerate() {
                    let partial_sum = square1 + square2;

                    for square3 in &squares2[0..i] {
                        let magic_sum = partial_sum + square3;
                        if magic_sum > 1_000_000_000 { return; }
                    }
                }
            // Otherwise, enumerate all pairwise combinations of the arrays.
            } else {
                let squares3 = &squares_by_residue[residue3 as usize];

                for square2 in squares2 {
                    let partial_sum = square1 + square2;

                    for square3 in squares3 {
                        let magic_sum = partial_sum + square3;
                        if magic_sum > 1_000_000_000 { return; }
                    }
                }
            }
        }

        // Skip duplicate perfect squares by pushing after the for loop.
        squares_by_residue[residue1 as usize].push(square1);
    }
}

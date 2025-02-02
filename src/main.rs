fn main() {
    let mut squares_by_residue_class = [vec![], vec![], vec![]];

    for number in 1_u64.. {
        // Skip perfect squares that are not congruent to 1 modulo 24.
        let square1 = number * number;
        if square1 % 24 != 1 { continue; }

        // Map the residues modulo 72 (1, 25, 49) to indexes (0, 1, 2).
        let residue_class1 = (square1 % 72 / 24) as usize;

        for (residue_class2, squares2) in squares_by_residue_class.iter().enumerate() {
            let residue_class3 = (6 - residue_class1 - residue_class2) % 3;

            // Skip enumerating symmetries in B + C when the arrays are equal.
            if residue_class2 == residue_class3 {
                for (i, square2) in squares2.iter().enumerate() {
                    let partial_sum = square1 + square2;

                    for square3 in &squares2[0..i] {
                        let magic_sum = partial_sum + square3;
                        if magic_sum > 1_000_000_000 { return; }
                    }
                }
            // Otherwise, enumerate all pairwise combinations of the arrays.
            } else {
                let squares3 = &squares_by_residue_class[residue_class3];

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
        squares_by_residue_class[residue_class1].push(square1);
    }
}

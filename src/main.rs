use ahash::AHashMap;

fn main() {
    let target_sum = read_target_sum_from_cli();
    let triples = pythagorean_triples(target_sum);
    let kernel = magic_square_kernel(triples);

    if let Some(kernel) = kernel {
        let filename = write_graph_to_file(target_sum, &kernel);
        println!("{}", filename);
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

fn magic_square_kernel(mut triples: Vec<u64>) -> Option<Vec<u64>> {
    // No magic square exists there are less than 8 ways to make the magic sum.
    if triples.len() < 8 * 3 { return None; }

    // Count the number of times each perfect square appears in a magic sum.
    let mut occurrences = AHashMap::<u64, u32>::new();
    for &square in &triples { *occurrences.entry(square).or_insert(0) += 1; }

    let mut kernel = vec![];

    loop {
        let mut center_candidates = 0;
        let mut corner_candidates = 0;
        let mut edge_candidates = 0;

        for &count in occurrences.values() {
            if count >= 4 { center_candidates += 1; }
            if count >= 3 { corner_candidates += 1; }
            if count >= 2 { edge_candidates += 1; }
        }

        // No magic square exists if there aren't enough of each candidate cell.
        if center_candidates < 1 { return None; }
        if corner_candidates < 5 { return None; } // Includes the center.
        if edge_candidates < 9 { return None; } // Includes the center and corners.

        kernel.clear();

        // Eliminate Pythagorean triples where any of the perfect squares appears
        // less than twice since those triples can't be part of the magic square.
        for squares in triples.chunks_exact(3) {
            let square1 = squares[0];
            let square2 = squares[1];
            let square3 = squares[2];

            if occurrences[&square1] < 2 || occurrences[&square2] < 2 || occurrences[&square3] < 2 {
                *occurrences.get_mut(&square1).unwrap() -= 1;
                *occurrences.get_mut(&square2).unwrap() -= 1;
                *occurrences.get_mut(&square3).unwrap() -= 1;
            } else {
                kernel.extend_from_slice(squares);
            }
        }

        // Re-check the initial condition after eliminating triples.
        let num_remaining = kernel.len();
        if num_remaining < 8 * 3 { return None; }

        // If we didn't manage to eliminate any triples, we've found the kernel.
        if num_remaining == triples.len() {
            return Some(kernel);
        } else {
            // Otherwise replace triples with the subset of triples and iterate again.
            std::mem::swap(&mut triples, &mut kernel);
            occurrences.retain(|_, &mut count| count != 0);
        }
    }
}

fn write_graph_to_file(target_sum: u64, kernel: &[u64]) -> String {
    let mut node_labels = AHashMap::<u64, u16>::new();
    let mut node_targets = AHashMap::<u16, Vec<u16>>::new();

    let mut label: i32 = -1;

    for squares in kernel.chunks_exact(3) {
        let square1 = squares[0];
        let square2 = squares[1];
        let square3 = squares[2];

        // Assign a unique label to each node, starting from 0.
        node_labels.entry(square1).or_insert_with(|| { label += 1; label as u16 });
        node_labels.entry(square2).or_insert_with(|| { label += 1; label as u16 });
        node_labels.entry(square3).or_insert_with(|| { label += 1; label as u16 });

        let label1 = node_labels[&square1];
        let label2 = node_labels[&square2];
        let label3 = node_labels[&square3];

        // Add one undirected edge for each pair of nodes.
        node_targets.entry(label1).or_default().push(label2);
        node_targets.entry(label1).or_default().push(label3);
        node_targets.entry(label2).or_default().push(label3);
    }

    // Write the graph to a file in 'vf' format as specified here:
    // https://github.com/MiviaLab/vf3lib?tab=readme-ov-file#text
    use std::io::Write;
    let mut buffer = vec![];

    // # Number of nodes
    let num_nodes = node_labels.len() as u16;
    writeln!(buffer, "{}\n", num_nodes).unwrap();

    // # Node attributes
    for label in 0..num_nodes {
        writeln!(buffer, "{} 0", label).unwrap();
    }

    // # Edges coming out of node {label}
    for label in 0..node_labels.len() as u16 {
        if let Some(targets) = node_targets.get(&label) {
            let num_targets = targets.len() as u16;
            writeln!(buffer, "\n{}", num_targets).unwrap();

            for target in targets {
                writeln!(buffer, "{} {}", label, target).unwrap();
            }
        } else {
            writeln!(buffer, "\n0").unwrap();
        }
    }

    let tmp = std::env::temp_dir();
    let path = tmp.join(format!("magic_sum_graph_{}.vf", target_sum));

    std::fs::File::create(&path).unwrap().write_all(&buffer).unwrap();
    path.display().to_string()
}

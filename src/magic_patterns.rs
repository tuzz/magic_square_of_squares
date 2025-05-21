use std::cell::RefCell;

// The patterns are from figure 5 of http://www.multimagie.com/Search.pdf#page=2

thread_local! {
    static SQUARES: RefCell<(Vec<u128>, Vec<u128>)> = const { RefCell::new((vec![], vec![])) };
}

pub fn check_magic_patterns(a_values: &[u64], b_values: &[u64], c: u64) {
    if crate::HIDE_KNOWN_SOLUTION && c % 425 == 0 { return; }

    let center = c as u128;
    let center_square = center * center;
    let center_sum = center_square + center_square;
    let magic_sum = center_sum + center_square;

    SQUARES.with_borrow_mut(|(a_squares, b_squares)| {
        a_squares.clear();
        a_squares.extend(a_values.iter().map(|&a| { let a = a as u128; a * a }));

        b_squares.clear();
        b_squares.extend(b_values.iter().map(|&b| { let b = b as u128; b * b }));

        for (i, (&a_square1, &b_square1)) in a_squares.iter().zip(b_squares.iter()).enumerate() {
            let other_a_squares = &a_squares[i + 1..];
            let other_b_squares = &b_squares[i + 1..];

            let a_remainder = magic_sum - a_square1;
            let b_remainder = magic_sum - b_square1;
            let b_minimum = a_square1 - center_square;

            let a_upto = other_a_squares.partition_point(|&s| s < a_remainder);
            let b_upto = other_b_squares.partition_point(|&s| s >= b_minimum);

            for (&a_square2, &b_square2) in other_a_squares[..a_upto].iter().zip(&other_b_squares[..a_upto]) {
                let aa_candidate = a_remainder - a_square2;
                if is_square(aa_candidate) {
                    check_pattern_2(aa_candidate, a_square1, b_square1, a_square2, b_square2, center_square);
                    check_pattern_3_and_4(aa_candidate, a_square1, a_square2, b_square1, b_square2, center_square, magic_sum, center_sum);
                }
            }

            for (&a_square2, &b_square2) in other_a_squares.iter().zip(other_b_squares) {
                let ab_candidate = a_remainder - b_square2;
                if is_square(ab_candidate) {
                    check_pattern_3_and_4(ab_candidate, a_square1, b_square2, b_square1, a_square2, center_square, magic_sum, center_sum);
                }

                let ba_candidate = b_remainder - a_square2;
                if is_square(ba_candidate) {
                    check_pattern_3_and_4(ba_candidate, b_square1, a_square2, a_square1, b_square2, center_square, magic_sum, center_sum);
                }
            }

            for (&a_square2, &b_square2) in other_a_squares[..b_upto].iter().zip(&other_b_squares[..b_upto]) {
                let bb_candidate = b_remainder - b_square2;
                if is_square(bb_candidate) {
                    check_pattern_2(bb_candidate, a_square1, b_square1, a_square2, b_square2, center_square);
                    check_pattern_3_and_4(bb_candidate, b_square1, b_square2, a_square1, a_square2, center_square, magic_sum, center_sum);
                }
            }

            // TODO: pattern 6
        }
    });
}

fn check_pattern_2(top_middle: u128, a_square1: u128, b_square1: u128, a_square2: u128, b_square2: u128, center_square: u128) {
    let middle_left = a_square1 - b_square1;
    if is_square(middle_left) {
        print(a_square1, top_middle, a_square2, middle_left, center_square, 0, b_square1, 0, b_square2);
    }

    let middle_left = a_square1 - b_square2;
    if is_square(middle_left) {
        print(a_square1, top_middle, a_square2, middle_left, center_square, 0, b_square1, 0, b_square2);
    }

    let middle_left = a_square2 - b_square1;
    if is_square(middle_left) {
        print(a_square1, top_middle, a_square2, middle_left, center_square, 0, b_square1, 0, b_square2);
    }

    let middle_left = a_square2 - b_square2;
    if is_square(middle_left) {
        print(a_square1, top_middle, a_square2, middle_left, center_square, 0, b_square1, 0, b_square2);
    }
}

fn check_pattern_3_and_4(top_left: u128, left_square1: u128, left_square2: u128, right_square1: u128, right_square2: u128, center_square: u128, magic_sum: u128, center_sum: u128) {
    let top_middle = magic_sum - top_left - right_square1;
    if is_square(top_middle) {
        print(top_left, top_middle, right_square1, left_square2, center_square, right_square2, left_square1, 0, 0);
    } else {
        let bottom_middle = center_sum - top_middle;
        if is_square(bottom_middle) {
            print(top_left, 0, right_square1, left_square2, center_square, right_square2, left_square1, bottom_middle, 0);
        }
    }

    let top_middle = magic_sum - top_left - right_square2;
    if is_square(top_middle) {
        print(top_left, top_middle, right_square2, left_square1, center_square, right_square1, left_square2, 0, 0);
    } else {
        let bottom_middle = center_sum - top_middle;
        if is_square(bottom_middle) {
            print(top_left, 0, right_square1, left_square2, center_square, right_square2, left_square1, bottom_middle, 0);
        }
    }
}

fn is_square(n: u128) -> bool {
    let root = n.isqrt();
    root * root == n
}

fn print(top_left: u128, top_middle: u128, top_right: u128, middle_left: u128, middle_middle: u128, middle_right: u128, bottom_left: u128, bottom_middle: u128, bottom_right: u128) {
    println!("----------------------------------------------------------------------------------------------------");
    println!("| {top_left:^30} | {top_middle:^30} | {top_right:^30} |");
    println!("|--------------------------------------------------------------------------------------------------|");
    println!("| {middle_left:^30} | {middle_middle:^30} | {middle_right:^30} |");
    println!("|--------------------------------------------------------------------------------------------------|");
    println!("| {bottom_left:^30} | {bottom_middle:^30} | {bottom_right:^30} |");
    println!("----------------------------------------------------------------------------------------------------\n\n");
}

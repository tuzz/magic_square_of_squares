use std::cell::RefCell;

// A "magic hourglass" is defined in this paper: http://www.multimagie.com/Buell.pdf

pub fn detect_magic_hourglass<F: Fn(u128, u128, u128, u128)>(primitive_start: usize, a_values: &[u64], b_values: &[u64], c: u64, callback: F) {
    let center = c as u128;
    let squared_center = center * center;
    let magic_sum = squared_center * 3;

    thread_local!(static SQUARES: RefCell<Vec<u128>> = const { RefCell::new(vec![]) });

    SQUARES.with_borrow_mut(|squares| {
        squares.clear();
        squares.extend(b_values[..primitive_start].iter().rev().map(|&b| { let b = b as u128; b * b }));
        squares.extend(a_values[..primitive_start].iter().map(|&a| { let a = a as u128; a * a }));

        let squares_primitive_start = squares.len();

        squares.extend(b_values[primitive_start..].iter().rev().map(|&b| { let b = b as u128; b * b }));
        squares.extend(a_values[primitive_start..].iter().map(|&a| { let a = a as u128; a * a }));

        for (i, &square1) in squares[squares_primitive_start..].iter().enumerate() {
            let remainder = magic_sum - square1;
            let upto_index1 = squares[..squares_primitive_start].partition_point(|&square| square < remainder);

            for (j, &square2) in squares[..upto_index1].iter().enumerate() {
                let target = remainder - square2;

                if squares[j + 1..squares_primitive_start].binary_search(&target).is_ok() {
                    callback(square1, square2, target, magic_sum);
                }
            }

            let upto_index2 = squares[squares_primitive_start..squares_primitive_start + i].partition_point(|&square| square < remainder);

            for (j, &square2) in squares[squares_primitive_start..squares_primitive_start + upto_index2].iter().enumerate() {
                let target = remainder - square2;

                if squares[j + 1..i].binary_search(&target).is_ok() {
                    callback(square1, square2, target, magic_sum);
                }
            }
        }
    });
}

use std::cell::RefCell;
use crate::patterns_234::print;

// The patterns are from figure 5 of http://www.multimagie.com/Search.pdf#page=2

pub fn check_patterns_1_and_6(primitive_start: usize, a_values: &[u64], b_values: &[u64], c: u64) {
    let center = c as u128;
    let squared_center = center * center;
    let magic_sum = squared_center * 3;

    thread_local! {
        static STATE: RefCell<(Vec<(u128, u128)>, Vec<(u128, u128)>, Vec<u128>, Vec<u128>)> = const { RefCell::new((vec![], vec![], vec![], vec![])) };
    }

    STATE.with_borrow_mut(|(non_primitive, primitive, non_primitive_flat, primitive_flat)| {
        non_primitive.clear();
        non_primitive.extend(a_values[..primitive_start].iter().zip(b_values[..primitive_start].iter()).map(|(&a, &b)| { let a = a as u128; let b = b as u128; (a * a, b * b) }));

        primitive.clear();
        primitive.extend(a_values[primitive_start..].iter().zip(b_values[primitive_start..].iter()).map(|(&a, &b)| { let a = a as u128; let b = b as u128; (a * a, b * b) }));

        non_primitive_flat.clear();
        non_primitive_flat.extend(non_primitive.iter().rev().map(|&(_, b)| b));
        non_primitive_flat.extend(non_primitive.iter().map(|&(a, _)| a));

        primitive_flat.clear();
        primitive_flat.extend(primitive.iter().rev().map(|&(_, b)| b));
        primitive_flat.extend(primitive.iter().map(|&(a, _)| a));

        check_pattern_1(primitive_flat, non_primitive_flat, magic_sum);
        check_pattern_6(primitive, non_primitive, squared_center, magic_sum);
    });
}

fn check_pattern_1(primitive_flat: &[u128], non_primitive_flat: &[u128], magic_sum: u128) {
    for (i, &square1) in primitive_flat.iter().enumerate() {
        let remainder = magic_sum - square1;

        let upto_index1 = primitive_flat[..i].partition_point(|&square| square < remainder);
        let upto_index2 = non_primitive_flat.partition_point(|&square| square < remainder);

        for &square2 in &primitive_flat[..upto_index1] {
            let target = remainder - square2;

            if primitive_flat[..upto_index1].binary_search(&target).is_ok() {
                println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
            };

            if non_primitive_flat[..upto_index2].binary_search(&target).is_ok() {
                println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
            }
        }

        for &square2 in &non_primitive_flat[..upto_index2] {
            let target = remainder - square2;

            if non_primitive_flat[..upto_index2].binary_search(&target).is_ok() {
                println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
            }
        }
    }
}

fn check_pattern_6(primitive: &[(u128, u128)], non_primitive: &[(u128, u128)], squared_center: u128, magic_sum: u128) {
    for (i, &(top_left, bottom_right)) in primitive.iter().enumerate() {
        let remainder1 = magic_sum - top_left;
        let remainder2 = magic_sum - bottom_right;

        let upto_index1 = primitive[..i].partition_point(|&(square, _)| square < remainder1);
        let upto_index2 = non_primitive.partition_point(|&(square, _)| square < remainder1);

        for &(middle_left, middle_right) in &primitive[..upto_index1] {
            let bottom_left = remainder1 - middle_left; // smaller
            let top_right = remainder2 - middle_right; // bigger

            let bottom_middle = remainder2 - bottom_left; // bigger
            let top_middle = remainder1 - top_right; // smaller
            let target = (bottom_middle, top_middle);

            if primitive[..upto_index1].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
            };

            if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
            }

            // Check the symmetrical case where the middle row is swapped.
            let bottom_left = remainder1 - middle_right; // bigger or smaller
            let top_right = remainder2 - middle_left; // bigger or smaller

            let bottom_middle = remainder2 - bottom_left; // bigger or smaller
            let top_middle = remainder1 - top_right; // bigger or smaller
            let target = if bottom_middle > top_middle { (bottom_middle, top_middle) } else { (top_middle, bottom_middle) };

            if primitive[..upto_index1].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
            };

            if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
            }
        }

        for &(middle_left, middle_right) in &non_primitive[..upto_index2] {
            let bottom_left = remainder1 - middle_left; // smaller
            let top_right = remainder2 - middle_right; // bigger

            let bottom_middle = remainder2 - bottom_left; // bigger
            let top_middle = remainder1 - top_right; // smaller
            let target = (bottom_middle, top_middle);

            if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
            }

            // Check the symmetrical case where the middle row is swapped.
            let bottom_left = remainder1 - middle_right; // bigger or smaller
            let top_right = remainder2 - middle_left; // bigger or smaller

            let bottom_middle = remainder2 - bottom_left; // bigger or smaller
            let top_middle = remainder1 - top_right; // bigger or smaller
            let target = if bottom_middle > top_middle { (bottom_middle, top_middle) } else { (top_middle, bottom_middle) };

            if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
            }
        }
    }
}

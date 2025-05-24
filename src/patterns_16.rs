use std::cell::RefCell;
use crate::patterns_234::print;

// The patterns are from figure 5 of http://www.multimagie.com/Search.pdf#page=2

pub fn check_patterns_1_and_6(primitive_start: usize, a_values: &[u64], b_values: &[u64], c: u64) {
    let center = c as u128;
    let squared_center = center * center;
    let magic_sum = squared_center * 3;

    thread_local! {
        static TRIPLES: RefCell<(Vec<(u128, u128)>, Vec<(u128, u128)>)> = const { RefCell::new((vec![], vec![])) };
    }

    TRIPLES.with_borrow_mut(|(non_primitive, primitive)| {
        non_primitive.clear();
        non_primitive.extend(a_values[..primitive_start].iter().zip(b_values[..primitive_start].iter()).map(|(&a, &b)| { let a = a as u128; let b = b as u128; (a * a, b * b) }));

        primitive.clear();
        primitive.extend(a_values[primitive_start..].iter().zip(b_values[primitive_start..].iter()).map(|(&a, &b)| { let a = a as u128; let b = b as u128; (a * a, b * b) }));

        for (i, &(top_left, bottom_right)) in primitive.iter().enumerate() {
            let remainder1 = magic_sum - top_left;
            let remainder2 = magic_sum - bottom_right;

            let upto_index1 = primitive[..i].partition_point(|&(square, _)| square < remainder1);
            let upto_index2 = non_primitive.partition_point(|&(square, _)| square < remainder1);

            for &(middle_left, middle_right) in &primitive[..upto_index1] {
                let bottom_left = remainder1 - middle_left; // smaller
                let top_right = remainder2 - middle_right; // bigger (increasing)
                let pattern_1_target = (top_right, bottom_left);

                let bottom_middle = remainder2 - bottom_left; // bigger
                let Some(top_middle) = remainder1.checked_sub(top_right) else { break }; // smaller
                let pattern_6_target = (bottom_middle, top_middle);

                if primitive[..i].binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                };

                if non_primitive.binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                }

                if primitive[..i].binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                };

                if non_primitive.binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                }
            }

            for &(middle_left, middle_right) in &non_primitive[..upto_index2] {
                let bottom_left = remainder1 - middle_left; // smaller
                let top_right = remainder2 - middle_right; // bigger (increasing)
                let pattern_1_target = (top_right, bottom_left);

                let bottom_middle = remainder2 - bottom_left; // bigger
                let Some(top_middle) = remainder1.checked_sub(top_right) else { break }; // smaller
                let pattern_6_target = (bottom_middle, top_middle);

                if non_primitive.binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                }

                if non_primitive.binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_left, squared_center, middle_right, bottom_left, bottom_middle, bottom_right);
                }
            }

            // Check the symmetrical case where (middle_left, middle_right) are swapped.
            for &(middle_right, middle_left) in &primitive[..i] {
                let bottom_left = remainder1 - middle_left; // bigger or smaller
                let top_right = remainder2 - middle_right; // bigger or smaller (decreasing)
                let pattern_1_target = if top_right > bottom_left { (top_right, bottom_left) } else { (bottom_left, top_right) };

                let bottom_middle = remainder2 - bottom_left; // bigger or smaller
                let Some(top_middle) = remainder1.checked_sub(top_right) else { continue }; // bigger or smaller
                let pattern_6_target = if bottom_middle > top_middle { (bottom_middle, top_middle) } else { (top_middle, bottom_middle) };

                if primitive[..i].binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                };

                if non_primitive.binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                }

                if primitive[..i].binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                };

                if non_primitive.binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                }
            }

            for &(middle_right, middle_left) in non_primitive.iter() {
                let bottom_left = remainder1 - middle_left; // bigger or smaller
                let top_right = remainder2 - middle_right; // bigger or smaller (decreasing)
                let pattern_1_target = if top_right > bottom_left { (top_right, bottom_left) } else { (bottom_left, top_right) };

                let bottom_middle = remainder2 - bottom_left; // bigger or smaller
                let Some(top_middle) = remainder1.checked_sub(top_right) else { continue }; // bigger or smaller
                let pattern_6_target = if bottom_middle > top_middle { (bottom_middle, top_middle) } else { (top_middle, bottom_middle) };

                if non_primitive.binary_search(&pattern_1_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                }

                if non_primitive.binary_search(&pattern_6_target).is_ok() {
                    print(top_left, top_middle, top_right, middle_right, squared_center, middle_left, bottom_left, bottom_middle, bottom_right);
                }
            }
        }
    });
}

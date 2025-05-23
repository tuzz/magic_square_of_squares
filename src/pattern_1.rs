use std::cell::RefCell;

// The patterns are from figure 5 of http://www.multimagie.com/Search.pdf#page=2
// A "magic hourglass" is defined in this paper: http://www.multimagie.com/Buell.pdf

pub fn check_pattern_1(primitive_start: usize, a_values: &[u64], b_values: &[u64], c: u64) {
    let center = c as u128;
    let squared_center = center * center;
    let magic_sum = squared_center * 3;

    thread_local! {
        static SQUARES: RefCell<(Vec<u128>, Vec<u128>)> = const { RefCell::new((vec![], vec![])) };
    }

    SQUARES.with_borrow_mut(|(non_primitive, primitive)| {
        non_primitive.clear();
        non_primitive.extend(b_values[..primitive_start].iter().rev().map(|&b| { let b = b as u128; b * b }));
        non_primitive.extend(a_values[..primitive_start].iter().map(|&a| { let a = a as u128; a * a }));

        primitive.clear();
        primitive.extend(b_values[primitive_start..].iter().rev().map(|&b| { let b = b as u128; b * b }));
        primitive.extend(a_values[primitive_start..].iter().map(|&a| { let a = a as u128; a * a }));

        for (i, &square1) in primitive.iter().enumerate() {
            let remainder = magic_sum - square1;

            let upto_index1 = primitive[..i].partition_point(|&square| square < remainder);
            let upto_index2 = non_primitive.partition_point(|&square| square < remainder);

            for &square2 in &primitive[..upto_index1] {
                let target = remainder - square2;

                if primitive[..upto_index1].binary_search(&target).is_ok() {
                    println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
                };

                if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                    println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
                }
            }

            for &square2 in &non_primitive[..upto_index2] {
                let target = remainder - square2;

                if non_primitive[..upto_index2].binary_search(&target).is_ok() {
                    println!("EUREKA! {square1} + {square2} + {target} = {magic_sum}")
                }
            }
        }
    });
}

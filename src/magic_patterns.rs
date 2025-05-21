use std::cell::RefCell;

// The patterns are from figure 5 of http://www.multimagie.com/Search.pdf#page=2

pub struct Cell {
    pub value: u128,
    pub magic_sum_minus_value: u128,
}

impl Cell {
    pub fn new(value: u128, magic_sum: u128) -> Self {
        Self { value, magic_sum_minus_value: magic_sum - value }
    }

    pub fn squared(n: u64, magic_sum: u128) -> Self {
        let n = n as u128;
        Self::new(n * n, magic_sum)
    }
}

thread_local! {
    static CELLS: RefCell<(Vec<Cell>, Vec<Cell>)> = const { RefCell::new((vec![], vec![])) };
}

pub fn check_magic_patterns(a_values: &[u64], b_values: &[u64], c: u64) {
    if crate::HIDE_KNOWN_SOLUTION && c % 425 == 0 { return; }

    let center = c as u128;
    let center_square = center * center;
    let center_sum = center_square + center_square;
    let magic_sum = center_sum + center_square;
    let center_cell = Cell::new(center_square, magic_sum);

    CELLS.with_borrow_mut(|(a_cells, b_cells)| {
        a_cells.clear();
        a_cells.extend(a_values.iter().map(|&a| Cell::squared(a, magic_sum)));

        b_cells.clear();
        b_cells.extend(b_values.iter().map(|&b| Cell::squared(b, magic_sum)));

        for (i, (a_cell1, b_cell1)) in a_cells.iter().zip(b_cells.iter()).enumerate() {
            let other_a_cells = &a_cells[i + 1..];
            let other_b_cells = &b_cells[i + 1..];

            let upto_index = other_a_cells.partition_point(|a_cell| a_cell.value < a_cell1.magic_sum_minus_value);

            // a_cell1 and a_cell2 may be on the same line:
            for (a_cell2, b_cell2) in other_a_cells[..upto_index].iter().zip(&other_b_cells[..upto_index]) {
                check_patterns_1_and_2(&center_cell, a_cell1, b_cell1, a_cell2, b_cell2);
                check_patterns_3_4_and_6(&center_cell, a_cell1, b_cell1, b_cell2, a_cell2);
                check_patterns_3_4_and_6(&center_cell, a_cell2, b_cell2, b_cell1, a_cell1);
            }

            // a_cell1 and b_cell2 may be on the same line:
            for (a_cell2, b_cell2) in other_a_cells.iter().zip(other_b_cells) {
                check_patterns_3_4_and_6(&center_cell, a_cell1, b_cell1, a_cell2, b_cell2);
                check_patterns_3_4_and_6(&center_cell, b_cell2, a_cell2, b_cell1, a_cell1);
            }
        }
    });
}

fn check_patterns_1_and_2(center: &Cell, top_left: &Cell, bottom_right: &Cell, top_right: &Cell, bottom_left: &Cell) {
    let top_middle = top_left.magic_sum_minus_value - top_right.value;
    let middle_left = top_left.magic_sum_minus_value - bottom_left.value;
    let Some(bottom_middle) = bottom_left.magic_sum_minus_value.checked_sub(bottom_right.value) else { return };

    let mut num_squares = 5;
    if is_square(top_middle) { num_squares += 1; }
    if is_square(middle_left) { num_squares += 1; }
    if is_square(bottom_middle) { num_squares += 1; }
    if num_squares < 6 { return; }

    let middle_right = center.magic_sum_minus_value - middle_left;
    if is_square(middle_right) { num_squares += 1; }
    if num_squares < 7 { return; }

    println!("{} | {} | {}", top_left.value, top_middle, top_right.value);
    println!("{} | {} | {}", middle_left, center.value, middle_right);
    println!("{} | {} | {}\n", bottom_left.value, bottom_middle, bottom_right.value);
}

fn check_patterns_3_4_and_6(center: &Cell, top_right: &Cell, bottom_left: &Cell, middle_left: &Cell, middle_right: &Cell) {
    let top_left = bottom_left.magic_sum_minus_value - middle_left.value;
    let bottom_right = top_right.magic_sum_minus_value - middle_right.value;
    let Some(top_middle) = top_right.magic_sum_minus_value.checked_sub(top_left) else { return };

    let mut num_squares = 5;
    if is_square(top_left) { num_squares += 1; }
    if is_square(bottom_right) { num_squares += 1; }
    if is_square(top_middle) { num_squares += 1; }
    if num_squares < 6 { return; }

    let Some(bottom_middle) = bottom_left.magic_sum_minus_value.checked_sub(bottom_right) else { return };
    if is_square(bottom_middle) { num_squares += 1; }
    if num_squares < 7 { return; }

    println!("{} | {} | {}", top_left, top_middle, top_right.value);
    println!("{} | {} | {}", middle_left.value, center.value, middle_right.value);
    println!("{} | {} | {}\n", bottom_left.value, bottom_middle, bottom_right);
}

fn is_square(n: u128) -> bool {
    let root = n.isqrt();
    root * root == n
}

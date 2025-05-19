use crate::{PythagoreanTriples, TemporaryBuffer};
use std::ops::{Range, RangeInclusive};

pub struct CompositeNumber {
    non_final_terms: Box<[NonFinalTerm]>,
    final_term_start_index: usize,
    final_term_end_index: usize,
    search_range: Range<u64>,
    pythagorean_triples: PythagoreanTriples,
    temporary_buffer: TemporaryBuffer,
}

struct NonFinalTerm {
    current_triple: (u64, u64, u64, u32),
    cumulative_product: u64,
    triples_powerset: PythagoreanTriples,
    next_index: usize,
    end_index: usize,
}

impl CompositeNumber {
    pub fn new(num_factors: RangeInclusive<usize>, start_range: Range<u64>, pythagorean_triples: PythagoreanTriples) -> Self {
        let min_factors = *num_factors.start();
        let max_factors = *num_factors.end();
        assert!(min_factors >= 2);

        let mut composite_number = Self {
            non_final_terms: (0..max_factors - 1).map(NonFinalTerm::new).collect(),
            final_term_start_index: 0,
            final_term_end_index: 0,
            search_range: start_range,
            pythagorean_triples,
            temporary_buffer: TemporaryBuffer::default(),
        };

        composite_number.next_non_final_term(max_factors - min_factors);
        composite_number
    }

    pub fn next_available_term(&mut self) -> bool {
        (0..self.non_final_terms.len()).rev().any(|i| self.next_non_final_term(i))
    }

    pub fn next_non_final_term(&mut self, term_index: usize) -> bool {
        let num_terms = self.non_final_terms.len() + 1;
        let max_value = self.search_range.end - 1;

        let (previous_terms, next_terms) = self.non_final_terms.split_at_mut(term_index);
        let previous_term = previous_terms.last();
        let (previous_product, previous_c, previous_f) = previous_term.map(|t| (t.cumulative_product, t.current_triple.2, t.current_triple.3)).unwrap_or((1, 1, 0));
        let previous_powerset = previous_term.map(|t| &t.triples_powerset);
        let current_term = next_terms.first_mut().unwrap();

        if current_term.next_index < current_term.end_index {
            let c = self.pythagorean_triples.c_values[current_term.next_index];
            let mut product = previous_product * c;

            let mut next_max = Self::max_value_for_term(term_index + 1, num_terms, product, max_value);
            if next_max < c { return false; }

            let a = self.pythagorean_triples.a_values[current_term.next_index];
            let b = self.pythagorean_triples.b_values[current_term.next_index];
            let f = if c == previous_c { previous_f } else { previous_f + 1 };

            current_term.current_triple = (a, b, c, f);
            current_term.cumulative_product = product;
            current_term.update_triples_powerset(previous_powerset);
            current_term.triples_powerset.sort_and_dedup_by_c_and_a(&mut self.temporary_buffer);
            current_term.next_index += 1;
            let next_index = current_term.next_index;

            for i in term_index + 1..self.non_final_terms.len() {
                let (previous_terms, next_terms) = self.non_final_terms.split_at_mut(i);
                let previous_powerset = previous_terms.last().map(|t| &t.triples_powerset);
                let next_term = next_terms.first_mut().unwrap();

                product *= c;

                next_term.current_triple = (a, b, c, f);
                next_term.cumulative_product = product;
                next_term.update_triples_powerset(previous_powerset);
                next_term.triples_powerset.sort_and_dedup_by_c_and_a(&mut self.temporary_buffer);
                next_term.next_index = next_index;
                next_term.end_index = self.pythagorean_triples.c_values.partition_point(|&c| c <= next_max);

                next_max = Self::max_value_for_term(i + 1, num_terms, product, max_value);
            }

            let next_min = c.max(self.search_range.start.div_ceil(product));

            self.final_term_end_index = self.pythagorean_triples.c_values.partition_point(|&c| c <= next_max);
            self.final_term_start_index = self.pythagorean_triples.c_values[..self.final_term_end_index].partition_point(|&c| c < next_min);
            self.non_final_terms.last_mut().unwrap().triples_powerset.remove_trivial(&mut self.temporary_buffer);

            true
        } else {
            false
        }
    }

    fn max_value_for_term(term_index: usize, num_terms: usize, previous_product: u64, max_value: u64) -> u64 {
        let remaining_multiple = max_value / previous_product;
        let remaining_terms = num_terms - term_index;

        match remaining_terms {
            1 => remaining_multiple,
            2 => remaining_multiple.isqrt(),
            _ => (remaining_multiple as f64).powf(1. / remaining_terms as f64).floor() as u64,
        }
    }

    #[cfg(test)]
    fn non_final_factors(&self) -> Vec<u64> {
        self.non_final_terms.iter().map(|t| t.current_triple.2).collect()
    }

    fn final_factors(&self) -> Vec<u64> {
        self.pythagorean_triples.c_values[self.final_term_start_index..self.final_term_end_index].to_vec()
    }
}

impl NonFinalTerm {
    fn new(_term_index: usize) -> Self {
        Self {
            current_triple: (0, 0, 1, 0),
            cumulative_product: 1,
            triples_powerset: PythagoreanTriples::new(0),
            next_index: 0,
            end_index: usize::MAX,
        }
    }

    pub fn update_triples_powerset(&mut self, previous_powerset: Option<&PythagoreanTriples>) {
        self.triples_powerset.clear();
        self.triples_powerset.push(self.current_triple);

        if let Some(previous_powerset) = previous_powerset {
            self.triples_powerset.extend(previous_powerset);
            previous_powerset.product(self.current_triple, &mut self.triples_powerset);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn it_can_advance_through_each_non_final_term_ensuring_lexical_ordering() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=3, 0..1000, pythagorean_triples);
        assert_eq!(composite_number.non_final_factors(), &[1, 5]);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[1, 13]);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[1, 17]);

        composite_number.next_non_final_term(0);
        assert_eq!(composite_number.non_final_factors(), &[5, 5]);
    }

    #[test]
    fn it_returns_false_when_the_search_range_has_been_exhausted() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=3, 0..1000, pythagorean_triples);
        assert_eq!(composite_number.non_final_factors(), &[1, 5]);
        assert!(5 * 5 < 1000);

        assert!(composite_number.next_non_final_term(1));
        assert_eq!(composite_number.non_final_factors(), &[1, 13]);
        assert!(13 * 13 < 1000);

        assert!(composite_number.next_non_final_term(1));
        assert_eq!(composite_number.non_final_factors(), &[1, 17]);
        assert!(17 * 17 < 1000);

        assert!(composite_number.next_non_final_term(1));
        assert_eq!(composite_number.non_final_factors(), &[1, 29]);
        assert!(29 * 29 < 1000);

        assert!(!composite_number.next_non_final_term(1));
        assert_eq!(composite_number.non_final_factors(), &[1, 29]);
        assert!(37 * 37 >= 1000);

        assert!(composite_number.next_non_final_term(0));
        assert_eq!(composite_number.non_final_factors(), &[5, 5]);
        assert!(5 * 5 * 5 < 1000);

        assert!(!composite_number.next_non_final_term(0));
        assert_eq!(composite_number.non_final_factors(), &[5, 5]);
        assert!(13 * 13 * 13 >= 1000);
    }

    #[test]
    fn it_sets_the_current_triple_of_each_non_final_term() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=3, 0..1000, pythagorean_triples);

        assert_eq!(composite_number.non_final_factors(), &[1, 5]);
        assert_eq!(composite_number.non_final_terms[0].current_triple, (0, 0, 1, 0));
        assert_eq!(composite_number.non_final_terms[1].current_triple, (3, 4, 5, 1));

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[1, 13]);
        assert_eq!(composite_number.non_final_terms[0].current_triple, (0, 0, 1, 0));
        assert_eq!(composite_number.non_final_terms[1].current_triple, (5, 12, 13, 1));

        composite_number.next_non_final_term(0);
        assert_eq!(composite_number.non_final_factors(), &[5, 5]);
        assert_eq!(composite_number.non_final_terms[0].current_triple, (3, 4, 5, 1));
        assert_eq!(composite_number.non_final_terms[1].current_triple, (3, 4, 5, 1));

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[5, 13]);
        assert_eq!(composite_number.non_final_terms[0].current_triple, (3, 4, 5, 1));
        assert_eq!(composite_number.non_final_terms[1].current_triple, (5, 12, 13, 2));
    }

    #[test]
    fn it_calculates_the_cumulative_product_of_each_non_final_term() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=3, 0..1000, pythagorean_triples);

        assert_eq!(composite_number.non_final_factors(), &[1, 5]);
        assert_eq!(composite_number.non_final_terms[0].cumulative_product, 1);
        assert_eq!(composite_number.non_final_terms[1].cumulative_product, 5);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[1, 13]);
        assert_eq!(composite_number.non_final_terms[0].cumulative_product, 1);
        assert_eq!(composite_number.non_final_terms[1].cumulative_product, 13);

        composite_number.next_non_final_term(0);
        assert_eq!(composite_number.non_final_factors(), &[5, 5]);
        assert_eq!(composite_number.non_final_terms[0].cumulative_product, 5);
        assert_eq!(composite_number.non_final_terms[1].cumulative_product, 25);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[5, 13]);
        assert_eq!(composite_number.non_final_terms[0].cumulative_product, 5);
        assert_eq!(composite_number.non_final_terms[1].cumulative_product, 65);
    }

    #[test]
    fn it_computes_the_triples_powerset_for_each_non_final_term() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=3, 0..1000, pythagorean_triples);

        assert_eq!(composite_number.non_final_factors(), &[1, 5]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.a_values, &[3]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.b_values, &[4]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.c_values, &[5]);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[1, 13]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.a_values, &[5]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.b_values, &[12]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.c_values, &[13]);

        composite_number.next_non_final_term(0);
        assert_eq!(composite_number.non_final_factors(), &[5, 5]); // Duplicate factors.
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.a_values, &[3]);
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.b_values, &[4]);
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.c_values, &[5]);

        assert_eq!(composite_number.non_final_terms[1].triples_powerset.a_values, &[3, 7]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.b_values, &[4, 24]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.c_values, &[5, 25]);

        composite_number.next_non_final_term(1);
        assert_eq!(composite_number.non_final_factors(), &[5, 13]); // Distinct factors.
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.a_values, &[3]);
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.b_values, &[4]);
        assert_eq!(composite_number.non_final_terms[0].triples_powerset.c_values, &[5]);

        assert_eq!(composite_number.non_final_terms[1].triples_powerset.a_values, &[3, 5, 33, 63]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.b_values, &[4, 12, 56, 16]);
        assert_eq!(composite_number.non_final_terms[1].triples_powerset.c_values, &[5, 13, 65, 65]);
    }

    #[test]
    fn it_can_fully_exhaust_the_search_range() {
        let pythagorean_triples = PythagoreanTriples::new(100);
        let mut composite_number = CompositeNumber::new(2..=4, 485..1000, pythagorean_triples);
        assert_eq!(composite_number.non_final_factors(), &[1, 1, 5]);
        assert_eq!(composite_number.final_factors(), &[97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[1, 1, 13]);
        assert_eq!(composite_number.final_factors(), &[41, 53, 61, 73]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[1, 1, 17]);
        assert_eq!(composite_number.final_factors(), &[29, 37, 41, 53]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[1, 1, 29]);
        assert_eq!(composite_number.final_factors(), &[29]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[1, 5, 5]);
        assert_eq!(composite_number.final_factors(), &[29, 37]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[1, 5, 13]);
        assert_eq!(composite_number.final_factors(), &[13]);

        assert!(composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[5, 5, 5]);
        assert_eq!(composite_number.final_factors(), &[5]);

        assert!(!composite_number.next_available_term());
        assert_eq!(composite_number.non_final_factors(), &[5, 5, 5]);
        assert_eq!(composite_number.final_factors(), &[5]);
    }
}

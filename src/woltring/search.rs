pub(crate) fn find_knot_interval(knots: &Vec<f64>, point: f64, knot_guess: usize) -> usize {
    if point < knots[0] {
        return 0;
    }
    let num_knots = knots.len();
    if point >= knots[num_knots - 1] {
        return num_knots;
    }

    let mut knot_interval = knot_guess.max(1);
    if knot_interval >= num_knots {
        knot_interval = num_knots - 1;
    }

    let (mut lower_index, mut upper_index);
    // Often L will be in an interval adjoining the interval found in a previous call to search
    if point >= knots[knot_interval - 1] {
        if point < knots[knot_interval] {
            return knot_interval;
        } else {
            knot_interval += 1;
            if point < knots[knot_interval] {
                return knot_interval;
            }
            lower_index = knot_interval + 1;
            upper_index = num_knots;
        }
    } else {
        knot_interval -= 1;
        if point >= knots[knot_interval - 1] {
            return knot_interval;
        } else {
            lower_index = 1;
            upper_index = knot_interval;
        }
    }

    // Binary search
    loop {
        knot_interval = (lower_index + upper_index) / 2;
        if upper_index - lower_index <= 1 {
            return knot_interval;
        }
        if point < knots[knot_interval - 1] {
            upper_index = knot_interval;
        } else {
            lower_index = knot_interval;
        }
    }
}
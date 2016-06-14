//! A type for specifying recurrences and iterating over the solutions.
//!
//! Recurrences are only returned by more complex hyperbolic cases of
//! diophantine equations.

use std::fmt;

/// A pair of recurrence equations which can be used to generate successive terms
/// in sequence.
///
/// The results returned a pairs of x, y coordinates which satisfy an equation.
///
/// ## Note
/// This iterator is infinite.
///
/// ## Examples
/// ```text
/// X0 = 1
/// Y0 = 0
/// ```
///
/// and
///
/// ```text
/// Yn+1 = 2 * Xn + 7 * Yn
/// Xn+1 = Xn + 3 * Yn + 1
/// ```
#[derive(Clone, PartialEq)]
pub struct RecurrenceEquation {
    /// Is this the first call
    // We need this flag so we are not pre-emptively generating the next
    // element required, which could cause off-by-one errors on overflow
    // for example.
    first_call: bool,

    /// Current value of x
    x: i64,

    /// Current value of y
    y: i64,

    /// The coefficient on the first recurrence term for X_{n+1}
    x1: i64,

    /// The coefficient on the second recurrence term for X_{n+1}
    x2: i64,

    /// A constant term added to the Xn recurrence term
    x3: i64,

    /// The coefficient on the first recurrence term for Y_{n+1}
    y1: i64,

    /// The coefficient on the second recurrence term for Y_{n+1}
    y2: i64,

    /// A constant term added to the Yn recurrence term
    y3: i64,
}

impl RecurrenceEquation {
    /// Generate a new recurrence equation with the specified coefficients.
    ///
    /// The input represents a pair of recurrences of the following form:
    ///
    /// ```text
    /// X0 = $x
    /// Y0 = $y
    ///
    /// and
    ///
    /// Yn = $x1 * Xn + $x2 * Yn + $x3
    /// Xn = $y1 * Xn + $y2 * Yn + $y3
    /// ```
    // A user does not need to know how to create this, but we still want to
    // document that is an iterator and how it can be used.
    #[doc(hidden)]
    pub fn new(x: i64, y: i64,
               (x1, x2, x3): (i64, i64, i64),
               (y1, y2, y3): (i64, i64, i64)) -> RecurrenceEquation {
        RecurrenceEquation {
            first_call: true, x: x, y: y,
            x1: x1, x2: x2, x3: x3, y1: y1, y2: y2, y3: y3
        }
    }
}

impl Iterator for RecurrenceEquation {
    type Item = (i64, i64);

    fn next(&mut self) -> Option<(i64, i64)> {
        if self.first_call {
            self.first_call = false;
        } else {
            // Generate the next value of x and y
            let xx = self.x;
            let yy = self.y;

            // The recurrences we deal with are all of a particular form, so we
            // do not really require boxing closures and the like for more
            // generic recurrence types.
            self.x = xx * self.x1 + yy * self.x2 + self.x3;
            self.y = xx * self.y1 + yy * self.y2 + self.y3;
        }

        Some((self.x, self.y))
    }
}

impl fmt::Debug for RecurrenceEquation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Recurrence({}, {}, {}; {}, {}, {})",
               self.x1, self.x2, self.x3, self.y1, self.y2, self.y3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recurrence_01() {
        // x0 = 1
        // y0 = 0
        //
        // xn+1 = 5*xn + 2*yn
        // yn+1 = 2*xn + 3*yn - 1
        let mut eq = RecurrenceEquation::new(1, 0, (5, 2, 0), (2, 3, -1));

        assert_eq!((1, 0), eq.next().unwrap());
        assert_eq!((5, 1), eq.next().unwrap());
        assert_eq!((27, 12), eq.next().unwrap());
    }
}

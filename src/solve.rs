//! The `solve` function is the main generic entry-point. This delegates
//! to the specific `solve_X` functions.

#![allow(unused_variables)]

use helper;
use recurrence::RecurrenceEquation;

/// Represents a solution that may be returned when evaluating an equation.
///
/// # Todo
/// Recurrences feel more natural as being a `Form`, however, since they
/// represent both the x and y state it doesn't feel natural to place them there.
#[derive(Clone, Debug, PartialEq)]
pub enum Solution {
    /// No solution
    None,

    /// Single Solution
    Single((Form, Form)),

    /// Multiple Solutions
    Multiple(Box<[(Form, Form)]>),

    /// Infinite solutions represented as a recurrence: `X_{n+1} = A*X_{n} + Q*Y_{n}`
    // Since a recurrence relates both x and y, it can not be specified as an
    // individual form.
    Recurrence(RecurrenceEquation),

    /// A series of recurrence relations.
    Recurrences(Box<[RecurrenceEquation]>),
}

/// A specific form that a single solution may be. A form represents a range of
/// solutions which satisfy a given equation.
///
/// These are copied by default.
///
/// ## Examples
///
/// Consider the equation `2x - 4 = 0`.
///
/// Clearly we must have `x = 2` so x is of the form `Point(2)`.
///
/// Y is a free variable and so is of the form `All`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Form {
    /// `(A, B, C)` => Satisfied by the parabola `At² + Bt + C` for `t` in ℤ
    Quadratic(i64, i64, i64),

    /// `(A, B)` => Satisfied by the line `At + B` for `t` in ℤ
    Linear(i64, i64),

    /// `(A)` => Satisfied by the singular point `A`
    Point(i64),

    /// Satisfied by all points
    All
}

impl Form {
    /// Does this form produce multiple equations
    pub fn can_evaluate(&self) -> bool {
        match self {
            &Form::Quadratic(..) | &Form::Linear(..) => true,
            _ => false
        }
    }

    /// Evaluate the equation at the specified position t in Z
    pub fn evaluate(&self, t: i64) -> i64 {
        match self {
            &Form::Quadratic(a, b, c) => a * t * t + b * t + c,
            &Form::Linear(a, b) => a * t + b,
            x @ _ => panic!("cannot evaluate form of type '{:?}'", x)
        }
    }

    /// Retrieve the first solution from a form. This is mainly useful for
    /// reducing infinite solutions when we are only interested in one.
    pub fn get_first(&self) -> i64 {
        match self {
            // Default to evaluation at `t = 0` or the first recurrence
            &Form::Quadratic(_, _, c) => c,
            &Form::Linear(_, b) => b,
            &Form::Point(a) => a,
            &Form::All => 0,
        }
    }
}

impl Solution {
    /// Return the first pair of solutions if they exist.
    /// This will evaluate linear and quadratic solutions and retrieve the
    /// first solution in a recurrence.
    ///
    /// This is a convenience function for specific cases. Consider using a
    /// match clause and extracting only what you need if you know the
    /// particular form and output you will receive.
    pub fn get_first(&self) -> Option<(i64, i64)> {
        match self {
            &Solution::None => None,
            &Solution::Single((ref x, ref y)) => Some((x.get_first(), y.get_first())),

            &Solution::Multiple(ref sols) => {
                assert!(sols.len() > 0);
                let (x, y) = sols[0];
                Some((x.get_first(), y.get_first()))
            }

            &Solution::Recurrence(ref eq) => {
                eq.clone().nth(0)
            }

            &Solution::Recurrences(ref rec) => {
                assert!(rec.len() > 0);
                rec[0].clone().nth(0)
            }
        }
    }
}

/// Generic interface for solving diophantine equations.
pub fn solve(a: i64, b: i64, c: i64, d: i64, e: i64, f: i64) -> Solution {
    if a == 0 && b == 0 && c == 0 {
        solve_linear(d, e, f)
    }
    else if a == 0 && c == 0 && b != 0 {
        solve_simple_hyperbolic(b, d, e, f)
    }
    else if b*b - 4*a*c < 0 {
        solve_elliptical(a, b, c, d, e, f)
    }
    else if b*b - 4*a*c == 0 {
        solve_parabolic(a, b, c, d, e, f)
    }
    else {
        solve_hyperbolic(a, b, c, d, e, f)
    }
}

/// Map a set of solutions to the most appropriate type. Sometimes we do
/// not know in advance if we will have none, single, or multiple solutions.
fn map_solution_set(mut sols: Vec<(Form, Form)>) -> Solution {
    match sols.len() {
        0 => Solution::None,
        1 => {
            let (x, y) = sols.pop().unwrap();
            Solution::Single((x, y))
        }
        _ => Solution::Multiple(sols.into_boxed_slice())
    }
}

/// Solve a linear equation of the form `Dx + Ey + F = 0`.
///
/// # Returns
/// - `Solution::None`
/// - `Solution::Single`
pub fn solve_linear(d: i64, e: i64, f: i64) -> Solution {
    match (d, e, f) {
        (0, 0, f) => {
            // F = 0 => solution iff F == 0
            if f == 0 {
                Solution::Single((Form::All, Form::All))
            } else {
                Solution::None
            }
        }

        (0, e, f) => {
            // Ey + F = 0 => solution if f | e
            if f % e == 0 {
                Solution::Single((Form::All, Form::Point(-f / e)))
            } else {
                Solution::None
            }
        }

        (d, 0, f) => {
            // Dx + F = 0 => solution if f | d
            if f % d == 0 {
                Solution::Single((Form::Point(-f / d), Form::All))
            } else {
                Solution::None
            }
        }

        (d, e, f) => {
            let g = helper::gcd(d, e);

            // D | g and E | g => Dx + Ey | g for x, y in Z so if f !| g there
            // are no solutions
            if f % g != 0 {
                Solution::None
            }
            // f | g
            else {
                let d = d / g;
                let e = e / g;
                let f = f / g;

                // Match the resulting sign
                let (a, b) = helper::extended_gcd(d, e);

                // Solutions are actually of form `x = et +- fa` and
                // `y = -dt +- fb'` for `t` in `Z`.
                Solution::Multiple(
                    Box::new([
                        (Form::Linear(e, f*a), Form::Linear(-d, f*b)),
                        (Form::Linear(e, -f*a), Form::Linear(-d, -f*b)),
                    ])
                )
            }
        }
    }
}

/// Solve a quadratic equation of the form `Bxy + Dx + Ey + F = 0`
///
/// # Returns
/// - `Solution::None`
/// - `Solution::Single`
/// - `Solution::Multiple`
pub fn solve_simple_hyperbolic(b: i64, d: i64, e: i64, f: i64) -> Solution {
    // Lines parallel to x and y axes
    if d*e - b*f == 0 {
        // Bx + E = 0 => E | B
        if b % e == 0 {
            Solution::Single((Form::Point(-e / b), Form::All))
        }
        // By + D = 0 => E | B
        else if b % d == 0 {
            Solution::Single((Form::All, Form::Point(-d / b)))
        }
        else {
            Solution::None
        }
    }
    // Hyperbola parallel to x and y axes
    else {
        let mut sols = Vec::new();

        // Find the divisors of DE - BF
        for divisor in helper::divisors(d*e - b*f) {
            // check if an integer solution is found and add it if so
            if (divisor - e) % b == 0 && ((d*e - b*f) / divisor - d) % b == 0 {
                // Push both negative and positive divisor value forms
                sols.push((Form::Point((divisor - e) / b),
                           Form::Point(((d*e - b*f) / divisor - d) / b)));
                sols.push((Form::Point((-divisor - e) / b),
                           Form::Point(((d*e - b*f) / -divisor - d) / b)));
            }
        }

        map_solution_set(sols)
    }
}

/// Solve a quadratic equation where `B^2 - 4AC < 0`.
///
/// # Returns
/// - `Solution::None`
/// - `Solution::Single`
pub fn solve_elliptical(a: i64, b: i64, c: i64, d: i64, e: i64, f: i64) -> Solution {

    let (lo, hi) = helper::quadratic_roots(b*b - 4*a*c, 2*(b*e - 2*c*d), e*e - 4*c*f).unwrap();
    let mut sols = Vec::new();

    // Iterate through integer portions in [lo, hi], return the first solution currently.
    for x in (lo as i64)..(hi as i64) {
        println!("{}", x);
        let discriminant = (b*x + e) * (b*x + e) - 4*c*(a*x*x + d*x + f);

        if helper::is_perfect_sq(discriminant) {
            let discriminant_sqrt = (discriminant as f64).sqrt() as i64;

            let positive = -(b*x + e) + discriminant_sqrt;
            if positive % (2*c) == 0 {
                sols.push((Form::Point(x), Form::Point(positive / (2*c))));
            }
            else {
                // Only check the negative y value if the positive failed
                let negative = -(b*x + e) - discriminant_sqrt;
                if negative % (2*c) == 0 {
                    sols.push((Form::Point(x), Form::Point(negative / (2*c))));
                }
            }
        }
    }

    map_solution_set(sols)
}

/// Solve a quadratic equation where `B^2 - 4AC == 0`.
///
/// # Returns
/// - `Solution::None`
/// - `Solution::Recurrence`
pub fn solve_parabolic(a: i64, b: i64, c: i64, d: i64, e: i64, f: i64) -> Solution {
    // match the sign of g with a
    let g = helper::gcd(a, c) * if a > 0 { 1 } else { -1 };

    // a and c > 0 or one == 0
    let a = a / g;
    let c = c / g;

    // b^2 - 4ac = 0 => b^2 / 4 == ac => a and c are perfect squares
    let a_sqrt = (a as f64).sqrt() as i64;

    // sqrt of c is negative if b / a < 0
    let c_sqrt = if b / a < 0 { -1 } else { 1 } * (c as f64).sqrt() as i64;

    // sqrt(c)D - sqrt(a)E == 0 => parallel lines
    if c_sqrt * d - a_sqrt * e == 0 {
        // we have sqrt(a)x + sqrt(c)y - u_1 == 0 and sqrt(a)x + sqrt(c)y - u_2 == 0
        // where u_1 and u_2 are roots of the equation sqrt(a)gu^2 + Du + sqrt(a)F == 0
        let (u1, u2) = helper::quadratic_roots(a_sqrt*g, d, f).unwrap();

        // reduces to a linear problem
        solve_linear(a_sqrt, c_sqrt, u1 as i64)
    }
    // sqrt(c)D - sqrt(a)E != 0 => parabola
    else {
        let mut sols = Vec::new();

        // we can compute directly using a closed-form (which is rather verbose)
        // which gives us a quadratic of solutions.
        let target = c_sqrt * d - a_sqrt * e;

        for u in 0..(c_sqrt*d - a_sqrt*e).abs() {
            // Find all u .st `sqrt(a)*g*u^2 + d*u + sqrt(a)*f` | `sqrt(c)*d - sqrt(a)*e`
            let eval_poly = a_sqrt*g*u*u + d*u + a_sqrt*f;

            if eval_poly % target == 0 {
                // if u is a divisor then we can construct a closed form solution
                let den = c_sqrt * d - a_sqrt * e;
                sols.push(
                    // equation in reference are wrong by a sign it appears
                    (Form::Quadratic(c_sqrt * g * (a_sqrt * e - c_sqrt * d),
                                     (e + 2 * c_sqrt * g * u),
                                     -(c_sqrt * g * u * u + e * u + c_sqrt * f) / den),
                     Form::Quadratic(a_sqrt * g * (c_sqrt * d - a_sqrt * e),
                                     -(d + 2 * a_sqrt * g * u),
                                     (a_sqrt * g * u * u + d * u + a_sqrt * f) / den)
                    )
                );
            }
        }

        map_solution_set(sols)
    }
}

/// Solve a quadratic equation where `B^2 - 4AC > 0`.
pub fn solve_hyperbolic(a: i64, b: i64, c: i64, d: i64, e: i64, f: i64) -> Solution {
    match (a, b, c, d, e, f) {
        (a, b, c, 0, 0, 0) => {
            let mut sols = vec![(Form::Point(0), Form::Point(0))];

            // determine any other solutions
            if helper::is_perfect_sq(b*b - 4*a*c) {
                // reduces into solving two linear equations
                let ks = ((b*b - 4*a*c) as f64).sqrt() as i64;

                // `2ax + (b - sqrt(b^2 - 4ac))y = 0`
                match solve_linear(2*a, b - ks, 0) {
                    Solution::Single((x, y)) => sols.push((x, y)),
                    _ => ()
                }

                // `2ax - (b - sqrt(b^2 - 4ac))y = 0`
                match solve_linear(2*a, b + ks, 0) {
                    Solution::Single((x, y)) => sols.push((x, y)),
                    _ => ()
                }
            }

            map_solution_set(sols)
        }

        (a, b, c, 0, 0, f) if helper::is_perfect_sq(b*b - 4*a*c) => {
            let mut sols = vec![(Form::Point(0), Form::Point(0))];
            let k = ((b*b - 4*a*c) as f64).sqrt() as i64;

            for divisor in helper::divisors(-4*a*f) {
                let y = divisor + 4*a*f / divisor;
                let x = divisor - (b + k) * y;

                if x % (2*a) == 0 && y % (2*k) == 0 {
                    sols.push((Form::Point(x / 2*a), Form::Point(y / 2*k)));
                }
            }

            map_solution_set(sols)
        }

        (a, b, c, 0, 0, f) if !helper::is_perfect_sq(b*b - 4*a*c) => {
            let g = gcd!(a, b, c);

            if f % g != 0 {
                Solution::None
            }
            else if 4*f*f >= b*b - 4*a*c {
                // reduce equation by gcd
                let a = a / g;
                let b = b / g;
                let c = c / g;
                let f = f / g;

                // as^2 + bs + c must be a multiple of f
                for s in 0..f.abs()-1 {
                    if (a*s*s + b*s + c) % f != 0 {
                        continue;
                    }

                    // substitute s for `[-(as^2 + bs + c)/f]y^2 + (2sa + b)yz - afz^2 = 1`
                    let ac = -(a*s*s + b*s + c)/f;
                    let bc = 2*s*a + b;
                    let cc = -a*f;

                    // determine continued fraction expansion of the roots of
                    // the specified equation.
                    let leading = -bc;
                    let root = bc*bc - 4*ac*cc;
                    let den = 2*ac;

                    let cfrac = helper::continued_fraction(leading, root, den);

                    // we need to iterate through 2 cycles of the periodic part,
                    // scanning the convergents for equality to 1
                    // TODO: may require a big number library for accuracy
                }

                Solution::None
            }
            else {
                Solution::None
            }
        }

        // general quadratic equation
        (a, b, c, d, e, f) => {
            let g = helper::gcd(4*a*c - b*b, 2*a*e - b*d);
            Solution::None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_01() {
        let solution = solve_linear(10, 84, 16);
        assert_eq!(solution, Solution::Multiple(Box::new([
                                (Form::Linear(42, -136), Form::Linear(-5, 16)),
                                (Form::Linear(42, 136), Form::Linear(-5, -16)),
                            ])));
    }

    #[test]
    fn linear_02() {
        let solution = solve_linear(3, 0, 6);
        assert_eq!(solution, Solution::Single((Form::Point(-2), Form::All)));
    }

    #[test]
    fn linear_03() {
        let solution = solve_linear(5, 0, 3);
        assert_eq!(solution, Solution::None);
    }

    #[test]
    fn simple_hyperbolic_01() {
        let solution = solve_simple_hyperbolic(2, 5, 56, 7);
        assert_eq!(solution,
                   Solution::Multiple(Box::new([
                       (Form::Point(-27), Form::Point(64)),
                       (Form::Point(-29), Form::Point(-69)),
                       (Form::Point(-21), Form::Point(7)),
                       (Form::Point(-35), Form::Point(-12)),
                       (Form::Point(-9), Form::Point(1)),
                       (Form::Point(-47), Form::Point(-6)),
                       (Form::Point(105), Form::Point(-2)),
                       (Form::Point(-161), Form::Point(-3)),
                   ]))
        );
    }

    #[test]
    fn elliptical_01() {
        let solution = solve_elliptical(42, 8, 15, 23, 17, -4915);
        assert_eq!(solution, Solution::Single((Form::Point(-11), Form::Point(-1))));
    }

    #[test]
    fn parabolic_01() {
        let solution = solve_parabolic(8, -24, 18, 5, 7, 16);
        assert_eq!(solution, Solution::Multiple(Box::new([
                                 (Form::Quadratic(-174, -17, -2), Form::Quadratic(-116, -21, -2)),
                                 (Form::Quadratic(-174, -41, -4), Form::Quadratic(-116, -37, -4))
                             ]))
        );
    }
}

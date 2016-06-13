//! Contains helper functions for use in `solve`.

#![allow(dead_code)]

/// Compute the gcd of two numbers
pub fn gcd(mut x: i64, mut y: i64) -> i64 {
    assert!(x >= 0 && y >= 0);

    if x == 0 {
        return y;
    } else if y == 0 {
        return x;
    }

    let mut result = 0;

    while x & 1 == 0 && y & 1 == 0 {
        x >>= 1;
        y >>= 1;
        result += 1;
    }

    while x != y {
        if x & 1 == 0 {
            x >>= 1;
        } else if y & 1 == 0 {
            y >>= 1;
        } else if x > y {
            x = (x - y) >> 1;
        } else {
            y = (y - x) >> 1;
        }
    }

    x * (1 << result)
}

macro_rules! gcd {
    ($($x:expr),*) => {
        {
            let mut result = 0; // gcd(0, x) == x
            $(
                result = ::helper::gcd(result, $x);
            )*
            result
        }
    };
}

/// Compute the lcm of two numbers
pub fn lcm(x: i64, y: i64) -> i64 {
    x * y / gcd(x, y)
}

macro_rules! lcm {
    ($($x:expr),*) => {
        {
            let mut result = 1; // lcm(1, x) == x
            $(
                result = ::helper::lcm(result, $x);
            )*
            result
        }
    };
}


/// Compute the extended gcd of two numbers.
///
/// In other words solve a*x + b*y = +-gcd(x, y) for the minimum a, b which
/// satisfy this constraint.
pub fn extended_gcd(x: i64, y: i64) -> (i64, i64) {
    let mut a = 0;
    let mut b = 1;
    let mut m = 1;
    let mut n = 0;
    let mut p = y;
    let mut q = x;

    while p != 0 {
        let r = q / p;

        let z1 = p;
        p = q - r * p;
        q = z1;

        let z2 = m;
        m = n - r * m;
        n = z2;

        let z3 = a;
        a = b - r * a;
        b = z3;
    }

    println!("{} {}: {} {}: {} {}", a, b, m, n, p, q);
    println!("{} {}", b, n);

    // Always return (-, +) parity if given the choice. This is what is
    // usually expected in the solve phase.
    if b >= 0 && a < 0 {
        (-b, -n)
    } else {
        (b, n)
    }
}

/// Return whether the specified integer is a perfect square
pub fn is_perfect_sq(n: i64) -> bool {
    if n & 2 == 2 || n & 7 == 5 || n & 11 == 8 || n & 31 == 20 ||
       n & 47 == 32 || n & 127 == 80 || n & 191 == 128 || n & 511 == 320 {
        false
    } else {
        let v = (n as f64).sqrt().trunc() as i64;
        v * v == n
    }
}

/// Return the roots of a quadratic equation of the form `ax^2 + bx + c = 0`.
/// Return None if no real solutions exist. If only one solution exists, then
/// it will be returned twice.
///
/// The tuple returned (a, b) will be ordered such that a <= b
pub fn quadratic_roots(a: i64, b: i64, c: i64) -> Option<(f64, f64)> {
    let discriminant = (b*b - 4*a*c) as f64;
    if discriminant < 0_f64 {
        None
    } else {
        Some(((-b as f64 + discriminant.sqrt()) / (2*a) as f64,
              (-b as f64 - discriminant.sqrt()) / (2*a) as f64))
    }
}

/// Return the continued fraction expansion of a value of the form
/// `(a + sqrt(b)) / c`. The returned value is of the form:
///
/// `(integer-part, non-periodic-part, periodic-part)`
///
/// # Note
/// b must not be a perfect square
///
/// The returned fraction may be infinitely long (non-zero periodic-part) or
/// will have a finite length (zero-length periodic-part).
#[allow(unused_assignments)]
pub fn continued_fraction(a: i64, b: i64, c: i64) -> (i64, Box<[i64]>, Box<[i64]>) {
    use std::collections::HashMap;

    // input b cannot be a perfect square
    assert!(!is_perfect_sq(b));

    let mut a0 = (((b as f64).sqrt() + a as f64) / c as f64).floor() as i64;
    let mut v0 = c;
    let mut u0 = -a + a0 * c;

    let integer = a0;
    let mut components = Vec::new();

    // we can match the start of the periodic part by hashing against a tuple
    // (a0, u0, v0) against the start index it was encountered.
    let mut sequence = HashMap::<(i64, i64, i64), usize>::new();

    let mut pidx = 0;

    loop {
        // An provides the fractional component
        v0 = (b - u0 * u0) / v0;
        a0 = (((b as f64).sqrt() + u0 as f64) / (v0 as f64)).floor() as i64;
        u0 = a0 * v0 - u0;

        let k = (a0, u0, v0);

        // determine if we have encountered this sequence before
        if sequence.contains_key(&k) {
            pidx = *sequence.get(&k).unwrap();
            break;
        } else {
            sequence.insert(k, components.len());
        }

        components.push(a0);
    }

    let periodic = components.split_off(pidx);
    (integer, components.into_boxed_slice(), periodic.into_boxed_slice())
}

pub struct DivisorIterator {
    /// Current divisor we are checking
    divisor: i64,

    /// Current value we are reducing
    value: i64,
}

impl Iterator for DivisorIterator {
    type Item = i64;

    fn next(&mut self) -> Option<i64> {
        if self.divisor == 1 {
            self.divisor += 1;
            return Some(1)
        }

        loop {
            if self.divisor > self.value {
                break;
            }

            if self.value % self.divisor == 0 {
                self.divisor += 1;
                return Some(self.divisor - 1);
            }
            else {
                self.divisor += 1;
            }
        }

        None
    }
}

/// Return an iterator over all the divisors of the specified number.
///
/// This is a naive algorithm for the moment, we can optimize this
/// specifically later if required.
pub fn divisors(x: i64) -> DivisorIterator {
    DivisorIterator {
        divisor: 1,
        value: x,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {
            {
                let (a, b) = (&$a, &$b);
                assert!((*a - *b).abs() < 1.0e-6,
                        "{} is not approx equal to {}", *a, *b);
            }
        };
    }

    #[test]
    fn gcd_function_01() {
        assert_eq!(gcd(64, 2348), 4);
        assert_eq!(gcd(64, 2347), 1);
        assert_eq!(gcd(0, 56), 56);
        assert_eq!(gcd(0, 0), 0);
    }

    #[test]
    #[should_panic]
    fn gcd_function_invalid_01() {
        gcd(-5, 54);
    }

    #[test]
    fn gcd_macro_01() {
        assert_eq!(gcd!(0, 64, 2348), 4);
        assert_eq!(gcd!(1, 2, 3, 4, 5, 6, 7, 8), 1);
    }

    #[test]
    fn divisors_01() {
        assert_eq!(divisors(266).collect::<Vec<_>>(), &[1, 2, 7, 14, 19, 38, 133, 266]);
        assert_eq!(divisors(36).collect::<Vec<_>>(), &[1, 2, 3, 4, 6, 9, 12, 18, 36]);
    }

    #[test]
    fn is_perfect_sq_01() {
        assert_eq!(is_perfect_sq(5329), true);
        assert_eq!(is_perfect_sq(5328), false);
        assert_eq!(is_perfect_sq(1532886657604), true);
        assert_eq!(is_perfect_sq(1532886657606), false);
        assert_eq!(is_perfect_sq(1), true);
    }

    #[test]
    fn quadratic_roots_01() {
        let (x, y) = quadratic_roots(1, 9, -7).unwrap();
        assert_approx_eq!(x, 0.7201533);
        assert_approx_eq!(y, -9.7201533);
        assert_eq!(quadratic_roots(-5, -9, -7), None);
    }

    #[test]
    fn ext_gcd_01() {
        assert_eq!(extended_gcd(77, -89), (37, 32));
        assert_eq!(extended_gcd(2, 3), (-1, 1))
    }

    #[test]
    fn continued_fraction_01() {
        let solution = continued_fraction(-725, 313, 608);
        assert_eq!(solution,
                   (-2,
                    vec![1, 5].into_boxed_slice(),
                    vec![8, 5, 1, 3, 1, 1, 2, 2, 1, 1, 3, 1, 5, 8, 1, 2, 17, 2, 1]
                    .into_boxed_slice())
        );
    }
}

//! A crate for solving quadratic diophantine equations.
//!
//! See the [solve](./solve/index.html) module for the main API.
//!
//! ## Examples
//!
//! ```
//! extern crate diophantine;
//!
//! use diophantine::Solution;
//!
//! fn main() {
//!     // 0x^2 + 0xy + 0y^2 + 0x + 6y + 12 = 0
//!     let sol1 = diophantine::solve(0, 0, 0, 0, 6, 12);
//!
//!     match sol1 {
//!         // single solution
//!         Solution::Single(x) => {
//!             println!("{:?}", x);
//!         }
//!
//!         // finite set of solutions
//!         Solution::Multiple(x) => {
//!             println!("{:?}", x);
//!         }
//!
//!         // an iterator over an infinite sequence of solutions
//!         Solution::Recurrence(eq) => {
//!             println!("{:?}", eq.take(5).collect::<Vec<_>>());
//!         }
//!
//!         // multiple recurrences can be found
//!         Solution::Recurrences(eq) => {
//!             println!("{:?}", eq);
//!         }
//!
//!         // always check the case where no solution exists
//!         Solution::None => {
//!             println!("no solution found");
//!         }
//!     }
//! }
//! ```

// Contains things like the gcd, euclid's algorithm etc.
#[macro_use]
mod helper;

// Contains code which can evaluate recurrence and store them in a consistent
// form.
mod recurrence;

// We expose the iterator over recurrence equations for documentation.
pub use recurrence::RecurrenceEquation;

// Contains actual solver code.
pub mod solve;

// Re-export all solve functions so the user only has to prefix with
// `diophantine`
pub use solve::*;

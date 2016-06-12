//! A crate for solving quadratic diophantine equations.
//!
//! See the [solve](./solve/index.html) module for the main API.
//!
//! ## Examples
//!
//! ```text
//! extern crate diophantine;
//!
//! use diophantine::{Solution, Form};
//!
//! fn main() {
//!     // 5x^2 + 0xy + 1y^2 + 0x + 2y + 1 = 0
//!     let solution = diophantine::solve(5, 0, 1, 0, 2, 1);
//!
//!     // Pattern match for the specific type of equation
//!     match diophantine::solve(0, 0, 0, 1, 1, 1) {
//!         // Some equations can return multiple solutions, these can be
//!         // one of many forms.
//!         Solution::Single(x) | Solution::Multiple(x) => {
//!             println!("{:?}", x)
//!         }
//!
//!         // Recurrences can be returned which produce an infinite
//!         // sequence of values.
//!         Solution::Recurrence(equation) => {
//!             // Print the first 5 solutions of this recurrence
//!             println!("{:?}", eq.take(5).collect::<Vec<_>>())
//!         }
//!
//!         // No solutions may be found for a form
//!         Solution::None => {
//!             println!("no solution found")
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

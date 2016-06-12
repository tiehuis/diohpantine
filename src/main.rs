extern crate diophantine;

fn main() {
    // Check extended_gcd calculation
    println!("{:#?}", diophantine::solve_simple_hyperbolic(2, 5, 56, 7));
    // println!("{:?}", diophantine::solve_elliptical(42, 8, 15, 23, 17, -4915));
}

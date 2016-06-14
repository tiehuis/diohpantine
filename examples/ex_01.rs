extern crate diophantine;

use diophantine::Solution;

fn main() {
    // 0x^2 + 0xy + 0y^2 + 0x + 6y + 12 = 0
    let sol1 = diophantine::solve(0, 0, 0, 0, 6, 12);

    match sol1 {
        // single solution
        Solution::Single(x) => {
            println!("{:?}", x);
        }

        // finite set of solutions
        Solution::Multiple(x) => {
            println!("{:?}", x);
        }

        // an iterator over an infinite sequence of solutions
        Solution::Recurrence(eq) => {
            println!("{:?}", eq.take(5).collect::<Vec<_>>());
        }

        // multiple recurrences can be found
        Solution::Recurrences(eq) => {
            println!("{:?}", eq);
        }

        // always check the case where no solution exists
        Solution::None => {
            println!("no solution found");
        }
    }

    // 5x + 7y = 39
    let sol2 = diophantine::solve_linear(5, 7, 39);

    match sol2 {
        Solution::Single((x, y)) => {
            println!("({:?}, {:?})", x, y);
        }

        Solution::Multiple(x) => {
            println!("{:?}", x);
        }

        // The documentation states that `solve_linear` will only return the
        // variants of type `Solution::Single` and `Solution::Multiple` so we
        // can discard any other types.
        _ => unreachable!()
    }
}

# Diophantine

Provides methods for solving quadratic and linear diophantine equations.

Diophantine equations are of equations of the form:

```
Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
```

The methods used to solve these are those specified [here](https://www.alpertron.com.ar/METHODS.HTM).


## Examples

```
extern crate diophantine;

use diophantine::{Solution, Form};

fn main() {
    // 5x^2 + 0xy + 1y^2 + 0x + 2y + 1 = 0
    let solution = diophantine::solve(5, 0, 1, 0, 2, 1);

    // Pattern match for the specific type of equation
    match diophantine::solve(0, 0, 0, 1, 1, 1) {
        // Some equations can return multiple solutions, these can be
        // one of many forms.
        Solution::Single(x) | Solution::Multiple(x) => {
            println!("{:?}", x)
        }

        // Recurrences can be returned which produce an infinite
        // sequence of values.
        Solution::Recurrence(equation) => {
            // Print the first 5 solutions of this recurrence
            println!("{:?}", eq.take(5).collect::<Vec<_>>())
        }

        // No solutions may be found for a form
        Solution::None => {
            println!("no solution found")
        }
    }
}
```

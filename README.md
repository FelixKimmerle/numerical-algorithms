# Numerical Algorithms

C++ implementations of algorithms from the course **Numerische Grundlagen** (Numerical Foundations) at the
**University of Stuttgart**.

## Educational focus

This project is purely for learning and exam preparation. The algorithms are implemented the way you compute them by
hand in the exam, not the way a numerics library would do it. Every algorithm prints its intermediate steps: pivot
choices, elimination factors, Givens rotations, the full Aitken-Neville and divided-difference schemes, the Horner-like
evaluation and the Fourier matrix. You can use the output to check your own calculations step by step.

So the code favors readability and matching the lecture over speed and numerical robustness. Matrices are plain
`std::vector`s, polynomials are multiplied out explicitly and nothing is optimized. Section and algorithm numbers in the
code and output refer to the lecture notes (Skript) of the summer semester 2022.

## Implemented algorithms

| Topic | Algorithm | Lecture notes | Function |
| --- | --- | --- | --- |
| Linear systems | LR decomposition with column pivoting (`P·A = L·R`) | Alg. 3.10 | `Matrix::LR` |
| Linear systems | LR decomposition, row-by-row variant | Alg. 3.8 | `Matrix::LR_Line` |
| Linear systems | QR decomposition with Givens rotations (also non-square) | Alg. 3.21 | `Matrix::QR` |
| Interpolation | Lagrange polynomials | 4.3.3 | `lagrange_interpolation` |
| Interpolation | Aitken-Neville scheme | Alg. 4.12 | `aitken_neville` |
| Interpolation | Newton's divided differences | Alg. 4.17 | `newton` |
| Interpolation | Newton representation in the monomial basis | (4.44) | `get_newton_polynom` |
| Interpolation | Evaluating the Newton representation | Alg. 4.16 | `evaluate_newton` |
| Trigonometric interpolation | Inverse discrete Fourier transform (`v = M·c`) | Def. 4.30 | `IDFT` |
| Trigonometric interpolation | Discrete Fourier transform (via the IDFT) | Def. 4.30, (4.151) | `DFT` |

Not implemented yet: linear least squares (3.5), Hermite interpolation (4.3.6), the fast Fourier transform (4.5.2).

## Building and running

Only a C++17 compiler and `make` are needed, there are no dependencies.

```bash
make run
```

This builds the release version to `build/release/num` and runs it. `make debug` builds a debug version to
`build/debug/num`.

`main.cpp` runs an example for every algorithm with exam-style inputs that give clean numbers, and checks each result
(e.g. `P·A - L·R = 0`, `QᵀQ = I`, the interpolation polynomial hits all points, `IDFT(DFT(v)) = v`). To work on your own
exercise, change the inputs in [src/main.cpp](src/main.cpp).

## Reading the output

- Matrix rows and columns are numbered from 1 like in the lecture (`l_21`, `G_12`, step `j = 1`), interpolation nodes
  from 0 (`x_0, ..., x_n`, `P_i,k`).
- Numbers are printed with 3 decimal places (`PRECISION` in [src/Format.hpp](src/Format.hpp)); rounding noise like
  `1e-17` is printed as `0`.
- In the Aitken-Neville and divided-difference schemes, the newly computed entry is **red**, the two entries it was
  computed from are **green** and the two nodes that were used are highlighted.

Excerpt of the divided differences for the points `(0, 1), (1, 0), (2, 5), (3, 22)`:

```
[x_0,...,x_2]f = (5.000 - (-1.000)) / (2.000 - 0.000) = 3.000
   0.000 |  1.000 |
   1.000 |  0.000 | -1.000 |
   2.000 |  5.000 |  5.000 | 3.000 |
   3.000 | 22.000 | 17.000 |

...

Newton representation (4.44):
p(x) = 1.000 - 1.000 * x + 3.000 * x * (x - 1.000) + 1.000 * x * (x - 1.000) * (x - 2.000)
     = x^3 - 2.000x + 1.000
```

## Project structure

```
src/main.cpp            examples for every algorithm
src/Matrix.*            matrix class with the LR and QR decompositions
src/Interpolation.*     polynomial interpolation and the (inverse) discrete Fourier transform
src/Polynom.*           polynomials in the monomial basis
src/Format.hpp          number formatting for the step output
src/Colors.hpp          ANSI color codes
```

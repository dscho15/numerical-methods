# Lecture Notes

- Read from page `1-52`
- It requires some Latex interpreter to read the math.
- Compile with: `pandoc test.md -o test.pdf`

---

- Data types can differ in the number of bits utilized (wordlength), which differs from the fundamental respect whether it is stored in ``fixed-point`` or ``floating-point``
- Floating point is represented> by: 
  $$ S \times M \times b^{E-e}$$
b is the `base` (binary: 2), e is the bias of `exponent` (fixed integer constant for a machine), S depends on the sign, M (23 for a float) and E depends upon the number (1-254 for a float).
- All modern processers share the same floating-point representation, namely `IEEE Standard 754-1985`.
- Big-endian is reverse group of bytes.
- `` std::numerical_limits<double>::min `` is used to get the minimum number for a given datatype within STL
- `` std::numerical_limits<double>::epsilon`` is used to get the error between addition of two double, and can be generalized.
- Different type of error 
  - Roundoff Error (addition, multiplication of floats etc, `hardware`)
  - Truncation Error (integral of some function, convert from float to int, `programmer`)
  - Stability (unstable reccurence)


- ``nr3.h`` has its own typedef, to avoid penalizaion a run-time.
  - Format for vectors: `VectInt, VecUint,VecChar,VecDoub`
  - Format for matrix: `MatInt...`
  - `_I` is defined as input, const.

---

# Solution of Linear Algebraic Equations

Given a set of linear algebraic equations 

$$ a_{00}x_0 + a_{01}x_1 +  a_{03}x_2 +  a_{04}x_3 = b_0$$
$$ a_{10}x_0 + a_{11}x_1 +  a_{13}x_2 +  a_{14}x_3 = b_1$$

Here `N` is unknowns $x_j$, j = 0,1,2,3, this can be written as and $b$ can be written likewise. If the amount of unknowns from matrix A is equal to N then there is a chance of finding a unique solution of x. `M` is the amount of rows.

$$ A \cdot x = b$$

The dot is equal to matrix multiplication and is a so called contraction operator, that represents the sum over a pair of indicies for example

$$C = A \cdot B \xrightarrow{} c_{ik} = \sum_j a_{ij}b_{jk}$$ 

## Nonsingular versus singular

Even though M = N there might not be a unique solution this is called a singularity, `row degeneracy` (række degeneration, en række er en linær kombination af de andre) or `column degeneracy` (kolonne degeneration, en kolonne er en linær kombination).

- Prevention of these two things are important
  - Some equations might be close linearly depedent that roundoff errors in the machine render them depedent
  - `Accumulated roundoff errors` can swamp the true solution (if `N` is large).

## Tasks of Computational Linear Algebra
- When M = N
  - Solution of matrix equation 
  - Solution of more than one matrix equation
  - Calculation of inverse matrix.
  - Calculation of determinant of square matrix.

If `M < N` or the same size, then there is effectively fewer unknowns than knowns. In this case there is usually no solution or else more than one solution. 

- The solution space consist of a particular solution denoted **$x_p$**  added to any linear combination of $N-M$ vectors (which are said to be in the nullspace of the matrix A).

If there are more equations than unknowns, M > N, there is in general no solution $x$ to equation $A \cdot x = b$.

## What is a subspace?

A subspace $S$ is defined by the 3 rules:

- $\vec{O} \; \epsilon \; S$ is just the null vector.
- $\vec{v}_1, \vec{v}_2 \; \epsilon \; S \xrightarrow{} \; \vec{v}_1 + \vec{v}_2 \; \epsilon \; S$ two vectors added should result in a new vector in the space.
- $c \; \epsilon \; \mathbb{R} \; , \; \vec{v_1} \; \epsilon \; S \xrightarrow{}c\vec{v}_1 \; \epsilon \; S$ obvious linear relationship.

A subspace, the `nullspace`, could be defined for $A\cdot x = b$, counterintuitively we call it the null space of A. 
Example say $A \; \epsilon \; \mathbb{R}^{m\times n}$ and $x \; \epsilon \; \mathbb{R}^{n\times1}$ and if m > n, then there might a null space. ``Gaussian elimination can be used to derive the nullspace``.



## Gauss-Jordan elimination


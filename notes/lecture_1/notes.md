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

---

## Nonsingular versus singular

---

Even though M = N there might not be a unique solution this is called a singularity, `row degeneracy` (række degeneration, en række er en linær kombination af de andre) or `column degeneracy` (kolonne degeneration, en kolonne er en linær kombination).

- Prevention of these two things are important
  - Some equations might be close linearly depedent that roundoff errors in the machine render them depedent
  - `Accumulated roundoff errors` can swamp the true solution (if `N` is large).

---

## Tasks of Computational Linear Algebra

---
- When M = N
  - Solution of matrix equation 
  - Solution of more than one matrix equation
  - Calculation of inverse matrix.
  - Calculation of determinant of square matrix.

If `M < N` or the same size, then there is effectively fewer unknowns than knowns. In this case there is usually no solution or else more than one solution. 

- The solution space consist of a particular solution denoted **$x_p$**  added to any linear combination of $N-M$ vectors (which are said to be in the nullspace of the matrix A).

If there are more equations than unknowns, M > N, there is in general no solution $x$ to equation $A \cdot x = b$.

---

## What is a subspace?

---

A subspace $S$ is defined by the 3 rules:

- $\vec{O} \; \epsilon \; S$ is just the null vector.
- $\vec{v}_1, \vec{v}_2 \; \epsilon \; S \xrightarrow{} \; \vec{v}_1 + \vec{v}_2 \; \epsilon \; S$ two vectors added should result in a new vector in the space.
- $c \; \epsilon \; \mathbb{R} \; , \; \vec{v_1} \; \epsilon \; S \xrightarrow{}c\vec{v}_1 \; \epsilon \; S$ obvious linear relationship.

A subspace, the `nullspace`, could be defined for $A\cdot x = b$, counterintuitively we call it the null space of A. 
Example say $A \; \epsilon \; \mathbb{R}^{m\times n}$ and $x \; \epsilon \; \mathbb{R}^{n\times1}$ and if m > n, then there might a null space. ``Gaussian elimination can be used to derive the nullspace``.

---

## Gauss-Jordan elimination

---

- The idea is to add or subtract a linear combination of the given equation until equation contains only one of the unknowns.
- The deficienies:
  - It requires all the right hand sides to be stored and manipulated at the same time.
  - It is slow, when the inverse matrix is not desired, it is `considered stable` however.

The principle can be applied to multiple systems, e.g.

$$ [\mathbf{A}] \cdot [\mathbf{x}_0 \; | \; \mathbf{x}_1 \; | \; \mathbf{x}_2 \; | \; Y \; ] = [\mathbf{b}_0 \; | \; \mathbf{b}_1 \; | \; \mathbf{b}_2 \; | \; I ]$$ 

notice | is the column augmentation operator and just shows that it is possible to seperate the matrixes and the solve the linear sets.

$$ A \cdot \mathbf{x}_0 = \mathbf{b}_0$$

There is exactly 3 operations that is possible

- Swapping two rows,
- Multiplying a row by a nonzero number,
- Adding a multiple of one row to another row.

The main idea is to reach `row echoleon form`, an matrix with values equal to 0 under the diagonal. This is perfomed by taking $a_{1,1}$ and dividing every row underneath e.g. $a_{2,1} - \lambda a_{1,1} = 0$ and then continueing like this.

The main reason why we want to avoid this is because of $a_{1,1}$ is 0 hence $\lambda$ would be infinite.

---

## Gaussian Elimination: Partial Pivot Method

---

Just gonna use the link for this one. Can be used on example 3x3 | 3x3 etc

https://ocw.mit.edu/courses/chemical-engineering/10-34-numerical-methods-applied-to-chemical-engineering-fall-2005/lecture-notes/lecturenotes123.pdf

**The main idea**


    Find the entry in the left column with the largest absolute value. This entry is called the pivot.

    a_{ij} = max_{k} | a_kj | (if all elements are 0 then we give up.)

    Perform row interchange (if necessary), so that the pivot is in the first row.

    Find new pivot (repeat)

    ...

    When finished perform back-substituion (the process of computing the unknowns from a system that is in upper-triangular form)

  ``Pivoting helps reduce rounding errors; you are less likely to add/subtract with very small number (or very large) numbers``

  - Why pick the largest pivot? It minimizes the round-off error.
  - Partial pivoting is the interchanging of rows
  - Full pivoting is the interchanging of both rows and columns in order to place a particularly "good" element in the diagonal position prior to a particular operation. 

---


## Other Methods?

---

Column operations corresponds to multiplying wth simple matrices called C and its inverse, they are only used for doing permutations (ændre rækkefølgen).

Visual:

$$ (... R_2 \cdot R_1 \cdot R_0 \cdot A ) \cdot x = (... R_2 \cdot R_1 \cdot R_0 \cdot A ) b $$

``More to be added...``

---

# Gaussian Elimination with Backsubstitution

After the upper triangular matrix is obtained from gaussian elimination, then backsubstituion is begun. 

The algorithm for this is fairly straight forward

$$ x_3 = b_3'/a_{33} $$

Then $x_3$ is used for $x_2$

$$ x_2 = \frac{1}{a_{22}} [ b_2'-x_3a_{23}']$$

This proceeds, and is a recursive algorithm

$$ x_i = \frac{1}{a_{ii}}[b_2' - \sum^{N-1}_{j=i+1}a'_{ij}x_j] $$

This algorithm is faster than gaussian jordan hence we do not find the full algorithm.

##
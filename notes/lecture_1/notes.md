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

## LU-Decomposition Method

The idea is to split matrix A, which has NxN dimension to two matrices: L and U:

$$ \mathbf{A} \cdot \mathbf{x} = \mathbf{b}$$

L always have ones in the diagonal is a lower triangular matrix and U is a upper triangular matrix. We only use 4x4 as examples

$$
\begin{bmatrix}
a_{0,0} & a_{0,1} & a_{0,2} & a_{0,3} \\
a_{1,0} & a_{1,1} & a_{1,2} & a_{1,3} \\
a_{2,0} & a_{2,1} & a_{2,2} & a_{2,3} \\
a_{3,0} & a_{3,1} & a_{3,2} & a_{3,3}
\end{bmatrix} = \begin{bmatrix}
1 & 0 & 0 & 0 \\
L_{1,0} & 1 & 0 & 0 \\
L_{2,0} & L_{2,1} & 1 & 0 \\
L_{3,0} & L_{3,1} & L_{3,2} & 1
\end{bmatrix} \cdot \begin{bmatrix}
U_{0,0} & U_{0,1} & U_{0,2} & U_{0,3} \\
0 & U_{1,1} & U_{1,2} & U_{1,3} \\
0 & 0 & U_{2,2} & U_{2,3} \\
0 & 0 & 0 & U_{3,3}
\end{bmatrix} 
$$

---

The LU algorithm can be achieved by running through each column, followed by row using iterator *i* and for each row having two mechanisms depending on *j*. An example is shown using 3x3.

---


$$
\begin{bmatrix}
a_{0,0} & a_{0,1} & a_{0,2}  \\
a_{1,0} & a_{1,1} & a_{1,2}  \\
a_{2,0} & a_{2,1} & a_{2,2}  \\
\end{bmatrix} = \begin{bmatrix}
1 & 0 & 0 \\
L_{1,0} & 1 & 0  \\
L_{2,0} & L_{2,1} & 1  
\end{bmatrix} \cdot \begin{bmatrix}
U_{0,0} & U_{0,1} & U_{0,2}  \\
0 & U_{1,1} & U_{1,2}  \\
0 & 0 & U_{2,2}  \\
\end{bmatrix} 
$$

---

$$
\begin{bmatrix}
a_{0,0} & a_{0,1} & a_{0,2}  \\
a_{1,0} & a_{1,1} & a_{1,2}  \\
a_{2,0} & a_{2,1} & a_{2,2}  \\
\end{bmatrix} = 
\begin{bmatrix}
U_{0,0} & \cdots &  \\
L_{1,0}\cdot U_{0,0} & \cdots &   \\
L_{2,0}\cdot U_{0,0} & \cdots &   
\end{bmatrix}
$$

---

This is first column iteration

$$
\begin{bmatrix}
a_{0,0} & a_{0,1} & a_{0,2}  \\
a_{1,0} & a_{1,1} & a_{1,2}  \\
a_{2,0} & a_{2,1} & a_{2,2}  \\
\end{bmatrix} = 
\begin{bmatrix}
U_{0,0} & U_{0,1}  & \cdots  \\
L_{1,0}\cdot U_{0,0} & L_{1,0} \cdot U_{0,1} + U_{1,1} & \cdots  \\
L_{2,0}\cdot U_{0,0} & L_{2,0} \cdot U_{0,1} + L_{2,1}U_{1,1} &   \cdots
\end{bmatrix}
$$

---

Now let's try to generalize this using null indexing

    for j = 0,1 ... N-1

      for i = 0,1 ... j
        
        U[i][j] = A[i][j]

        for k = 0,1 ... j-1

          U[i][j] -= L[i][j] * U[i][j]

        end

      end

      for i = j+1, j+2 ... N-1

        U[i][j] = A[i][j]

        for k = j+1, j+1 ... N-1

          U[i][j] -= L[i][j] * U[i][j]

        end

        U[i][j] /= U[j][j]

      end

    end

---

Where $j$ runs from $0$ ... $N-1$

Where $i$ runs from $0$ ... $j$

$$ U_{ij} = a_{ij} - \sum_{k=0}^{j-1} L_{i,k}U_{k,j} $$

and where $i$ runs from $j+1$ ... $N$

$$ L_{ij} = \frac{a_{ij} - \sum_{k=0}^{j-1} L_{i,k}U_{k,j}}{U_{j,j}} $$

---

This is just decomposition and to solve a set of linear equations this has to be performed firstly. Then it is possible to calculate the set of equations. This is done by first analyzing the linear equation

$$ \mathbf{A} \cdot \mathbf{x} = \mathbf{b} \xrightarrow{} ( \mathbf{L} \cdot \mathbf{U})\cdot \mathbf{x} = \mathbf{b} $$
$$ \mathbf{L} \cdot ( \mathbf{U}\cdot \mathbf{x} ) = \mathbf{b} \xrightarrow{} \mathbf{L} \cdot \mathbf{y} = \mathbf{b} $$
$$ \mathbf{U} \cdot \mathbf{x} = \mathbf{y} $$

This is a recursive method and be can done for multiple e.g. $C=L\cdot U\cdot L' \cdot U'$. So let's figure out how to solve this, and actually it is rather easy using ``forward substitution for L and backward substitution U``. This is makes it relatively easy to solve for ``y`` and afterwards ``x``. 

## Forward substitution by an example

$$
\begin{bmatrix}
1 & 0 & 0 & 0 \\
L_{1,0} & 1 & 0 & 0 \\
L_{2,0} & L_{2,1} & 1 & 0 \\
L_{3,0} & L_{3,1} & L_{3,2} & 1
\end{bmatrix} \cdot \begin{bmatrix}
y_0 \\
y_1 \\
y_2 \\
\end{bmatrix} = \begin{bmatrix}
b_0 \\
b_1 \\
b_2 \\
\end{bmatrix}
$$

$$ b_0 = y_1$$
$$ b_1 = L_{1,0} \cdot y_0 + y_1$$
$$ b_2 = L_{2,0} \cdot y_0 + L_{1,0} \cdot y_1 + y_2$$

Now to find $\mathbf{y}$

$$ y_i = b_i - \sum_{k=0}^i L_{i,k} \cdot y_k$$

## Backward substitution by an example

$$
\begin{bmatrix}
U_{0,0} & U_{0,1} & U_{0,2}  \\
0 & U_{1,1} & U_{1,2}  \\
0 & 0 & U_{2,2}  \\
\end{bmatrix}  \cdot \begin{bmatrix}
x_0 \\
x_1 \\
x_2 \\
\end{bmatrix} = \begin{bmatrix}
y_0 \\
y_1 \\
y_2 \\
\end{bmatrix}
$$

$$ y_2 = U_{2,2} x_2$$
$$ y_1 = U_{1,1} x_1 + U_{1,2} \cdot x_2 $$
$$ y_0 = U_{0,0} x_0 + U_{0,1} x_1 + U_{0,2} \cdot x_2 $$

Now to find $\mathbf{x}$, starting from the end say $N-1$

$$ x_i = \frac{y_i - \sum_{k=i+1}^{N-1} U_{i,k} \cdot y_k}{U_{i,i}}$$

## Why is it usefull?
We can store the matrix A in L U after use, this has computation complexity $O(N^3)$ and for every time we want to solve a linear system we got $2*O(N^2)$.

We if want to determine the determinant of A we simple multiply the product of ``U``

If we want to add more system we rember that is it recursive as stated earlier.

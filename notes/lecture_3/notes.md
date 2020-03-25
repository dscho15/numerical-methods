# Tridiagonal and Band-diagonal systems of Equations

This sort of matrix has values on the the diagonal and 1 above it and below it. Back and forward substituion only takes O(N) in this case.

$$\begin{bmatrix}
b_0 & c_0 & 0 \\
a_1 & b_1 & \cdots \\
\cdots & \cdots
\end{bmatrix} * \begin{bmatrix} u_0 \\ u_1 \\ \cdots \end{bmatrix} = \begin{bmatrix} r_0 \\ r_1 \end{bmatrix} $$ 

The idea is then to split it up in even and odd elements e.g.

$$
\begin{bmatrix}
b_0 & c_0 & 0 & 0 \\
a_1 & b_1 & c_1 & 0 \\
0 & a_1 & b_2 & c_2
\end{bmatrix} \begin{bmatrix} u_0 \\ u_2 \\ u_4 \\ u_6 \\ u_1 \\ u_3 \\ u_5 \end{bmatrix}$$

, this accomplished by permutating rows and columns. The idea is then to subtract all the elements under the diagonal and equal them to 0, thus ending up with only have elemenets above the diagonal. Which easily can be solved in parallel

$$\begin{bmatrix}
b_0 & 0 & 0 & 0  & c_0  \\
0 & b_1 & 0 & 0  & a_2 & c_2 \\
0 & 0 & b_2 & 0  & 0 &
\end{bmatrix} \begin{bmatrix} u_0 \\ u_2 \\ u_4 \\ u_6 \\ u_1 \\ u_3 \\ u_5 \end{bmatrix}$$


# Band-Diagonal Systems

This matrix is a special form of the tridiagonal matrix and is basically a way of solving so, so a special setup.

# Iterative Improvement of a Solution to Linear Equations

For large sets of linear equations it is not always easy to obtain precision equal to or even comparable to the computer's limit hence roundoff errors accumulate. This method assumes that the matrix was not close to being singular e.g.

$$ A \cdot x = b$$

all was know is some slighly wrong solution $x + \delta x$ 

$$ A \cdot (x+\delta x) = (b+\delta b)$$

So the basic idea is subtracting b and then seeing how big the error is, thereby figuring out what x is supposed to be. However what if A is not exact? Let's say we ran LU-decomposition on A, then we have to define some residual of A. So say we can find some matrix B, which is close to the inverse of A e.g.

$$ R = I - B \cdot A $$
$$ B A = I-R $$
$$ A^{-1} = (I-R)^{-1}B $$
$$ B_n = (1 + R + R^2...) B_0 $$

where R can be expanded as a partial sum, now there can be made a new recursive method using this B.

$$x_{n+1} = x_n + B_0 * (R - A*x_n)$$

The only requirement is that the residual is small.

# Linear Independency

Vectors are linear independent if and only if 

$$c_1 * x_1 + c_2 * x_2 ... c_k * x_k = 0$$

where $x$ is a vector, and they are $R^n$, which implies that 

$$c_1 = c_2 = 0..$$

**A be  an  arbitrarym×nmatrix  where m≥n.   Now perceive  the  columns  of A as  vectors  in Rm.   The  maximum  number  of linearly independent columns in A is denoted the rank of A.  Furthermore, if then columns of A all are linearly independent A is said to have full rank, i.e.  the rank is n.**

# Testing Singularity

Let $A$ be an $m\times n$ matrix over some field $F$. Recall that $Ax=0$ always has the tuple of 0's as a solution. This solution is called the trivial solution. All other solutions are called nontrivial. We can apply **row -echolean** form to determine the non-trivial solutions.

What is a singular vector?

## Why is a singular matrix not invertible?

say, some square matrix $ A \epsilon R^{nxn}$

$$ Ax=0$$
$$ A^{-1} * A x = A^{-1} 0$$
$$ I x = 0$$

surely this can not be the case if x is a none zero vector.

# Orthonormal

Vectors are orthogonal if the dot product is equal to 0, furthermore if the vectors are unit length then they considered as a *orthonormal* vector pair.

They have a special property if we consider som orthonormal matrix, then 

    It's inverse is its transpose
    It is always full rank.

# Column Space

The column space contains all linear combinations of the columns of A.  It is a subspace of $R^m$. 

We can describe all combinations of the two columns geometrically: $Ax=b$ can be solved if and only if b lies in the plane that is spanned by the two column vectors. This is the thin set of attain able $b$. If $b$ lies off the plane, then it is not a combination of the two columns. In that case $Ax=b$ has no solution

The plane is a subset of $R^3$ and the subspace is denoted C(A).

# All spaces we are concerned with

The column space of A is denoted by C(A). Its dimension is the rank r.
  
The nullspace of A is denoted by N(A). Its dimension is n-r 

The row space of A is the column space of $A^T$. It is $C(A^T)$, and it is spanned by the rows of A. Its dimension is also r.

The left null space of A is the nullspace of AT.  It contains all vectors y such that $A^Ty=0$, and it is written $N(A^T)$. Its dimension is n-r.

The row-space and nullspace are subset of $R^n$, takes some vector as input.
The column-space is the transposed $(m\times n)$ and is therefore $R^m$ vector as input.

# Online stuff

In this lecture we learn what it means for vectors, bases and subspaces to be orthogonal. The symbol for this is ⊥. The “big picture” of this course is that the row space of a matrix’ is orthog­onal to its nullspace, and its column space is orthogonal to its left nullspace. 

https://ocw.mit.edu/courses/mathematics/18-06sc-linear-algebra-fall-2011/least-squares-determinants-and-eigenvalues/orthogonal-vectors-and-subspaces/MIT18_06SCF11_Ses2.1sum.pdf

---

$$ N(A^T A) = N(A)$$

Example

$$ A = \begin{bmatrix} 1 & 2 & 5 \\ 2 & 4 & 10 \end{bmatrix}$$


The row space has dimension 1 and the same goes for the column space. The nullspace has dimension n-r, where r is the dimension of row space. Therefore the null space has dimension 2.

$$ \begin{bmatrix} 1 & 2 & 5 \\ 2 & 4 & 10 \end{bmatrix} u = 0$$
$$ u_1 + 2u_2 + 5u_3= 0 $$

This has two free variables, say u1 and u2 thereby we define the nullspace as


$$ \bar{u} = N(A) =  u1 * \begin{bmatrix}  1 \\ 
                                    0 \\ -1/5 \end{bmatrix} + u2 * \begin{bmatrix}  0 \\ 
                                    1 \\ -2/5 \end{bmatrix} $$

Now let's try to draw this


# Independence

Two vectors are independent if they not the same vector

# Basis vector

#





# SVD topic

https://blog.statsbot.co/singular-value-decomposition-tutorial-52c695315254

# Orthonormal Matrix

An orthonormal matrix can be considered as a rotation matrix, it has the properties of unit length on each column or row and the transpose is it's invert.

$$RR^{T} = I$$

# Back to SVD

We consider an arbitrary matrix A, which we write as $U W V^T$, where U is a orthonormal column matrix ($m\times n$), V is a column orthonormal matrix ($n\times n$) and when transposed is converted to a row matrix, W is a diagonal matrix, *non-negative!!*.

V is a special matrix, hence it typically composes information about the null-space.

## Explaination of Thereom 7

We consider matrix A and we assume that $W$ has k elements, given the rank of $A$ is $k$. The rest of the elements (singular values / eigenvalues) are 0.

$N(A)$, the null space of A must then have dimension $n-k$, where the $n-k$ V columns forms a basis for the null-space. 

$B(A)$ is the set of $y$'s that forms a basis for $B(A)$, which is determined by the first $k$ orthogonal basis vectors.

The SVD solution $x=V \tilde{W}^{-1}U^{T}b$ is the SVD solution, which is merely taking the tranpose of V, invert W and transpose U. The diagonal elements of $W$ could however be 0 and the invert of a diagonal matrix is all the elements inverted.

# Example of inverting a diagonal matrix

Is this possible? Well, no, it does not have a invert matrix.

$$ A = \begin{matrix}
2 & 0 \\
0 & 0 
\end{matrix} $$



# Normal Equation - what it is

First of all, we face some linear system of equation problem

$$ Ax = b$$

The normal equation is that, which minimizes the sum of squared differences on each side

$$ A^T Ax = A^T b $$

# Normal Equation - why it is bad

#  SVD topic





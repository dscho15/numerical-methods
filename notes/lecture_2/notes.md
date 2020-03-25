# General Linear Least Square

The idea is to fit a set of data points, say $(x_i,y_i)$ to a model that is not just a linear combination of $x$, but rather a linear combination of M specified functions of $x$

Example:

$$y(x) = a_o+a_1x+a_2x^2...$$

which, is some $M-1$ order polynomium

It can be written in the form:

$$y(x) = \sum^{M-1}_{k=0} a_k X_k(x)$$

where $X_k(x)$ is some $k$'th *basis function* and $a$ is the unknown parameter. 

Then we define the merit function, **a function that  measures the difference between the parametrized function y and the measurement**

$$\chi^2 = \sum^{N-1}_{i=0} [ \frac{y_i-\sum^{M-1}_{k=0}a_kX_k(x_i)}{\sigma_i}]^2$$

, where $\sigma_i$ is the measurement error, if unknown set them to 1. The idea is to minimize $\chi^2$, finding the minima. 

# First idea

Define a matrix $A$, whose components $N \times M$ are constructed from the M basis functions, evalulated at $N$ abscissas(X-axes) $x_i$, and from the N measurement errors $\sigma_i$, this is formulated by the following prescription

$$ A_{ij} = \frac{X_j(x_i)}{\sigma_i}$$

This matrix is called the **Design Matrix: A**, e.g. the columns are the basis functions $X_k$ and the **M** rows are the measurements. We also define a vector **b** of length **N** 

$$ b_i = y_i/\sigma_i$$

# Solution / Finding minima

The minima of the merit function is found when the derivate with respect to all M parameter $a_k$ vanishes, which yields the following equation

$$ 0 = \sum^{N-1}_{i=0} \frac{1}{\sigma_i^2} [ y_i - \sum^{M-1}_{j=0} a_j X_j(x_i) ]X_k$$

We would like to transform it into matrix-form and then solve it, firstly I will define A

$$ A = \begin{bmatrix}
\frac{X_0(x_0)}{\sigma_0} & \frac{X_1(x_0)}{\sigma_0} & \cdots\\
\frac{X_0(x_1)}{\sigma_1} & \frac{X_1(x_1)}{\sigma_1} & \cdots
\end{bmatrix}$$

$$ \beta_k = \sum_{j=0}^{M-1} \alpha_{kj} a_j$$

where

$$ \beta_k = \sum_{i=0}^{N-1} \frac{y_iX_k(x_i)}{\sigma_i}$$

which is equivalent to $\beta = A^T * b$ and then

$$ \alpha_{kj} = \sum_{i=0}^{N-1} \frac{X_j(x_i)X_k(x_i)}{\sigma_i}$$

which is equivalent to $\alpha = A^T * A$, now we have a system of linear equations abbreviated to

$$ (A^T*A)*a=A^T * b$$

the inverse matrix of $\alpha$ denoted $C$ is the covariance matrix and describes the untercainties of the estimated parameters of a. This is derived by firstly considering

$$ a_j = \sum^{M-1}_{k=0} \alpha_{jk}^{-1}\beta_k$$

, and the variance associated with the estimate of $a_j$ is

$$\sigma^2(a_j) = \sum^{N-1}_{i=0} \sigma_i^2 (\frac{\partial a_j}{\partial y_i})^2$$

and by doing some calculus it is possible to derive that $\sigma^2(a_j) = C_{jj}$. Since $A_t * A$ is positve-definite, Cholesky decomposition allows most effiecient way to solve the normail equations. 

# Cholesky Decomposition

This method requires the matrix to be positive-definite and symmetric usually 

$$ A^T*A $$
is consindered. Symmetric means $a_{ij} = a_{ji}$ so basically the transpose and that $v * A * v > 0$ for all vectors v besides the zero vector. Besides this all the eigenvalues are positive.

Cholesk is about 2 times faster than LU and the special thing is that only a lower triangle matrix is used and then it is transposed.

$$ A = L*L^T$$

It is derived quite simple and the square root is used on the diagonal.

$$ L_{ii} = (a_{ii} - \sum^{i-1}_{k=0} L_{ik})^{0.5}$$

and

$$ L_{ji} = \frac{1}{L_{ii}} (a_{ij} - \sum^{i-1}_{k=0}L_{ik}L_{jk}).
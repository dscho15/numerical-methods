# What is root finding?

What value of $x$ does $f(x)=0$?

If we let $f(x)=ax^2+bx+c=0$

that can happend multiple times, so we have multiple solutions, some however have analytical solution and some have not such as 

$f(x) = e^x-x=0$

it could also be that it isn't a function.

# Bracketing Methods

Finding a single root that falls within a known range. They are robust, and must know something ahead of time.

$x_l <= x_r <= x_u$

IT IS REQUIRED THAT THERE IS ONLY ONE ROOT.

# Open Methods
They don't have upper nor lower bounds. They require an initial guess. And they are considred more effiecient than bracketing methods. They can be unstable and not find a solution (such as local minima). 

# Root of Polynomials
Algorithm is specific to polynomials. The problem is fit to a polynomail and find the roots of this.

# Generalized Root Finding Algorithms
We wish to find all values of x such that $f(x) = a $, that is simple, we just compute the offset and make a new function.

# Use case, need optimized algorithms

In many circumstances a single computation of f(x) may take hours, days or weeks! In these cases it is higly desired to minimize the total number of compuations of f(x)!

# Convergence Properties

If a method converges to a solution $\hat{x}$ we write the error after k steps as
$$\epsilon_k = x_k - \hat{\epsilon_k}$$
if then $\frac{|\epsilon_{k+1}|}{|\epsilon_k|^{\alpha}} = C$ then we say that the convergence has order $\alpha$ with convergence constant C.


# Open Methods or iterative - Newton Raphson

Given a initial point, $x_0$ it will try to find a root:

$$x_1 = x_0 - \frac{f(x_0)}{f'(x_0)}$$

We may picture the algorithm as we the tangent line by setting y equal to 0, $y = f'(x)(x-x_n)+f(x_n)$. We stop when we hit a desired accuracy.

So the methods requires to calculate the derivative. This is a weak property. Besides this what occurs when a local minima occurs? And it might never find a solution if it is between a local minima and maxima.


# Bisection (Bracketing Method)
Given the root limits or the so called brackets, we know where the root is located. Now we keep dividing until we find the root.

$c = \frac{a+b}{2}$

if $f(c)f(b) < 0 $
$ a = c$ 
else
$ b = c$
$|(dx)| < \epsilon_{acc}$ 

- Very Robuist Root Finder
- Least efiecient Root Finder
- Guaranted to find a root as long as the bounds span a crossing
- Sometimes good to check sign change of bounds

# Secant Method (Open or Bracket.)
Given the root limits or the so called brackets, we would like to know where the root is located.

$$x_{i+1} = x_i - \frac{f(x_i)(x_i-x_{i-1})}{f(x_i)-f(x_{i-1})}$$

We start with an initial guess, we use the formula, changes either the upper or lower bound, and we look if the new value is close to our accuracy. 

Then we find

$$|\epsilon_a| = |\frac{x_{i+1}-x_i}{x_{i+1}}|*100$$

Now we stop if it is low enough.

Where we estimate the relative approximate error. The idea behind it is so draw a triangle 

$$ \frac{AB}{AE} = \frac{DC}{DE}$$ shown in the pdf.

http://mathforcollege.com/nm/mws/gen/03nle/mws_gen_nle_txt_secant.pdf

The convergence rate is 

$$ \alpha = 1/2(1+\sqrt(5))$$

# False Position
It is a bracketing method, inspired by bisection

https://www.youtube.com/watch?v=pg1I8AG59Ik

Given some interval [a,b] then we compute the secant between the points. When you compute the point c, where it crosses the x-axis. Depending on the sign, a or b is replaced and this repeats until the root is found.

The algorithm, or update is

$$x_{i+1} = x_i - \frac{x_i-y_i}{f(x_i)-f(y_i)}f(x_i)$$

The convergence properties, if the function is twice differentiable then $\alpha = 1$ and the convergence constant is $C=1-f'(\hat{x})\frac{y_m-\hat{x}}{f(y_m)}$

# Ridders Method

http://physics.wm.edu/~evmik/classes/matlab_book/ch_root_finding/ch_root_finding.pdf
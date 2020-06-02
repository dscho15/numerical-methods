# Important things

We are to give a numerical solution to ordinary differential equations, where we define the stepsize h 

$$x_n = x_0 + n\cdot h$$
$$y_n \approx y(x_n) $$
where $y_n$ is a numerical approximation to the true value $y(x_n)$

The aim is always to evaluate as few f computations as possible to integrate.

# Runge Kutta Methods (one step)

## 1st order
$$y_{n+1} = y_n + h\cdot f(x_n,y_n)$$
, where the order is defined as $O(h^2)$, however is a first order method to the general solution 

## 2nd order
$$k_1 = f(x_n,y_n)$$
$$k_2 = f(x_n+1/2h, y_n+h/2 \cdot k_1)$$
$$y_{n+1} = y_n + h \cdot k_2$$
, where the order is defined as $O(h^3)$, however is a second order method to the general solution 

## 4th order
$$k_1 = f(x_n,y_n)$$
$$k_2 = f(x_n+h/2 , y_n + h/2 \cdot k_1)$$
$$k_3 = f(x_n+h/2,y_n + h/2 \cdot k_2)$$
$$k_4 = f(x_n+h, y_n+h \cdot k_4)$$
$$y_{n+1} = y_n + h/6(k_1+2k_2+2k_3+k_4)$$
, where the order is defined as $O(h^5)$, however is a fourth order method to the general solution 

## Trapezoidal Method

This is a implicit one step method, aka  $y_{n+1} = F(x,h,y_n,y_{n+1},f)$

$$y_{n+1} = y_n + \frac{h}{2} \left( f(x_n,y_n) + f(x_{n+1},y_{n+1}) \right) + O(h^3)$$
, but how do we know the next step? We perform an euler step

$$y=y^{*}_{n+1} = y_n + hf(x_n,y_n)$$

and we then use Newtons method to solve the system of non-linear equations

$$ O(y) = y-y_n-\frac{h}{2}(f(x_n,y_n) + f(x_{n+1},y)) = 0 $$
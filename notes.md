
## Polynomial multiplication using FFT

Given polynomials f(x) and g(x) of degree n, one can output product polynomial fg(x) in n^2 time. But turns out using FFT one can do this in N log N time. 

A simple observation plus the fact that one can evaluate any polynomial at different powers of $2n^{th}$ root of unity in N log N time may hint something. 

The observation is "interpolation of n data points into a polynomial of n-1 degree give a unique polynomial".  With this observation and the fact that resulting polynomial will be degree 2n - 1 we can first FFT transform f and g with 2nth root of unity (i.e. evaluate f(k) and g(k) for $k \in [0, \omega_{2n}, ..., \omega_{2n}^{2n-1}]$). as FFT(f) and FFT(g). Then we take the Hadamard product $FFT(f) \cdot FFT(g)$ to produce evaluation 2n-1 data points for polynomial $fg(x)$. Observe that the resulting 2n-1 points are $fg(x)$ evaluated at $k \in [0, \omega_{2n}, ..., \omega_{2n}^{2n-1}]$. 

To get the coefficient of fg(x) we simply perform FFT inverse as $FFT^{-1}(FFT(f) \cdot FFT(g))$. To understand why FFT inverse results in coefficients (i.e. is equivalent to interpolation) one may view FFT as matrix multiplication:

$$
\begin{matrix}  
1 & 1 & 1 & ... & 1 \\  
\omega_{2n} & \omega_{2n} & \omega_{2n} & ... & \omega_{2n} \\  
\omega_{2n}^2 & \omega_{2n}^2 & \omega_{2n}^2 & ... & \omega_{2n}^2 \\  
\vdots \\
\omega_{2n}^{2n-1} & \omega_{2n}^{2n-1} & \omega_{2n}^{2n-1} & ... & \omega_{2n}^{2n-1}
\end{matrix}
\cdot
\begin{matrix}
a_0 \\
a_1 \\
a_2 \\
\vdots \\
a_{2n-1}
\end{matrix}
=
\begin{matrix}
\hat{a_0} \\
\hat{a_1} \\
\hat{a_2} \\
\vdots \\
\hat{a_{2n-1}} \\
\end{matrix}
$$
where \omega_{2n} is 2n^{th} root of unity, $a_0, a_1,...a_{2n-1}$ are coefficients of polynomial fg(x) and resulting column vector is contains fg(\omega_{2n}^k) at row k.

The left most matrix has a special name called Vandemont matrix. Let us denote it with M. Given that $FFT(f) FFT(g)$ gives us column vector on RHS, we can calculate the column vector consisting of coefficients by calculating $M^{-1} \cdot (FFT(f) FFT(g))$. 

Due to "Inversion theorem":
$$M(\omega_{2n})^{-1} = \frac{1}{n}M(\omega_{2n}^{-1})$$
Thus to calculate the coefficients we can calculate: 
$$
\begin{matrix}  
1 & 1 & 1 & ... & 1 \\  
\omega_{2n}^{-1} & \omega_{2n}^{-1} & \omega_{2n}^{-1} & ... & \omega_{2n}^{-1} \\
\omega_{2n}^{-2} & \omega_{2n}^{-2} & \omega_{2n}^{-2} & ... & \omega_{2n}^{-2} \\
\vdots \\
\omega_{2n}^{-(2n-1)} & \omega_{2n}^{-(2n-1)} & \omega_{2n}^{-(2n-1)} & ... & \omega_{2n}^{(-2n-1)}
\end{matrix}
\cdot
\begin{matrix}
\hat{a_0} \\
\hat{a_1} \\
\hat{a_2} \\
\vdots \\
\hat{a_{2n-1}} \\
\end{matrix}
=
\begin{matrix}
a_0 \\
a_1 \\
a_2 \\
\vdots \\
a_{2n-1}
\end{matrix}
$$

Notice that this corresponds to FFT operation on vector FFT(f)FFT(g) with $\omega_{2n}^{-1}$ instead of $\omega_{2n}$. Thus one can find the coefficients of resulting polynomial by evaluating $FFT^{-1}(FFT(f)FFT(g))$ where $FFT^{-1}$ is FFT algorithm with $\omega_{2n}^{-1}$.

Reference: 
https://www.cs.toronto.edu/~denisp/csc373/docs/tutorial3-adv-writeup.pdf

## Cyclic multiplication using FFT

Let f and g be vectors, then we denote the result $FFT^{-1}(FFT(f)FFT(g))$ as convolution of vectors f and g. 

The resulting vector contains: 
$$\sum_{k=0}^{N} a_kb_{j-k\mod{N}}$$
at j^th row, where $a_k$ are coefficients of f and b_{k} are coefficients of g. 

This is equivalent to cyclic multiplication modulo $X^N - 1$.


## Negacylic multiplication using FFT

For modulo $X^N + 1$ (i.e. negacylic multiplication) one cannot directly use FFT because the terms wrap around and negate. 
### Double trick

Notice that $(X^{2N} - 1)$ can be factorised into (X^N-1) and X^N+1. This implies, due to chinese remainder theorem, that each value in modulo X^{2N} -1 has unique representation in modulo X^N+1. 

If one can re-write (pre-image map) a value \in Z[X]/X^N+1 as a value \in Z[X]/X^{2N}-1, then they can perform cyclic multiplication and reduce the result to Z[X]/X^N+1

Let's define the pre-image mapping as: f(1-x^N). 

Reducing f(1-x^N) in Z[X]/X^{2N}-1 to Z[X]/X^{N}+1 will negate the upper coefficients of f (X^N = -1) and add same lower coefficient of f, thus resulting in 2f. 

To perform nega-cyclic multiplication we pre-image map f and g as $f* = f - fx^N$ and $g* = g - gx^N$ $\in Z[X]/X^{2N-1}$. Pre-image mapping can be done easily by negating coefficient vector of polynomial and append it to itself. 

Then we perform cyclic multiplication $f*g*$ using FFT. Then we reduce $f*g*$ modulo $X^N+1$. Notice that the resulting value is
$$f* \cdot g* \mod X^N+1 = 2f \cdot 2g = 4fg$$
Thus we divide the reduce polynomial in $Z[X]/X^N+1$ by 4. 

Reduction at the end is trivial if one observers the following equation: 
$$f*g*= -f(x^N-1) \cdot -g(x^N-1)$$
$$f*g*= fg(x^N-1)^2 = fg(x^{2N} + 1 - 2x^N) = 2(fg - fgx^N)$$
This implies the resulting product in bigger ring is pre-image mapping of product in smaller ring. Thus  to reduce we can do the opposite of pre-image. That is truncate upper half of coefficient vector. Then we divide by 2 to get $fg$.

### Twisting trick




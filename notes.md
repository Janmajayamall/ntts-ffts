
## Polynomial multiplication using FFT

Given polynomials $f(X)$ and $g(X)$ of degree $N$, one can output product polynomial $fg(X)$ in $N^2$ time. But turns out using FFT one can do this in $N log N$ time. 

A simple observation plus the fact that one can evaluate any polynomial at different powers of $2N^{th}$ root of unity in $N log N$ time may hint something. 

The observation is "interpolation of $N$ data points into a polynomial of degree $N-1$ is unique".  With this observation and the fact that resulting polynomial will be degree $2N - 1$ we can first apply FFT transform to $f$ and $g$ with $2n^{th}$ root of unity (i.e. evaluate $f(k)$ and $g(k)$ for $k \in [0, \omega_{2n}, ..., \omega_{2n}^{2n-1}]$). Then we take the Hadamard product $FFT(f) \odot FFT(g)$ to produce $2N-1$ evaluation points for polynomial $fg(X)$. Observe that the resulting $2N-1$ points equal $fg(x)$ evaluated at $k \in [0, \omega_{2n}, ..., \omega_{2n}^{2n-1}]$. 

To find coefficients of $fg(X)$ one can inverse FFT transform as $FFT^{-1}(FFT(f) \odot FFT(g))$. To understand why inverse transformation results in coefficients (i.e. is equivalent to interpolation) one may view FFT transform as matrix multiplication:

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
where $\omega_{2N}$ is $2N^{th}$ root of unity, $a_0, a_1,...a_{2n-1}$ are coefficients of polynomial $fg(X)$ and resulting column vector contains $fg(\omega_{2N}^k)$ at row $k$.

The left most matrix has a special name called Vandermonde matrix. Let us denote it with $M$. Given that $FFT(f) \odot FFT(g)$ results in column vector on RHS, we can calculate the column vector consisting of coefficients by  $M^{-1} \cdot (FFT(f) \odot FFT(g))$. 

Due to "Inversion theorem":
$$M(\omega_{2N})^{-1} = \frac{1}{n}M(\omega_{2N}^{-1})$$
Thus to calculate the coefficients one may calculate: 
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

The structure above corresponds to FFT transformation on vector $FFT(f) \odot FFT(g)$ with $\omega_{2n}^{-1}$ instead of $\omega_{2n}$. Thus one can find the coefficients of resulting polynomial by evaluating $FFT^{-1}(FFT(f) \odot FFT(g))$ where $FFT^{-1}$ is FFT algorithm with $\omega_{2n}^{-1}$.

Reference: 
https://www.cs.toronto.edu/~denisp/csc373/docs/tutorial3-adv-writeup.pdf

## Cyclic multiplication using FFT

Let $f$ and $g$ be vectors. We denote the result $FFT^{-1}(FFT(f) \odot FFT(g))$ as convolution of vectors $f$ and $g$. 

The resulting vector equals: 
$$\sum_{k=0}^{N} a_kb_{j-k\mod{N}}$$
at $j^{th}$ row, where $a_k$ are coefficients of $f$ and $b_{k}$ are coefficients of $g$. 

This is equivalent to cyclic multiplication $\mod {X^N - 1}$.

## Negacyclic multiplication using FFT

For modulo $X^N + 1$ (i.e. negacyclic multiplication) one cannot directly use FFT because the terms wrap around and negate. 
### Double trick

Notice that $(X^{2N} - 1)$ can be factorised into $(X^N-1)$ and $(X^N+1)$. This implies, due to Chinese remainder theorem, that each value in modulo $X^{2N} -1$ has unique representation in modulo $X^N+1$. 

If one can re-write (pre-image map) a value $\in Z[X]/X^N+1$ as a unique value $\in Z[X]/X^{2N}-1$, then they can perform cyclic multiplication and reduce the result to $Z[X]/X^N+1$.

Let's define the pre-image mapping as: $f(X)(1-X^N)$, where $f(X) \in Z[X]/X^N+1$

Reducing $f(X)-X^Nf(x)) \in Z[X]/X^{2N}-1$ to $Z[X]/X^{N}+1$ will negate the upper half coefficients of polynomial (since $X^N = -1$ and lower half coefficients a variable exponents smaller than $N$), thus resulting in $f(X) + f(X) = 2f(X)$. 

To calculate $f(X)g(X) \mod{X^N+1}$, we pre-image map $f(X)$ and $g(X)$ as $f'(X) = f(X) - f(X)X^N, \space g'(X) = g(X) - g(X)X^N$ $\in Z[X]/X^{2N-1}$. Pre-image mapping can be done easily by negating coefficient vector of polynomial and append it to itself. 

We calculate $f'(X)g'(X)$ using FFT and reduce the result $\mod{X^N+1}$. Now notice that the resulting value equals
$$f'(X)g'(X) \mod X^N+1 = 2f(X) \cdot 2g(X) = 4fg$$
Thus we divide the reduced polynomial in $Z[X]/X^N+1$ by 4. 

Reduction at the end is trivial if one observers the following: 
$$f'(X)g'(X)= -f(X)(X^N-1) \cdot -g(X)(X^N-1)$$
$$f(X)g(X)(X^N-1)^2 = f(X)g(X)(X^{2N} + 1 - 2X^N) = 2(f(X)g(X) - f(X)g(X)X^N)$$
This implies the resulting product $\in Z[X]/X^{2N}-1$ is pre-image mapping of product $\in Z[X]/X^N+1$. Thus  to reduce. one can invert the per-image image. That is truncate upper half of coefficient vector. Then we divide by 2 to get $f(X)g(X)$.
### Twisting trick

$X^N+1$ can be factorised into co-prime factors $(X^{N/2}+i) (X^{N/2}-i)$. This implies every value in $X^{N}+1$ has unique representation $X^{N/2}-i$. Hence, one can calculate negacyclic polynomial multiplication in ring $\mod {X^{N/2}-i}$ and map it back to ring $\mod{X^N+1}$.

But FFTs work directly only on rings with structure $Z[X]/X^k-1$ (i.e. cyclic polynomials). Thus we need a map from ring $\mod {X^{N/2}} + i$ to ring $\mod X^{N/2}-1$. Turns a trivial map exists if we map $X$ of $g(X) \in Z[X]/X^{N/2}-i$ to $\omega_{4(N/2)}X$. 

To see why, let: 
$$
f(X) = g(X) + h(X)(X^{N/2}-i)
$$
where $h(X)$ is quotient and $f(X)$ is original polynomial. 

Mapping $X$ to $\omega_{4(N/2)}X$:
$$f(\omega_{4(N/2)}X) = g(\omega_{4(N/2)}X) + h(\omega_{4(N/2)}X) ((\omega_{4(N/2)}X)^{N/2}-i)$$
Since $\omega_{4(N/2)} = e^{\frac{2i\pi}{4(N/2)}}$,
$$f(\omega_{4(N/2)}X) = g(\omega_{4(N/2)}X) + h(\omega_{4(N/2)}X) ((e^{\frac{\pi i}{2}}X^{N/2}-i)$$
Since $e^{\frac{\pi i}{2}} = i$,
$$f(\omega_{4(N/2)}X) = g(\omega_{4(N/2)}X) + ih(\omega_{4(N/2)}X) (X^{N/2}-1)$$
$$f(\omega_{4(N/2)}X) = g(\omega_{4(N/2)}X) \mod{(X^{N/2}-1)}$$

Mapping $X^N+1$ to $X^{N/2}-i$ can be done by interpreting $X^{N/2} = i$. This implies, to reduce an element $\in Z[X]/X^N+1$ to an element in $Z[X]/X^{N/2} - i$, multiply all coefficients in upper half by $i$ and then halve their exponents and add them to their corresponding coefficients in lower half. Thus, the upper half coefficient vector becomes the imaginary parts of lower half coefficient vector. 

For second mapping, observe that 
$$g(\omega_{4(N/2)}X) = a_0 + a_1(\omega_{4(N/2)}X) + a_2(\omega_{4(N/2)}X)^2 + ... + a_{(N/2)-1}(\omega_{4(N/2)}X)^{(N/2)-1}$$
$$= a_0 + \omega_{4(N/2)}a_1(X) + \omega_{4(N/2)}^2a_2(X)^2 + ... + \omega_{4(N/2)}^{(N/2)-1}a_{(N/2)-1}(X)^{(N/2)-1}$$
One can simply map $g(X)$ to $g(\omega_{4(N/2)}X)$ by multiplying all coefficients element-wise with vector: 
$$[1, \omega_{4(N/2)},\omega_{4(N/2)}^2, ..., \omega_{4(N/2)}^{(N/2)-1}]$$

After mapping $f,g \in Z[X]/X^N+1$ to element in $Z[X]/X^{N/2}-1$, one can perform FFT to get the product in $Z[X]/X^{N/2}-1$. Then they can invert the maps to get product in $Z[X]/X^{N}+1$. 

Inverse map for the second mapping is as simple as multiplying the coefficients of element with inverse of $\omega_{4(N/2)}$ powers vector. And inverse mapping from $Z[X]/X^{N/2}-i$ to $Z[X]/X^N+1$ is reverse of mapping other way around. That is, take the imaginary part of coefficients and append them to the coefficient array as real parts for upper half. 




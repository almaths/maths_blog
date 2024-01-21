---
layout: post
title:  "A Riemann-Hilbert derivation of the Christoffel-Darboux formula"
date:   2023-12-26 18:50:11 +0100
permalink: /christoffel-darboux/
katex: True
---

In my [recent work](https://arxiv.org/abs/2306.14107) on skew-orthogonal polynomials I was interested in deriving a Christoffel-Darboux type formula. To this end Thomas Bothner referred me to a [paper](https://arxiv.org/abs/1407.2597) of his in collaboration with Marco Bertola which contains a nice derivation of the Christoffel-Darboux formula using only complex analysis arguments (Theorem 2.8). I will present this argument in the simpler context of orthogonal polynomials rather than that of biorthogonal polynomials considered by the authors.

Let us recall the Riemann-Hilbert problem for orthogonal polynomials. Given a measurable function $$w : \mathbb{R} \to [0,+\infty]$$ such that $$\int_\mathbb{R} w(x) \lvert x \rvert ^k \, dx < +\infty$$ for all $$ k \in \mathbb{N}$$ and where $$w(x) > 0$$ on some open set we construct a sequence of monic polynomials $$P_n(x) = x^n + \mathcal{O}(x^{n-1})$$ such that

$$ \int_\mathbb{R} P_n(x)x^k w(x)\, dx = 0 $$

for all $$k=0, \dots, n-1$$. We call $$P_n$$ the $$n$$th monic orthogonal polynomial with respect to $$w$$.

We let $$h_n = \int_\mathbb{R} P_n(x)x^n w(x)\, dx$$ be the (squared) $$L^2(w)$$ norm of $$P_n$$. We may reformulate this as a Riemann-Hilbert problem.

**Riemann-Hilbert problem for orthogonal polynomials ([Fokas-Its-Kitaev, 1992](https://link.springer.com/article/10.1007/BF02096594)):**

Find a matrix valued function $$X_n : \mathbb{C} \setminus \mathbb{R} \to \mathbb{C}^{2 \times 2}$$ such that

1. $$X_n$$ is analytic (entry-wise) on $$\mathbb{C} \setminus \mathbb{R}$$.
2. $$X_n$$ has continuous non-tangential boundary values up to $$\mathbb{R}$$ from above ($$+$$) and below ($$-$$). We label these $$X_n^\pm (x) = \lim_{\epsilon \downarrow 0} X_n(x\pm i \epsilon)$$ for $$ x \in \mathbb{R}$$.
3. These boundary values are related by the jump condition
$$X_n^+(x) = X_n^-(x) \left( \begin{matrix} 1 & w(x) \\ 0 & 1\end{matrix} \right)$$
4. Finally, $$X_n$$ is normalised at infinity by the scaling as $$z \to \infty$$

$$X_n(z) = \left( \mathbb{I}+\mathcal{O}(z^{-1}) \right) \left( \begin{matrix} z^n & 0 \\ 0 & z^{-n} \end{matrix} \right)$$

**Proposition:** The above RHP has a unique solution given by (for $$n\geq 1$$)

$$ X_n(z) = \left( \begin{matrix} P_n(z) & C_\mathbb{R} \left( P_n w\right)(z) \\ - 2\pi i h_{n-1}^{-1} P_{n-1}(z) & - 2\pi i h_{n-1}^{-1} C_\mathbb{R} \left( P_{n-1} w\right)(z) \end{matrix} \right)$$

where

$$C_\mathbb{R}(f)(z) = \frac{1}{2\pi i} \int_\mathbb{R} \frac{f(x)}{x-z}\, dx$$

is the Cauchy transform of the function $$f$$. Furthermore $$\det X_n(z) = 1$$ identically. If $$n=0$$ the solution is $$ X_0(z) = \left( \begin{matrix} 1 & C_\mathbb{R} \left(  w\right)(z) \\ 0 & 1 \end{matrix} \right)$$. $$\triangle$$

That $$\det X_n(z) = 1$$ can be seen directly from the RHP since $$\det X_n(z)$$ has no jump across $$\mathbb{R}$$ and has continuous boundary values, and so by Morera's theorem is entire. $$\det X_n(z) \to 1$$ as $$z \to \infty$$ and so by Liouville's theorem is identically $$1$$. This implies that the RHP can have at most one solution. The reader may then verify that the above is a solution.

Because $$\det X_n(z) = 1$$ we can introduce a "dual" Riemann-Hilbert problem $$\widehat{X_n}(z) = X_n(z)^{-\mathsf{T}}$$ (inverse transpose). $$\widehat{X_n}$$ solves the following RHP.

**Dual Riemann-Hilbert problem for orthogonal polynomials:**

Find a matrix valued function $$\widehat{X_n} : \mathbb{C} \setminus \mathbb{R} \to \mathbb{C}^{2 \times 2}$$ such that

1. $$\widehat{X_n}$$ is analytic (entry-wise) on $$\mathbb{C} \setminus \mathbb{R}$$.
2. $$\widehat{X_n}$$ has continuous non-tangential boundary values up to $$\mathbb{R}$$ from above ($$+$$) and below ($$-$$). We label these $$\widehat{X_n}^\pm (x) = \lim_{\epsilon \downarrow 0} \widehat{X_n}(x\pm i \epsilon)$$ for $$ x \in \mathbb{R}$$.
3. These boundary values are related by the jump condition
$$\widehat{X_n}^+(x) = \widehat{X_n}^-(x) \left( \begin{matrix} 1 & 0 \\ -w(x) & 1\end{matrix} \right)$$
4. Finally, $$\widehat{X_n}$$ is normalised at infinity by the scaling as $$z \to \infty$$

$$\widehat{X_n}(z) = \left( \mathbb{I}+\mathcal{O}(z^{-1}) \right) \left( \begin{matrix} z^{-n} & 0 \\ 0 & z^n \end{matrix} \right)$$

Indeed, we know the unique solution of the dual RHP, 

$$ \widehat{X_n}(z) = \left( \begin{matrix} -2\pi i h_{n-1}^{-1} C_\mathbb{R} \left( P_{n-1} w\right)(z)  & 2\pi i h_{n-1}^{-1} P_{n-1}(z)  \\ -C_\mathbb{R} \left( P_n w\right)(z)  & P_n(z)   \end{matrix} \right)$$

Let us now derive a pair of recursion relations. We note that $$X_{n+1}$$ satisfies properties 1-3, differing only on property 4. Thus if we let $$\Delta_n(z) = X_{n+1}(z) X_n(z)^{-1}$$ we see that $$\Delta_n$$ has no jump across the real axis, has continuous boundary values, and is analytic on $$\mathbb{C}\setminus \mathbb{R}$$. It is thus entire by Morera's theorem. Let us expand this at infinity. Let $$X_n(z) = \left( \mathbb{I}+A_n z^{-1} + \mathcal{O}(z^{-2}) \right)\left( \begin{matrix} z^n & 0 \\ 0 & z^{-n}\end{matrix} \right)$$. Then

$$\Delta_n(z)  = z  E_1 +  A_{n+1}E_1 - A_n E_1 + \mathcal{O}(z^{-1})$$

where $$E_1 = \left( \begin{matrix} 1 & 0 \\ 0 & 0\end{matrix} \right)$$. However since $$\Delta_n$$ is entire the $$\mathcal{O}(z^{-1})$$ term is identically zero, so we have

$$X_{n+1}(z) = \left( z  E_1 +  A_{n+1}E_1 - E_1 A_n  \right)X_{n}(z) $$

This is the famous "three term recurrence" for orthogonal polynomials, derived by complex analysis arguments. By a similar argument we find

$$\widehat{X_{n+1}}(z) = \left( z  E_2 -  A_{n+1}^\mathsf{T} E_2 + E_2 A_n^\mathsf{T}  \right)\widehat{X_{n}}(z) $$

where $$E_2 = \left( \begin{matrix} 0 & 0 \\ 0 & 1\end{matrix} \right)$$. Let us now consider quantity

$$Y_n(z,w) := X_n(z)^{-1}X_n(w) = \widehat{X_n}(z)^{\mathsf{T}}X_n(w)$$

Using our recursion relations for $$\widehat{X_n}$$ and $$X_n$$ we may relate $$Y_n$$ and $$Y_{n+1}$$ by $$Y_{n+1}(z) = X_n^{-1}(z) \Delta_n(z)^{-1} \Delta_n(w) X_n(w)$$. $$\Delta_n(z)^{-1} \Delta_n(w)$$ is a polynomial in two variables, moreover

$$\Delta_n(z)^{-1} \Delta_n(w) = (z-w) (A_{n+1})_{21} E_{21} + \mathbb{I}$$

where $$E_{21} = \left( \begin{matrix} 0 & 0 \\ 1 & 0\end{matrix} \right)$$. By our formula for solution $$X_n$$ we find $$(A_{n+1})_{21} = - 2\pi i h_n^{-1}$$, and so

$$Y_{n+1}(z,w) = Y_n(z,w) - \frac{2\pi i}{h_n} (z-w) X_n^{-1}(z) E_{21} X_n(w)$$

If we now take the $$(2,1)$$ matrix element of both sides we find

$$Y_{n+1}(z,w)_{21} = Y_n(z,w)_{21} - 2\pi i \frac{P_n(z)P_n(w)}{h_n} (z-w) $$

If we now sum both sides from $$n=0$$ to $$n=N-1$$ we find

$$Y_{N+1}(z,w)_{21} = Y_0(z,w)_{21} - 2\pi i (z-w) \sum_{n=0}^{N-1} \frac{P_n(z)P_n(w)}{h_n}  $$

Then from our solution $$Y_0(z,w)_{21} = 0$$ we find 

$$\boxed{ \sum_{n=0}^{N-1} \frac{P_n(z)P_n(w)}{h_n} = -\frac{1}{2\pi  i} \frac{(X_N(z)^{-1}X_N(w))_{21}}{z-w}}$$

which is the Christoffel-Darboux formula.

**Remark:** This way of writing the Christoffel-Darboux formula mirrors nicely with what happens for $$\beta = 4$$. Here the relevant quantity that encodes eigenvalue correlation functions is the "pre-kernel," written as a sum over *skew-orthogonal polynomials*. Namely,

$$\sum_{k=0}^{n-1} \frac{ P_{2k}(x) e^{-V(x)} \frac{d}{dy}\left( P_{2k+1}(y) e^{-V(y)} \right) - P_{2k+1}(x) e^{-V(x)} \frac{d}{dy}\left( P_{2k}(y) e^{-V(y)} \right)}{2 h_k} $$
$$= - \frac{e^{-V(x)-V(y)}}{4\pi i} \frac{(A_n(x)^{-1}A_n(y))_{21}}{x-y}$$

where $$P_k$$ is the $$k$$th monic skew-orthogonal polynomial, $$h_k$$ is the skew-norm, and $$A_n$$ is a Riemann-Hilbert problem introduced in my [recent paper](https://arxiv.org/abs/2306.14107).
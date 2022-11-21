---
layout: post
title:  "Estimating generalised hypergeometric functions"
date:   2022-08-14 18:50:11 +0100
permalink: /hypergeometric/
katex: True
---

In my [work](https://arxiv.org/abs/2102.08842) on products of truncated orthogonal matrices (in collaboration [N. Simm](https://profiles.sussex.ac.uk/p435611-nicholas-simm) and [F. Mezzadri](https://www.bristol.ac.uk/people/person/Francesco-Mezzadri-66ca5240-8f45-4ffc-a838-d1f68827bd23/)) it became important to estimate the following function for $$x \in (-1,1)$$ for $$L,N \to \infty$$.

$$f_{N-2,L}(x) = \sum_{k=0}^{N-2} \binom{L+k}{k}^m x^k $$

where $$m \in \mathbb{N}$$ is fixed. In this post I will show how this estimate is carried out. In my opinion this is the key estimate of the paper so I would like to draw attention to it. This method can be used for the case $$N = +\infty$$ and so can be used to estimate the generalised hypergeometric function

$$f_{\infty,L}(x) = \sum_{k=0}^{\infty} \binom{L+k}{k}^m x^k = {}_m F_{m-1}\left( \begin{matrix} L+1 & \dots & L+1 \\
  1 & \dots & 1 \end{matrix} \, \bigg| \, x \right)$$

<h2>Background</h2>

Before discussing this, let me briefly outline the context of the problem. Let

$$U_1, \dots, U_m \in O(N+L)$$

be $$m$$ independently sampled matrices from the orthogonal group $$O(N+L)$$ according to Haar measure. We call the $$N \times N$$ upper left corner of $$\tilde{U}_i$$ a *truncated orthogonal matrix*. $$\tilde{U}_i$$ is thus a random matrix with real matrix elements, whose randomness is inherited from Haar measure.

We are interested in the spectrum of the product

$$
X = \tilde{U}_1 \tilde{U}_2 \dots \tilde{U}_m
$$

as $$N,L \to \infty$$ and where $$\frac{L}{N} \to \gamma > 0$$. Define the *real spectral density* as the function $$\rho_{N,L} : \mathbb{R} \to [0, \infty)$$ such that

$$\mathbb{E}[| \sigma(X) \cap A |] = \int_A \rho_{N,L}(x) \, dx $$

where $$\sigma(X)$$ is the spectrum of $$X$$, $$A \subset \mathbb{R}$$ is a Lebesgue measurable set and $$\lvert \cdot \rvert$$ denotes cardinality. From this we see that

$$\mathbb{E}[ |\sigma(X) \cap \mathbb{R}|] = \int_\mathbb{R} \rho_{N,L}(x) \, dx$$

gives the expected number of real eigenvalues. [Forrester, Ipsen and Kumar](https://arxiv.org/abs/1708.00967) supply the following exact formula for the real spectral density.

$$\rho_{N,L}(x)= \int_{[-1,1]} |x-y| w_L(x) w_L(y) f_{N-2,L}(xy) \, dy $$

$$w_L$$ is the so-called "weight function," which we will not write out and can be found in our paper. To estimate $$\rho_{N,L}$$ we wish to estimate both $$w_L$$ and $$f_{N-2,L}$$. The former turns out to be a straightforward application of the Laplace method; it is less obvious how to carry out the latter and is the subject of this post.

<h2>The estimate</h2>

**Lemma:** Let $$g_{K}(z) = \sum_{k=0}^{K} a_k z^k $$ for $$K \in \mathbb{N}$$ and suppose $$\lim_{K \to \infty} g_{K}(z) = g_\infty(z)$$ converges on some neighbourhood $$U$$ of $$0 \in \mathbb{C}$$.

Then

$$\sum_{k=0}^K a_k^m x^k = \frac{1}{(2\pi i)^{m-1}} \oint_{\Gamma^{m-1}} g_{K}\left( \frac{x}{z_1 \dots z_{m-1}}\right) g_{\infty} (z_1) \dots g_{\infty}(z_{m-1}) \frac{dz_1}{z_1} \dots \frac{dz_{m-1}}{z_{m-1}}$$

where $$\Gamma \subset U \setminus \{ 0 \}$$ is a closed contour enclosing $$0$$. This formula is also valid for $$K = +\infty$$ so long as $$\frac{x}{z_1 \dots z_{m-1}} \in U$$ throughout the contour $$\Gamma$$.

**Proof:** If we expand the product

$$g_{K}\left( \frac{x}{z_1 \dots z_{m-1}}\right) g_{\infty} (z_1) \dots g_{\infty}(z_{m-1}) = \sum_{k_1, \dots k_{m-1}=0}^{\infty} \sum_{n=0}^K a_n a_{k_1} \dots a_{k_{m-1}} x^n z_1^{k_1 - n} \dots z_{m-1}^{k_{m-1} - n}$$

Clearly $$\sum_{k=0}^K a_k^m x^k$$ is the coefficient of $$z_1^0 \dots z_{m-1}^0$$ in the above series, which can be picked out by the residue theorem. The above series is uniformly convergent on compact sets within the radius of convergence so term by term integration is justified.  $$\square $$

This means so long as we have good estimates on the case of $$m=1$$ we can extract good estimates in the case of general $$m$$.

**Remark:** Let $$f$$ and $$g$$ be two analytic functions defined in a neighbourhood of $$0$$. Define the convolution

$$(f \ast g)(x) = \frac{1}{2\pi i} \oint_\Gamma f(z) g\left( \frac{x}{z} \right) \frac{dz}{z}$$

where $$\Gamma$$ is a positively oriented contour that encloses $$0$$. Then our above lemma states that

$$\sum_{k=0}^K a_k^m x^k = g_K \ast \underbrace{g_\infty \ast \dots \ast g_\infty}_{m-1 \text{ times}} (x)$$

One thus sees that our above lemma is nothing other than an instance of the convolution theorem (sometimes going under the name of [Hadamard products](https://en.wikipedia.org/wiki/Generating_function_transformation#Hadamard_products_and_diagonal_generating_functions)). $$\triangle$$

The following is well known but we include a proof for completeness.

**Lemma:** $$\sum_{k=0}^\infty \binom{L+k}{k} x^k = \frac{1}{(1-x)^{L+1}} $$

**Proof:** Using the Cauchy residue theorem write

$$\binom{L+k}{k} = \frac{1}{2\pi i} \oint_\Gamma \frac{(1+z)^{L+k}}{z^{k+1}} \, dz$$
where $$\Gamma$$ is a positively oriented contour enclosing $$0$$. Then

$$\sum_{k=0}^\infty \binom{L+k}{k} x^k =  \frac{1}{2\pi i} \oint_\Gamma \frac{(1+z)^{L}}{z} \sum_{k=0}^\infty \left( \frac{x(1+z)}{z} \right)^k \, dz =  \frac{1}{1-x}\frac{1}{2\pi i} \oint_\Gamma (1+z)^{L} \frac{1}{z- \frac{x}{1-x}} \, dz $$

where $$x$$ is chosen sufficiently small that $$\left\lvert \frac{x(1+z)}{z} \right\rvert < 1$$ on the contour. This implies that the pole at $$z = \frac{x}{1-x}$$ is enclosed. $$\square $$

This immediately gives a formula for the $$N = +\infty$$ case,

$$f_{\infty,L}(x) = \frac{1}{(2\pi i )^{m-1}} \oint_{\Gamma^{m-1}} \frac{1}{\left( 1 - \frac{x}{z_1 \dots z_{m-1}} \right)^{L+1}} \prod_{k=1}^{m-1}\frac{dz_k}{(1-z_k)^{L+1} z_k} $$

Applying the method of steepest descent allows one to immediately obtain $$L \to +\infty$$ asymptotics of $$f_{\infty,L}$$.

**Proposition:** Let $$x \in (0,1)$$. Then we have the following asymptotics pointwise as $$L \to +\infty$$

$$f_{\infty,L}(x) \sim \frac{1}{\sqrt{m} (2\pi L)^{\frac{m-1}{2}}} \frac{1}{x^\frac{m-1}{2m} \left( 1 - x^\frac{1}{m} \right)^{mL+1}}$$

We also have the following bounds for $$x \in [0,1)$$

$$|f_{\infty,L}(x)| \leq \left( \frac{\pi}{2L} \right)^{\frac{m-1}{2}} \frac{1}{ x^\frac{m-1}{2m} \left( 1 - x^\frac{1}{m} \right)^{m(L+1)} }$$

$$|f_{\infty,L}(-x)| \leq  \frac{1}{  \left( 1 - x^\frac{1}{m} \right)^{m(L+1)}} e^{-\frac{Lx^\frac{1}{m}}{2m}}$$

**Proof:** The first follows from an application of the steepest descent method, the second follows from the inequality in equation 6.47 of our paper, the third is equation 2.16 of our paper. $$\square$$

Notice that the second estimate is quite good, it differs from the pointwise asymptotics by a $$\mathcal{O}(1)$$ factor.

Let $$g_{N-2}(x) = \sum_{k=0}^{N-2} \binom{L+k}{k}x^k$$. There are a variety integral representations of this, e.g. in terms of an incomplete beta function (see page 3 of [Khoruzhenko, Sommers and Zyczkowski](https://arxiv.org/abs/1008.2075)). If we write the coefficient $$\binom{L+k}{k} = \frac{1}{2\pi i } \oint_{\Gamma} \frac{1}{z^{N-1}(1-z)^{L+1}} \, dz $$ and sum, we find

$$g_{N-2}(x) = \frac{1}{(1-x)^{L+1}} \chi_{R> |x|} -\frac{x^{N-1 }}{2\pi i} \oint_{|z|=R} \frac{1}{z^{N-1}(1-z)^{L+1}} \frac{dz}{z-x} $$

for any $$R > 0$$. A calculation shows that the steepest descent contour for the integral contained in the second term is $$R = \frac{1}{1+\gamma}$$. Putting this all together yields an integral represention of $$f_{N-2,L}$$.

**Remark:** The technique discussed in this post can also be used to study the asymptotics of other generalised hypergeometric functions. For example, it allows one to obtain asymptotics of

$$\sum_{k=0}^\infty \frac{x^k}{(k!)^m} = {}_0 F_{m-1}\left( \begin{matrix}  & - &  \\
  1 & \dots & 1 \end{matrix} \, \bigg| \, x \right)$$

for any fixed $$m \in \mathbb{N}$$ in the régime $$x \to \infty$$ in any direction in the complex plane. $$\triangle$$

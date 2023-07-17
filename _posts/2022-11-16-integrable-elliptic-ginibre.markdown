---
layout: post
title:  "Integrable structure in the elliptic Ginibre ensemble"
date:   2022-11-16 12:00:11 +0100
permalink: /integrable-elliptic-ginibre/
katex: True
---

This post will provide an accessible introduction to work contained in my [joint paper](https://arxiv.org/abs/2208.04684) with Thomas Bothner. In this work we find integrable structure in the elliptic Ginibre ensemble and find interesting connections to the theory of *finite temperature* Airy processes. Our paper also performs a nonlinear steepest descent analysis of the cumulative distribution function of the rightmost particle, but this is left out of this post to keep the discussion nontechnical.

<h2>Background and motivation</h2>

Suppose you have an ecosystem with $$n$$ species, with populations $$x_1, \dots, x_n$$. Is this ecosystem stable? This question was raised by Robert May in his article ["Will a Large Complex System be Stable?"](https://www.gwern.net/docs/sociology/1972-may.pdf) (1972). Suppose as a first approximation that each species is isolated and the only thing limiting population is competition over some non-organic resources. You could model this by the system of uncoupled ODEs

$$\dot{x}_i = - \mu x_i$$

for $$\mu > 0$$. We have shifted the populations so that $$x_i = 0$$ is the equilibrium. Now if we wanted to model interactions between species we could write

$$\dot{x}_i = - \mu x_i + \sum_{j=1}^n A_{ij} x_j$$

where $$A$$ is some interaction matrix. $$A_{ij}$$ represents the effect of the population of species $$j$$ on the growth rate of species $$i$$. (Note that there is no reason to think this matrix is symmetric.) This could be written in vector form

$$\dot{\mathbf{x}} = (- \mu \mathbb{I} + A)\mathbf{x} $$

It is clear that this system is stable when $$\mu > \max_i \Re \lambda_i$$ and unstable when $$\mu < \max_i \Re \lambda_i$$, where $$\lambda_1, \dots, \lambda_n$$ are the eigenvalues of $$A$$. Given that we expect $$n$$ to be very large and the interactions between the species very complex, we cannot expect to model $$A$$ exactly. Instead we take $$A$$ to be random. Note that the probability of stability, $$\mathbb{P}(\max_i \Re \lambda_i < \mu)$$, is exactly the cumulative distribution function for $$\max_i \Re \lambda_i$$.


In our work we are interested in $$A \in \mathcal{M}_n(\mathbb{C})$$, so it models stability of a linear system over the field of *complex* numbers. I will list three candidate distribution functions for $$A$$. This list is by no means exhaustive -- you can give many others -- however these three cases are "exactly solvable" in that one can compute their correlation functions explicitly.

<h3>The (complex) Ginibre Ensemble</h3>

This ensemble was introduced by [Ginibre in 1965](https://aip.scitation.org/doi/10.1063/1.1704292). This is an ensemble of $$n \times n$$ matrices with independent, identically distributed matrix elements, each given by complex Gaussians. One could write the density as

$$P_{\mathrm{Gin}}(X) \, dX = \pi^{-n^2} e^{-\mathrm{tr}\, X X^\ast} \, d X$$

The eigenvalues have the joint distribution

$$\varrho_n(z_1, \dots, z_n) = C_n e^{-\sum_{k=1}^n |z_k|^2} \prod_{1 \leq i < j \leq n } |z_i - z_j|^2 $$

The Ginibre ensemble obeys a "circular law." If we define the 1-point density as

$$\varrho(x) = \frac{1}{n} \int_{(\mathbb{C})^{n-1}} \varrho_n(x, z_2, \dots, z_n) \, d^2 z_2 \dots d^2 z_n$$

then $$\varrho(\sqrt{n}x ) \to \frac{1}{\pi} \chi_{\lvert x \rvert < 1}$$ as $$n \to\infty$$ in the weak-$$\ast$$ sense. That is, the spectral density tends towards the unit disk. Furthermore, the eigenvalue with largest real part is asymptotically Gumbel distributed. More precisely

$$\mathbb{P}\left( \max_i \Re \lambda_i \leq \sqrt{n} + \sqrt{\frac{\gamma_n}{4}} + \frac{t}{\sqrt{4 \gamma_n}} \right) \to e^{-e^{-t}}$$

as $$n \to \infty$$, where

$$\gamma_n = \frac{1}{2}(\ln n - 5 \ln \ln n - \ln (2\pi^4))$$

This result can be found in [G. Cipolloni, L. Erdős, D. Schröder, and Y. Xu (2022)](https://arxiv.org/abs/2206.04443). Furthermore, we show in work forthcoming on the arXiv, that the real parts of the eigenvalues converge locally, in the bulk of the spectrum, to a Poisson point process as $$n\to \infty$$.

<h3>The Gaussian Unitary Ensemble</h3>

The Gaussian Unitary Ensemble (GUE) is an ensemble of $$n \times n$$ complex Hermitian random matrices. The "unitary" in the name comes from the fact that the ensemble has a unitary symmetry. The density is formally very similar to the Ginibre case

$$P_{\mathrm{GUE}}(X) \, dX = C^\prime_n e^{-\mathrm{tr}(X^2)} \, dX$$

The difference is that this density is restricted to the subspace of matrices such that $$X=X^\ast$$. In this case all the eigenvalues lie on $$\mathbb{R}$$ and their joint density is given by

$$\varrho_n(x_1, \dots , x_n) = C^{\prime \prime}_n e^{-\sum_{k=1}^n x_i^2} \prod_{1\leq i < j \leq n }|x_i - x_j|^2 $$

This looks practically identical to the case of Ginibre, but in fact the restriction to $$\mathbb{R}$$ rather than $$\mathbb{C}$$ changes things considerably. As is well known, the 1-point density [converges to a semicircular distribution](https://mathworld.wolfram.com/WignersSemicircleLaw.html). Furthermore, its rightmost eigenvalue asymptotically obeys a [Tracy-Widom (1994)](https://arxiv.org/abs/hep-th/9211141) law. That is

$$\mathbb{P}\left( \max_{i=1, \dots, n} \lambda_i \leq \sqrt{2n} + \frac{t}{\sqrt{2} n^\frac{1}{6}} \right) \to \exp\left(-\int_t^\infty (s-t) q(s)^2 \, ds\right)$$

as $$n \to \infty$$, and where $$q$$ satisfies Painlevé II,

$$q^{\prime \prime}(t) = t q(t) + 2 q(t)^3$$

with boundary condition $$q(t) \sim \mathrm{Ai}(t)$$ as $$t \to +\infty$$. This Tracy-Widom distribution is a highly universal object appearing in many different contexts such as random matrix theory, Ulam's problem, the KPZ equation, asymmetric exclusion processes, Aztec diamonds etc. This result is also celebrated because it connects random matrix theory to integrable systems (though it was not the first to do so).

<h3>The Elliptic Ginibre Ensemble</h3>

We start with the following observation. A Ginibre matrix can be generated by

$$X = \frac{1}{\sqrt{2}}(H_1 + i H_2)$$

where $$H_1, H_2$$ are independently sampled GUE matrices. Motivated by this, consider the random matrix

$$X =  \sqrt{\frac{1+\tau}{2}}H_1 + i \sqrt{\frac{1-\tau}{2}} H_2$$

for $$\tau \in [0,1]$$. $$\tau = 0$$ yields the Ginibre ensemble, $$\tau = 1$$ yields the GUE. Thus we can interpolate between these two ensembles. This is called the "elliptic Ginibre ensemble." Its density is given by

$$P_{\mathrm{eGin}}(X)\, dX = C_n^{\prime \prime \prime} e^{- \frac{1}{1-\tau^2}\mathrm{tr}\, (X X^\ast - \tau \Re (X^2))} \, dX$$

The ensemble derives its name because, up to a $$\sqrt{n}$$ scaling, the spectral density tends (as $$n \to \infty$$) towards a constant on the ellipse

$$\left\{ (x,y) \in \mathbb{R}^2 \, : \frac{x^2}{(1+\tau)^2} + \frac{y^2}{(1-\tau)^2} < 1  \right\} \subset \mathbb{R}^2 \simeq \mathbb{C}$$

and zero outside.

In order to see something interesting for the extremal eigenvalue we need to take a double scaling limit, where $$\tau_n = 1 - \frac{\sigma^2}{n^\frac{1}{3}}$$, for $$\sigma > 0$$ a fixed parameter. This limit where the ellipse is very "flat" is called "weak non-Hermicity." We are interested in the distribution of the largest real part in the limit. This scaling (due to [Bender, 2010](https://link.springer.com/article/10.1007/s00440-009-0207-9)) is carefully chosen. Let $$\tau_n$$ tend to $$1$$ too fast and local correlations look like the GUE, too slow and local correlations look like Ginibre.

<h2>Results</h2>

[Bender (2010)](https://link.springer.com/article/10.1007/s00440-009-0207-9) defines a rescaled point process at the edge of the ellipse. One zooms in on the edge of ellipse in a specified way.

$$x_j := \frac{\Re \lambda_j - c_n}{a_n}$$

$$y_j := \frac{\Im \lambda_j}{b_n}$$

The scalings $$a_n, b_n, c_n$$ can be found in Bender's paper, and I will just state the result: that the rescaled point process $$\{ z_j := (x_j , y_j) \}_{j=1}^n$$ converges as $$n\to \infty$$. This yields a determinantal point process parametrised by $$\sigma > 0$$. We now arrive at our main result.

**Theorem (Bothner-L.):**

$$\boxed{F_\sigma(t) := \lim_{n \to \infty} \mathbb{P}\left( \max_{j=1,\dots, n} x_j \leq t \right) = \exp\left( - \int_t^\infty (s-t) \left[ \int_\mathbb{R} p_\sigma(s,y)^2 \, d\nu_\sigma(y) \right] ds \right)}$$

where $$d\nu_\sigma(\lambda) = \frac{1}{\sigma \sqrt{\pi}} e^{-\left( \frac{\lambda}{\sigma}\right)^2}\, d\lambda$$ and where $$p_\sigma$$ satifies *integro-differential Painlevé II*,

$$ \boxed{\frac{\partial^2}{\partial t^2} p_\sigma(t,y) = \left[ t+y+2\int_\mathbb{R} p_\sigma(t,\lambda)^2 \, d \nu_\sigma(\lambda) \right] p_\sigma(t,y)}$$

with boundary condition $$p_\sigma(t,y) \sim \mathrm{Ai}(t+y)$$ as $$t \to +\infty$$, pointwise in $$y \in \mathbb{R}$$. $$\triangle$$

How does one see that this generalises the Tracy-Widom result? If one takes $$\sigma \downarrow 0$$ one expects to reduce to the GUE. One sees that $$d\nu_\sigma(\lambda)  \to \delta(\lambda) \, d\lambda$$, and so

$$\lim_{\sigma \downarrow 0} F_\sigma(t) = \exp\left( - \int_t^\infty (s-t) p_\sigma(s,0)^2 \, ds \right)$$

and our integro-differential equation, when evaluated at $$y =0$$, reduces to Painlevé II, with the right boundary condition.

On the other hand, recovering the Gumbel distribution as $$\sigma \to \infty$$ is actually harder and is done in our paper. This is best done by studying the asymptotics of the corresponding Fredholm determinant rather than working with the integro-differential equation. We also look at asymptotics of $$F_\sigma(t)$$ under various scalings of $$t$$ and $$\sigma$$ by means of nonlinear steepest descent techniques.

**Remark:** There is work left to be done here, especially for the left tail $$t \to -\infty$$, since we only obtain asymptotics in a scaling régime where $$\sigma$$ is very small.

<h2>Sketch of proof</h2>

As the point process defined by Bender is determinantal, gap probabilities may be expressed in terms of Fredholm determinants.

$$F_\sigma(t) = \det(1- K^{\sigma}_{\mathrm{Ai}})_{L^2((t,\infty)\times \mathbb{R})}$$

The kernel of $$K^{\sigma}_{\mathrm{Ai}}$$ (found by [Bender, 2010](https://link.springer.com/article/10.1007/s00440-009-0207-9)) is complicated and is given in Equation 1.10 of [our paper](https://arxiv.org/abs/2208.04684). We make the following observation. Let $$\mathcal{F} : L^2(\mathbb{R}_+ \times \mathbb{R}) \to L^2(\mathbb{R}_+ \times \mathbb{R})$$

$$\mathcal{F}f(x,\eta) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^\infty e^{-i y \eta} f(x,y) \, dy$$

be the Fourier transform in the vertical (imaginary) direction. Then we have

$$K_\sigma = \mathcal{F} K^{\sigma}_{\mathrm{Ai}} \mathcal{F}^{-1}$$

where $$K_{\sigma}$$ is an integral operator with kernel

$$K_{\sigma}(z_1, z_2) = \int_0^\infty \phi_\sigma(x_1 + s, y_1 ) \phi_\sigma(x_2 + s, y_2 ) \, ds$$

where $$\phi(x,y) = \pi^{-\frac{1}{4}} e^{- \frac{y^2}{2}}\mathrm{Ai}(x+\sigma y)$$ and $$z_i \equiv (x_i, y_i)$$. By the invariance of the determinant, we then have

$$F_\sigma(t) = \det(1- K_{\sigma})_{L^2((t,\infty)\times \mathbb{R})}$$

Notice that $$K_{\sigma}(z_1, z_2)$$ takes the form of a Hankel composition, and so one can use the method given in [my previous post](https://almaths.github.io/maths_blog/tracy-widom/) to derive a Tracy-Widom-like result.

<h3>Finite temperature Airy processes</h3>

Our result is curious because the integro-differential equation we derive also appears in the study of *finite temperature Airy processes* on $$\mathbb{R}$$, whereas we are looking at a point process in $$\mathbb{R}^2$$. What is the relationship?

Begin by transferring the $$t$$ dependence to the operator.

$$F_\sigma(t) = \det(1- K_{\sigma}^t)_{L^2(\mathbb{R}_+\times \mathbb{R})}$$

where

$$K_{\sigma}^t((x_1, y_1),(x_2, y_2)) = K_{\sigma}((x_1+t, y_1),(x_2+t, y_2))$$

Now define the operator $$P_{t,\sigma} : L^2(\mathbb{R}_+ \times \mathbb{R}) \to L^2(\mathbb{R}_+)$$ with kernel

$$P_{t,\sigma}(s, (x,y)) = \pi^{-\frac{1}{4}} e^{-\frac{1}{2}y^2} \mathrm{Ai}(s+x+\sigma y + t)$$

Then a short calculation shows that

$$P_{t,\sigma}^\ast P_{t,\sigma} = K_\sigma^t$$

Now consider $$P_{t,\sigma} P_{t,\sigma}^\ast : L^2(\mathbb{R}_+) \to L^2(\mathbb{R}_+)$$. A short calculation shows that this operator has kernel

$$P_{t,\sigma} P_{t,\sigma}^\ast(s_1, s_2)=\int_\mathbb{R} \Phi\left( \frac{y}{\sigma} \right) \mathrm{Ai}(s_1+y+t) \mathrm{Ai}(s_2+y+t) \, dy =: N_\sigma^t(s_1, s_2)$$

where $$\Phi(x) = \frac{1}{\sqrt{\pi}}\int_{-\infty}^x e^{-t^2}\, dt$$. Undoing the shift by letting $$N_\sigma (s_1 +t , s_2 + t) = N_\sigma^t (s_1  , s_2)$$ and using Sylvester's identity we find

$$F_\sigma(t) = \det(1-N_\sigma)_{L^2(t,\infty)}$$

This means that the largest particle of the finite temperature Airy process with function $$\Phi$$ is identically distributed to the largest *real* part in the elliptic edge process.

**Remark:** Indeed, this observation may be generalised. Let us consider any finite temperature Airy process, with kernel

$$N(x,y) = \int_\mathbb{R} \psi(s) \mathrm{Ai}(x+s) \mathrm{Ai}(y+s) \, ds$$

where $$\psi(x) \geq 0$$ is an increasing, continuously differentiable function such that $$\psi(-\infty) = 0$$ and $$\psi(+\infty) = 1$$. By an identical argument to that presented above we have the identity

$$\det(1-N)_{L^2(t,\infty)} = \det(1-K)_{L^2((t,\infty) \times \mathbb{R})}$$

where

$$K((x_1,y_1), (x_2,y_2)) = \int_0^\infty \phi(x_1 +s, y_1) \phi(x_2 +s, y_2) \, ds$$

where $$\phi(x,y) = \sqrt{\psi^\prime(y)}\mathrm{Ai}(x+y)$$. Thus $$K$$ has the form of a Hankel composition and so the Fredholm determinant has an associated integro-differential representation (by the method of [my previous post](https://almaths.github.io/maths_blog/tracy-widom/)). That such Fredholm determinants have an integro-differential representation is already known in the literature but, to the author's knowledge, this method of relating finite temperature kernels to Hankel compositions is new. $$\triangle$$


**Remark:** Here we make a remark that will be well known to experts, namely the relationship between finite temperature kernels and Riemann-Hilbert problems. It is clear that if we let

$$A : L^2(\mathbb{R}) \to L^2(t,\infty)$$

with kernel

$$A(x,x^\prime) = \mathrm{Ai}(x+x^\prime) \sqrt{\psi(x^\prime)}$$

Then clearly $$A A^\ast = N : L^2(t,\infty) \to L^2(t,\infty)$$. A simple calculation shows that $$M_t := A^\ast A : L^2(\mathbb{R}) \to L^2(\mathbb{R})$$ has kernel

$$M_t (x,y) = \sqrt{\psi(x)}\sqrt{\psi(y)} \frac{\mathrm{Ai}(x+t)\mathrm{Ai}^\prime (y+t) - \mathrm{Ai}^\prime (x+t)\mathrm{Ai}(y+t)}{x-y}$$

This is an integrable Its-Izergin-Korepin-Slavnov type operator (see this [nice introduction](https://www.ams.org/books/trans2/189/) to integrable operators by Deift). The IIKS theory relates Fredholm determinants of such operators to Riemann-Hilbert problems. Thus

$$\det(1-N)_{L^2(t,\infty)}= \det(1-M_t)_{L^2(\mathbb{R})}$$

may be related to a Riemann-Hilbert problem. This then allows for a Deift-Zhou steepest descent analysis.

**Remark:** As a final remark, these factorisations are useful at proving that the corresponding operator is trace class, since one need only show that the factors are respectively Hilbert-Schmidt.

<h3>The bulk of the spectrum (addendum)</h3>

The techniques developed above can be used to derive a similar integro-differential representation for gaps between real parts in the *bulk* of the spectrum (the reader is referred to [our paper on the arXiv](https://arxiv.org/abs/2212.00525) for details). In particular, the correlation kernel in the bulk also satisfies a curious factorisation property. Recall that we are considering the régime where $$\tau_n \to 1$$ and so "the bulk" is approximately the set $$(-2,2)$$.

"Weak non-Hermiticity" in the bulk requires a somewhat different scaling than at the edge. In particular, if $$\lambda_0 \in (-2,2)$$ is the point we're going to "zoom in" on and $$\rho_1(x) = \frac{1}{\pi}\sqrt{\left(1-\frac{x^2}{4}\right)_+}$$, then we must take

$$\tau_n = 1- \frac{1}{n} \left( \frac{\sigma}{\rho_1(\lambda_0)}\right)^2$$

for $$\sigma \geq 0$$. Then, after a suitable rescaling, the point process around the point $$\lambda_0 \in (-2,2)$$ converges to a determinantal point process with kernel given by

$$K_{\mathrm{sin}}^{\sigma}(z_1,z_2) = \frac{1}{\sigma\sqrt{\pi}} e^{-\frac{y_1^2 + y_2^2}{2\sigma^2}} \frac{1}{2\pi}\int_{-\pi}^\pi e^{-(\sigma u)^2} \cos(u(z_1 - \overline{z_2})) \, du$$

where $$y_k = \Im z_k$$ for $$k=1,2$$. This kernel describes a determinantal point process in the plane $$\mathbb{R}^2 \simeq \mathbb{C}$$.

Now let $$J_t = (-t,t)\times \mathbb{R} \subset \mathbb{R}^2$$ and suppose we are interested in the gap probability given by

$$\det(1-K_{\mathrm{sin}}^{\sigma})_{L^2(J_t)}$$

This corresponds to looking at gaps between real parts (note also that the point process is horizontally translation invariant). We now observe the following factorisation. Let

$$A_\sigma : L^2(J_t) \to L^2(-\pi,\pi)$$

with kernel

$$A_\sigma(a,z) = \frac{1}{\sqrt{2\pi^\frac{3}{2}\sigma}} \exp\left(-\frac{y^2}{2\sigma^2} - \frac{1}{2}(\sigma a)^2 - ia \overline{z}\right)$$

A simple calculation then shows that

$$K_{\mathrm{sin}}^{\sigma} = A_\sigma^\ast A_\sigma : L^2(J_t) \to L^2(J_t)$$

But then, by Sylvester's identity, we have

$$\det(1-K_{\mathrm{sin}}^{\sigma})_{L^2(J_t)} = \det(1-A_\sigma A_\sigma^\ast)_{L^2(-\pi,\pi)} = \det(1-S_\sigma^t )_{L^2(-t,t)}$$

where $$S_\sigma^t$$ is a rescaled version of $$A_\sigma A_\sigma^\ast$$ (so that its domain is now $$L^2(-t,t)$$). A simple calculation (see our paper) shows that one can write the kernel of $$S_\sigma^t$$ as

$$S_\sigma^t(a,b) = \int_0^\infty \left[ \Phi\left( \frac{t}{\sigma}(z+1) \right) - \Phi\left( \frac{t}{\sigma}(z-1)\right) \right]\cos(\pi(a-b)z) \, dz$$

This is exactly a *finite temperature sine kernel*. In our paper we show that for any determinantal point process on $$\mathbb{R}$$ with kernel of the form

$$K(a,b) = \int_0^\infty w(z)\cos(\pi(a-b)z) \, dz$$

for $$w: \mathbb{R}_+ \to [0,1)$$ and tending to $$0$$ exponentially fast at $$+\infty$$, the gap probability $$\det(1-K)_{L^2(-t,t)}$$ can be represented in terms of a solution to an *integro-differential Painlevé V equation*. The generalises the famous result of Jimbo-Miwa-Môri-Sato (1980).

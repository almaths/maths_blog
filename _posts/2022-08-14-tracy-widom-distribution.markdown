---
layout: post
title:  "A derivation of the Tracy-Widom distribution by the Hankel composition method"
date:   2022-09-08 18:50:11 +0100
categories: jekyll update
permalink: /tracy-widom/
katex: True
---

A celebrated result in the theory of random matrices is the connection between the extremal eigenvalue of a random matrix sampled from the Gaussian Unitary Ensemble and the Painlevé II equation. In this post I will give a derivation of this that uses the Hankel composition method. This method is given in great generality in the [paper](https://arxiv.org/abs/2205.15007) of Bothner and was used in our [joint work](https://arxiv.org/abs/2208.04684) to relate the distribution function of the extremal eigenvalue in the elliptic Ginibre ensemble to an integro-differential equation. By showing how this method works in the case of the GUE it should shed some light on how it works in other cases.

<h2>Background</h2>

The Gaussian Unitary Ensemble (GUE) is an ensemble of $$n \times n$$ Hermitian random matrices with probability density

$$\frac{1}{Z_{\mathrm{GUE}}} e^{- \frac{1}{2} \mathrm{tr}(H^2)} $$

$$Z_{\mathrm{GUE}}$$ is a normalisation constant. We are interested in the distribution of the extremal (rightmost) eigenvalue. A famous result (see Chapter 24 of Mehta's *Random Matrices*) shows that the cumulative distribution converges, under an appropriate scaling, to the Fredholm determinant of the *Airy kernel.* Let $$\lambda_n$$ be the rightmost eigenvalue.

$$F(t) \equiv \lim_{n \to \infty} \mathbb{P}\left(\lambda_n \leq \sqrt{2n} + \frac{t}{\sqrt{2}n^\frac{1}{6}}\right) = \det(1 - K)_{L^2(t,\infty)}$$

where $$K: L^2(t,\infty) \to L^2(t,\infty)$$ is the operator with kernel

$$K(x,y) = \frac{\mathrm{Ai}(x) \mathrm{Ai}^\prime(y) - \mathrm{Ai}^\prime(x)\mathrm{Ai}(y) }{x-y} = \int_0^\infty \mathrm{Ai}(x+s)\mathrm{Ai}(y+s) \, ds$$

The motivation for studying this is not simply that the GUE is an easy model to study, but also that this Airy kernel is *universal* (see Deift's *Orthogonal Polynomials and Random Matrices: A Riemann-Hilbert Approach*). That is, suppose we have an ensemble of Hermitian matrices with probability density

$$\frac{1}{Z_{V}} e^{- n \mathrm{tr} V(H)} $$

for some entire function $$V$$ which grows sufficiently rapidly at $$\pm \infty$$, e.g. a polynomial. The eigenvalues will asymptotically ($$n \to \infty$$) concentrate on disjoint intervals $$[\alpha_1, \beta_1], \dots , [\alpha_m, \beta_m]$$. The distribution of the extremal eigenvalue at these endpoints $$\alpha_1, \beta_1, \dots, \alpha_m , \beta_m$$ will converge after a suitable rescaling to $$\det(1 - K)_{L^2(t,\infty)}$$. There is a similar universality in the bulk where the "universal" kernel is the sine kernel. The latter was actually discovered first by Jimbo, Miwa, Môri and Sato in 1980 (see [here](https://core.ac.uk/download/pdf/25350076.pdf) for an accessible introduction to this work).

<h2>The Connection to Painlevé II</h2>

**Theorem** ([Tracy and Widom 1993](https://arxiv.org/abs/hep-th/9210074))**:**

$$\boxed{F(t) = \exp\left(-\int_t^\infty (s-t)q(s)^2 \, ds\right)}$$

where $$q$$ solves Painlevé II

$$\boxed{q^{\prime \prime}(t) = t q(t)+ 2 q(t)^3}$$

and we have the boundary condition $$q(t) \sim \mathrm{Ai}(t)$$ as $$t \to +\infty$$. $$\triangle$$

We demonstrate the above by showing that

$$\frac{d^2}{dt^2} \log F(t) = -q(t)^2 $$

We then obtain the above formula by integrating twice. To justify this requires showing that $$\log F(t)$$ and $$\frac{d}{dt} \log F(t)$$ tend to zero at $$t = +\infty$$. Showing this requires an asymptotic analysis of the Fredholm determinant $$\det(1 - K)_{L^2(t,\infty)}$$ which is beyond the scope of this post.

The first step is to bring the $$t$$ dependence into the operator. Let

$$K_t(x,y) = K(x+t,y+t) = \int_t^\infty \mathrm{Ai}(x+s) \mathrm{Ai}(y+s) \, ds$$

Then $$F(t)=\det(1-K)_{L^2(t,\infty)} = \det(1-K_t)_{L^2(\mathbb{R}_+)}$$.

**Notation:** We let $$\tau_t$$ be the shift operator, so that $$(\tau_t \phi)(x) = \phi(x+t)$$ and $$D$$ be the derivative operator, $$(D\phi)(x) = \phi^\prime(x)$$. We shall be somewhat careless and not specify not what spaces these operators act on. Let us also denote the Airy function $$\mathrm{Ai} = A$$. $$\triangle$$

We see that $$\frac{d}{dt} K_t(x,y) = - A(x+t)A(y+t)$$. Thus

$$ \frac{d}{dt} K_t = - \tau_t A \otimes \tau_t A$$

**Remark:** I should signpost a point of rigour. We have calculated the derivative with respect to $$t$$ pointwise on the kernel, but in fact what we'd like is for the limit implicit in the derivative to exist in the *trace norm*, and a complete proof would show this. $$\triangle$$

Then we see by Jacobi's formula  

$$\frac{d}{dt} \log F(t) = -\mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} \frac{dK_t}{dt} \right) = \mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} \tau_t A \otimes \tau_t A \right)$$

**Remark:** Note that $$\mathrm{tr}_{L^2(\mathbb{R}_+)} \left(\psi \otimes \phi\right) = \langle \psi, \phi \rangle_{L^2(\mathbb{R}_+)} = \int_{\mathbb{R}_+} \psi(x) \phi(x) \, dx$$ (we will always take functions to be real valued). Note also that $$K_t$$, and hence $$(1-K_t)^{-1}$$, are symmetric operators with respect to this inner product. $$\triangle$$

Next we use the identity

$$\frac{d}{dt} (1-K_t)^{-1} = (1-K_t)^{-1} \frac{dK_t}{dt}(1-K_t)^{-1} = - (1-K_t)^{-1} (\tau_t A \otimes \tau_t A) (1-K_t)^{-1} $$

Observe that $$\mathrm{tr}((\alpha \otimes \beta)(\gamma \otimes \delta)) = \mathrm{tr}(\alpha \otimes \delta) \mathrm{tr}(\beta \otimes \gamma)$$.

$$\frac{d^2}{dt^2} \log F(t) = 2 \mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} D \tau_t A \otimes \tau_t A\right) - \left(\mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} \tau_t A \otimes \tau_t A\right)\right)^2$$

<h3>A hierarchy of coupled ODEs</h3>

Introduce the following notation,

$$ q_n(t) = ((1-K_t)^{-1} D^n \tau_t A)(0) $$

$$p_n(t) = \mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} D^n \tau_t A \otimes \tau_t A\right)$$

In this notation we have $$\frac{d^2}{dt^2} \log F(t) = 2 p_1(t) - p_0(t)^2$$.

We now compute

$$ \frac{d}{dt} q_n(t) = q_{n+1}(t)- q_0(t)p_n(t) $$

$$ \frac{d}{dt} p_n(t) = p_{n+1}(t)- p_0(t)p_n(t) + \underbrace{\mathrm{tr}_{L^2(\mathbb{R}_+)} \left( (1-K_t)^{-1} D^n \tau_t A \otimes  D\tau_t A\right)}_{(\ast)}$$

To compute $$(\ast)$$ we integrate by parts

$$(\ast)  = -q_n(t)(\tau_t A)(0) - \underbrace{\mathrm{tr}_{L^2(\mathbb{R}_+)} \left( D(1-K_t)^{-1} D^n \tau_t A \otimes  \tau_t A\right)}_{(\ast \ast)}$$

Next we use the identity $$[D,(1-K_t)^{-1}] = (1-K_t)^{-1} [D,K_t] (1-K_t)^{-1}$$. We now give the following important exercise.

**Exercise:** Let $$\phi \in L^2(\mathbb{R}_+)$$ be a sufficiently nice function (e.g. continuously differentiable and $$\phi^\prime \in L^2(\mathbb{R}_+)$$). Then

$$([D,K_t]\phi)(x) = - ((\tau_t A \otimes \tau_t A)\phi)(x) - \phi(0) K_t(x,0)$$

This yields the formula

$$(\ast \ast) = p_{n+1}(t)- p_0(t)p_n(t)-q_n(t)((1-K_t)^{-1} K_t \tau_t A)(0) $$

Using that $$(1-K_t)^{-1} K_t = (1-K_t)^{-1} - 1$$, we obtain a formula for $$\frac{d}{dt} p_n(t)$$. We thus obtain an infinite hierarchy of coupled ODEs, $$n \in \mathbb{N}$$,

$$ \boxed{\frac{d}{dt} q_n(t) = q_{n+1}(t)- q_0(t)p_n(t) }$$

$$ \boxed{\frac{d}{dt} p_n(t) = -q_{n}(t)q_0(t)}$$

**Exercise:** The quantity $$C = p_0(t)^2 - q_0(t)^2 - 2 p_1(t)$$ is conserved. (There are actually infinitely many such conserved quantities but we only need this one.)

**Corollary:** It seems reasonable that since the Airy function decreases rapidly at $$+\infty$$ that $$q_n$$ and $$p_n$$ should tend to zero at $$t\to +\infty$$. It therefore follows that $$C=0$$. From this it follows that

$$\frac{d^2}{dt^2} \log F(t) = - q_0(t)^2 $$

**Remark:** It is "obvious" that since $$K_t$$ is "small" for $$t\to +\infty$$

$$q_0(t) \approx (\tau_t A)(0) = \mathrm{Ai}(t)$$

This explains the boundary condition. $$\triangle$$

<h3>Closing up the system</h3>

Everything up until now has been "universal" -- in that we haven't used any properties of $$A$$ -- we have only used the Hankel composition structure of $$K$$. In particular, we haven't used that $$A$$ solves the *Airy equation*, $$D^2 A = M A$$, where $$M$$ is the operator such that $$(M \phi)(x) = x \phi(x)$$. Such a "non-universal" property allows us to close up the system and obtain an ODE for $$q_0$$. Note that $$(M \phi)(0) = 0$$.

From this we get

$$q_2(t) = t q_0(t) + ([(1-K_t)^{-1},M] \tau_t A)(0)$$

As before $$[(1-K_t)^{-1},M] = (1-K_t)^{-1} [ K_t, M](1-K_t)^{-1}$$. If we recall our two equivalent formulae for the Airy kernel we see that

$$[ K_t, M] = - \tau_t A \otimes D\tau_t A + D\tau_t A \otimes \tau_t A$$

This gives

$$q_2(t) = t q_0(t) + q_1(t)p_0(t)- q_0(t)p_1(t)$$

If we combine this formula with our relation $$p_0(t)^2 - q_0(t)^2 - 2 p_1(t) = 0$$ we find

$$\boxed{q_0^{\prime \prime}(t) = tq_0(t) + 2 q_0(t)^3}$$

which is Painlevé II.

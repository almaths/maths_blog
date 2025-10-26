---
layout: post
title:  "Asymptotic expansion of a partition function with a complex potential"
date:   2025-08-15 12:00:00 +0100
permalink: /complex-partition/
katex: True
---

The blog post will be an accessible introduction to my recent paper, joint with A. Guionnet and K. Kozlowski, [Asymptotic expansion of the partition function for β-ensembles with complex potentials](https://arxiv.org/abs/2411.10610). This post is based off presentations I have given of this paper. However because of time constraints many "tricks" from the paper are left out from the presentation. This post gives me the opportunity to explain such tricks.

<h2>Real integrals</h2>

In many problems in mathematics one encounters integrals of the form

$$I_N = \int_a^b \mathrm{e}^{N \varphi(x)}\, \mathrm{d}x$$

for some smooth (real-valued) function $$\varphi \in C^\infty([a,b])$$ and $$a<b$$ being finite real numbers. We are interested in the behaviour of $$I_N$$ for large $$N$$. By the triangle inequality we can make the trivial bound

$$\left| I_N \right| \leq (b-a) \, \mathrm{e}^{N \sup_{[a,b]}\varphi} \, .$$

Then, in particular, we have

$$\limsup_{N \to +\infty} \frac{\ln I_N}{N} \leq \sup_{[a,b]}\varphi \, .$$

In fact, the limit exists and we actually have equality here. To see this, observe that, by continuity of $$\varphi$$, for any $$\epsilon > 0$$ there must exist a subinterval $$J_\epsilon \subset I$$, of positive length $$ \lvert J_\epsilon \rvert  > 0$$, such that $$\sup_{[a,b]}\varphi \leq \varphi(x) + \epsilon$$ for all $$x \in J_\epsilon$$. Then we have 

$$I_N \geq \int_{J_\epsilon} \mathrm{e}^{N \varphi(x)}\, \mathrm{d}x \geq |J_\epsilon| \,  \mathrm{e}^{N \sup_{[a,b]}\varphi - N \epsilon} \, .$$

Then

$$\liminf_{N \to +\infty} \frac{\ln I_N}{N} \geq \sup_{[a,b]}\varphi - \epsilon \, .$$

Since $$\epsilon > 0$$ was arbitrary, we conclude that $$\frac{\ln I_N}{N}$$ converges and

$$\lim_{N \to +\infty} \frac{\ln I_N}{N} = \sup_{[a,b]}\varphi \, .$$

This analysis only gives the leading order of $$\frac{\ln I_N}{N}$$. To go further we need further assumptions. Let us suppose that the supremum is attained at a unique point in the interior of the interval, $$x^\ast \in (a,b)$$, and suppose that $$\varphi^{\prime \prime}(x^\ast) < 0$$. Then there exists a neighbourhood $$U$$ of $$x^\ast$$ and a smooth function $$\psi \in C^\infty(U)$$ such that

$$\varphi(x) = \varphi(x^\ast) - \psi(x)^2 \, .$$

Indeed, we can take this equation as *defining* $$\psi$$ and then find a sufficiently small neighbourhood of $$x^\ast$$ such that it is smooth and real valued. Furthermore $$\psi^\prime(x) > 0$$ for all $$x \in U$$. Then

$$I_N =  \left(  1 +\mathcal{O}(\mathrm{e}^{-NC}) \right) \int_{U} \mathrm{e}^{N \varphi(x)}\, \mathrm{d}x $$

for some $$C > 0$$ and then 

$$I_N = \left(  1 +\mathcal{O}(\mathrm{e}^{-NC}) \right) \mathrm{e}^{N \varphi(x^\ast)} \int_{\psi^{-1}(U)} \mathrm{e}^{-N u^2} (\psi^{-1})^\prime(u) \, \mathrm{d}u\, .$$

Taylor expanding $$(\psi^{-1})^\prime$$ at $$0$$ and integrating term by term (after taking the limits of integration to $$\pm \infty$$) we find an asymptotic series

$$I_N \sim \sqrt{\frac{2\pi}{-N \varphi^{\prime \prime}(x^\ast)}} \mathrm{e}^{N \varphi(x^\ast)} \left( 1+ \frac{A_1}{N} + \frac{A_2}{N^2} + \dots \right) \, .$$

Note that that we obtain an asymptotic series in $$\frac{1}{N}$$ because all the odd integrals vanish. This is the *Laplace method* in brief.

**Remark:** Note that we have a kind of "central limit theorem" happening inside the integral, where if we think of the integrand as representing a distribution function, then $$\sqrt{N}(x - x^\ast)$$ is asymptotically Gaussian with mean $$0$$ and variance $$\frac{1}{\sqrt{-  \varphi^{\prime \prime}(x^\ast)}}$$. 

<h2>Contour integrals</h2>

Let us now repeat the analysis, where now $$a, b \in \mathbb{C}$$ and $$\varphi$$ is analytic on some domain containing $$a$$ and $$b$$. Then consider

$$I_N = \int_a^b \mathrm{e}^{N \varphi(z)} \, \mathrm{d} z \, .$$

If we choose a path $$\Gamma[a,b]$$ between the endpoints, we have

$$|I_N| \leq |\Gamma[a,b]| \, \mathrm{e}^{N \sup_{\Gamma[a,b]}  \Re \varphi} $$

where $$\lvert \Gamma[a,b]\rvert $$ is the arc-length of the contour. Hence 

$$ \limsup_{N \to +\infty} \frac{\ln|I_N|}{N} \leq  \sup_{z \in \Gamma[a,b]}  \Re \varphi(z)  \, .$$

However by Cauchy's theorem, the original integral is unchanged under deformations of the contour. Hence we may optimise our bound over some appropriate class of homotopic contours $$\mathcal{T}$$.

$$ \limsup_{N \to +\infty} \frac{\ln|I_N|}{N} \leq \inf_{\Gamma[a,b] \in \mathcal{T}} \sup_{z \in \Gamma[a,b]}  \Re \varphi(z)  \, .$$

This is the essence of the *saddle point method*. Suppose we found a curve $$\Gamma[a,b]$$ and a point $$z^\ast \in \Gamma[a,b]$$ such that the above inf-sup is actually attained---what would it look like? Well, *along* the curve $$\Re \varphi(z)$$ would attain a *maximum* at $$z = z^\ast$$. However *perpendicular* to the curve $$\Re \varphi(z)$$ must attain a *minimum*, since otherwise a lateral deformation of the contour could reduce the value of $$\sup_{z \in \Gamma[a,b]}  \Re \varphi(z)$$. Hence $$z= z^\ast$$ must be a saddle point of $$\Re \varphi(z)$$.

Next, we note that since $$\varphi$$ is analytic, $$\Re \varphi$$ is a harmonic function, and so has no local maxima or minima. (Alternatively one could apply the maximum modulus principle to $$\mathrm{e}^{\varphi}$$.) That is, the only stationary points of $$\Re \varphi$$ (points where the two-dimensional gradient vanishes) are saddle points. And by the Cauchy-Riemann equations the gradient of $$\Re \varphi$$ vanishes if and only if the holomorphic derivative vanishes, $$\varphi^\prime(z) = 0$$. Hence the inf-sup problem becomes one of finding these stationary points and appropriately traversing them. Indeed, we could choose the contour such that it traverses $$z^\ast$$ such that $$\Im \varphi$$ is constant in a neighbourhood of $$z^\ast$$ along the contour, and this would essentially reduce the problem to a real integral discussed previously. Such a contour is called a *steepest descent contour* since it is the direction of steepest descent of $$\Re \varphi$$. In this way one can show that $$\frac{\ln\lvert I_N\rvert }{N}$$ really converges to $$\inf_{\Gamma[a,b] \in \mathcal{T}} \sup_{z \in \Gamma[a,b]}  \Re \varphi(z)$$ as $$N \to +\infty$$. (This way of deriving the saddle point method by trying to optimise the triangle inequality I learnt from de Bruijn's book *Asymptotic Methods in Analysis*.)

<h2>Partition functions of β-ensembles</h2>

A β-ensemble is a random collection of $$N$$ particles on the real line, $$x_1, \dots, x_N \in \mathbb{R}$$, with joint density (with respect to Lebesgue measure on $$\mathbb{R}^N$$)

$$\varrho_N(x_1, \dots, x_N) = \frac{1}{\mathsf{Z}_N[V]} \prod_{1 \leq i < j \leq N}|x_i - x_j|^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta V(x_k)}\, .$$

$$\mathsf{Z}_N[V]$$ is a normalisation constant know as the *partition function*,

$$\mathsf{Z}_N[V] \overset{\mathrm{def}}{=} \int_{\mathbb{R}^N} \prod_{1 \leq i < j \leq N}|x_i - x_j|^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta V(x_k)}\, \mathrm{d}x . $$

$$\beta$$ is a fixed positive parameter, which can be interpreteted as the inverse of the "temperature" $$\beta = \frac{1}{T} > 0$$. This is because $$\varrho_N$$ can be interpreted as the thermal distribution of $$N$$ particles with Hamiltonian

$$\mathcal{H}(x_1, \dots, x_N) = \sum_{1 \leq i < j \leq N} \ln \frac{1}{|x_i - x_j|} + N \sum_{k=1}^N V(x_k)\, .$$

$$V$$ is a confining potential which grows sufficiently fast such that the density is normalisable. For "special" potentials $$V$$ there is a matrix-model representation of β-ensembles ([Dumitriu-Edelman, 2002](https://arxiv.org/abs/math-ph/0206043)) and also for more general $$V$$ if $$\beta = 1,2,4$$. The asymptotic expansion of $$\ln \mathsf{Z}_N[V]$$, also known as the free energy, is a central question of interest in statistical mechanics. A full asymptotic expansion was achieved, in a certain off-critical regime, by Borot and Guionnet (see [here](https://arxiv.org/abs/1107.1167) and [here](https://arxiv.org/abs/1303.1045)). The integrand of $$\mathsf{Z}_N[V]$$ depends on $$N$$ but so does the number of integrations, so the Laplace method breaks down. Instead Borot and Guionnet use logarithmic potential theory and the method of Dyson-Schwinger equations to obtain their asymptotic expansion. Thus these methods could be regarded as an ∞-dimensional version of the Laplace method.

Let me outline how the logarithmic potential theory arises. Define the empirical measure

$$L_N^{(\mathbf{x})} = \frac{1}{N} \sum_{k=1}^N \delta_{x_k} \, .$$

Then 

$$\mathsf{Z}_N[V] = \int_{\mathbb{R}^N} \exp\left( - \frac{\beta N^2}{2}\int_{\substack{\mathbb{R}^2 \\ x \neq y}} \ln \frac{1}{|x-y|} \, \mathrm{d}L_N^{(\mathbf{x})}(x) \otimes \mathrm{d}L_N^{(\mathbf{x})}(y) - \beta N^2 \int_{\mathbb{R}} V(x) \, \mathrm{d}L_N^{(\mathbf{x})}(x) \right)\, \mathrm{d}x . $$

We remark that $$\{ L_N^{(\mathbf{x})} \}_{\substack{\mathbf{x \in \mathbb{R}^N} \\ N \geq 1}}$$ forms a dense subset of $$\mathcal{M}_1(\mathbb{R})$$, the space of Borel probability measures on $$\mathbb{R}$$. Thus if we ignore the diagonal we expect that 

$$  \limsup_{N \to +\infty}\frac{\ln \mathsf{Z}_N[V]}{N^2} \leq - \frac{\beta}{2} \inf_{\mu \in \mathcal{M}_1(\mathbb{R})} \mathtt{I}_V[\mu] $$

where 

$$\mathtt{I}_V[\mu] = \int_{\mathbb{R}^2} \left( \ln \frac{1}{|x-y|} + V(x) + V(y) \right) \, \mathrm{d}\mu(x) \otimes \mathrm{d}\mu(y) \, .$$

In fact $$\frac{\ln \mathsf{Z}_N[V]}{N^2} \to - \frac{\beta}{2} \inf_{\mu \in \mathcal{M}_1(\mathbb{R})} \mathtt{I}_V[\mu]$$.

<h2>A complex partition function</h2>

We now turn to the problem with which my paper with Guionnet and Kozlowski was concerned. This involves studying a "complexification" of the β-ensemble partition function $$\mathsf{Z}_N[V]$$. The first way we could complexify is to make the potential $$V$$ complex. However, since we want to develop an analogue of the saddle point/steepest descent method, we want the integrand to be analytic. Thus let us take $$V$$ to be a polynomial with complex coefficients, say of degree $$\kappa \geq 2$$. By rescaling we can choose the leading coefficient to be any positive number, e.g. $$V(z) = \frac{1}{\kappa} z^{\kappa} + \mathcal{O}(z^{\kappa-1})$$.

Let us then take

$$\mathcal{Z}_{N , \Gamma}[V] \overset{\mathrm{def}}{=} \int_{\Gamma^N} \prod_{1 \leq i < j \leq N}(z_i - z_j)^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta V(z_k)} \, \mathrm{d}\mathbf{z} \, .$$

So that the integrand is analytic we should take $$\beta$$ to be a positive integer. In fact, it should be an even integer because if $$\beta$$ was odd then the interchange of any two integration variables would produce a minus sign, and by symmetry of the integration region $$\mathcal{Z}_{N , \Gamma}[V] = 0$$. Thus $$\beta$$ must be even for the problem to be nontrivial.

Finally, what conditions should we put on $$\Gamma$$?

**Definition:** Fix $$\alpha, \alpha^\prime \in [\![ 1, \kappa -1 ]\!]$$ to be two distinct integers. We say a contour $$\Gamma$$ is "admissible" (relative to $$(\alpha, \alpha^\prime)$$) if

(1) It consists of a finite number of $$C^1$$ arcs.

(2) It is connected.

(3) There exists an $$R > 0$$ sufficiently large so that

$$\Gamma \setminus D_R(0) = \mathrm{e}^{\frac{2\pi i \alpha}{\kappa}} [R,+\infty) \cup \mathrm{e}^{\frac{2\pi i \alpha^\prime}{\kappa}} [R,+\infty)$$

where $$D_R(0)$$ is the open disk of radius $$R$$ centred at $$0$$. We require incoming orientation on $$\mathrm{e}^{\frac{2\pi i \alpha}{\kappa}} [R,+\infty)$$ and outgoing orientation on $$\mathrm{e}^{\frac{2\pi i \alpha^\prime}{\kappa}} [R,+\infty)$$. $$\triangle$$

This means that outside a large compact set $$\Gamma$$ consists of two parts: an incoming ray and an outgoing ray, and these rays should lie along the line proportional to a $$\kappa$$th root of unity. A curve that satisfies all of the above is said to be "admissible." These properties mean that $$\lvert \mathrm{e}^{- V(z)}\rvert \to 0$$ rapidly along these rays. We could, of course, be less restrictive and allow contours that run "close" to the rays $$\mathrm{e}^{\frac{2\pi i \alpha}{\kappa}} \mathbb{R}_+ \cup \mathrm{e}^{\frac{2\pi i \alpha^\prime}{\kappa}} \mathbb{R}_+$$, however what happens outside of a sufficiently large compact set will make no contribution to the asymptotic series, and we can also deform our contour to lie exactly upon these rays. Thus we should really think of $$\mathcal{Z}_{N , \Gamma}[V]$$ as being a function of the homotopy class of $$\Gamma$$, which is labelled by $$(\alpha, \alpha^\prime)$$. Note that if we took $$\alpha = \alpha^\prime$$ then $$\mathcal{Z}_{N , \Gamma}[V] = 0$$, hence we exclude this trivial case.

**Remark:** Note that interchanging $$\alpha$$ and $$\alpha^\prime$$ changes $$\mathcal{Z}_{N , \Gamma}[V]$$ by a factor of $$(-1)^N$$. $$\triangle$$

The central question of the paper is how $$\ln \mathcal{Z}_{N , \Gamma}[V]$$ behaves as $$N \to +\infty$$. Following a similar reasoning to case of the partition function of a real β-ensemble, we have

$$  \limsup_{N \to +\infty}\frac{\ln |\mathcal{Z}_{N , \Gamma}[V]|}{N^2} \leq - \frac{\beta}{2} \inf_{\mu \in \mathcal{M}_1(\Gamma)} \mathtt{I}_V[\mu] $$

where 

$$\mathtt{I}_V[\mu] \overset{\mathrm{def}}{=} \int_{\mathbb{C}^2} \left( \ln \frac{1}{|z-w|} + \varphi(z) + \varphi(w) \right) \, \mathrm{d}\mu(z) \otimes \mathrm{d}\mu(w)$$

for $$\varphi = \Re V$$, and where $$\mathcal{M}_1(\Gamma)$$ is the space of Borel probability measures on $$\Gamma$$. Then, optimising, we find

$$  \limsup_{N \to +\infty}\frac{\ln |\mathcal{Z}_{N , \Gamma}[V]|}{N^2} \leq - \frac{\beta}{2} \sup_{\tilde{\Gamma} \in \mathcal{T}} \inf_{\mu \in \mathcal{M}_1(\tilde{\Gamma})} \mathtt{I}_V[\mu] \,  $$

where $$\mathcal{T}$$ is the space of admissible contours (we suppress the $$\alpha, \alpha^\prime$$ dependence since it is fixed throughout). If we believe this bound to be optimal then we have a kind of ∞-dimensional saddle point.

**Fact:** By potential theoretic arguments one can show that if $$\Gamma$$ is admissible then $$\mathtt{I}_V$$ has unique minimiser $$\mu^\Gamma \in \mathcal{M}_1(\Gamma)$$, which we call the *equilibrium measure* associated to $$\Gamma$$.

**Definition:** $$\Gamma_{\mathrm{eq}}$$ is said to solve the "max-min energy problem" if 

$$\sup_{\tilde{\Gamma} \in \mathcal{T}} \inf_{\mu \in \mathcal{M}_1(\tilde{\Gamma})} \mathtt{I}_V[\mu] = \mathtt{I}_V[\mu^{\Gamma_{\mathrm{eq}}}] \, .$$

**Theorem ([Silva-Kuijlaars, 2015](https://arxiv.org/abs/1311.7026)):** Let $$V$$ be a polynomial (with complex coefficients). Then there exists an admissible $$\Gamma_{\mathrm{eq}}$$ solving the max-min energy problem. $$\Gamma_{\mathrm{eq}}$$ is not unique however $$\mu^{\Gamma_{\mathrm{eq}}}$$ is unique (so that all solutions have the same equilibrium measure). $$\triangle$$

**Definition:** (1) Define the logarithmic potential associated to a probability measure $$\mu$$ as

$$U[\mu](z) \overset{\mathrm{def}}{=} \int_{\mathbb{C}} \ln \frac{1}{|z-w|} \, \mathrm{d}\mu(w)$$

wherever this makes sense. Of course, if $$z \not\in \mathrm{supp}\, \mu$$ and $$\mu$$ is compactly supported then the integral converges.

(2) Similarly, define the Cauchy transform of probability $$\mu$$

$$C[\mu](z) \overset{\mathrm{def}}{=} \frac{1}{2\pi i} \int_{\mathbb{C}} \frac{1}{w-z} \, \mathrm{d}\mu(w)$$

wherever this makes sense.

(3) We say that $$\Gamma$$ is an $$S$$-curve in the external field $$\varphi$$ if there is a set of zero capacity $$E$$ such that for every $$z \in \mathrm{supp}\, \mu^\Gamma \setminus E$$ there is a neighbourhood $$U$$ of $$z$$ such that $$\mathrm{supp}\, \mu^\Gamma \cap U$$ is an analytic arc, and 

$$\frac{\partial}{\partial n^+} \left( U[\mu^\Gamma]+ \varphi \right)(z) = \frac{\partial}{\partial n^-} \left( U[\mu^\Gamma]+ \varphi \right)(z)$$

where $$\frac{\partial}{\partial n^+}$$ and $$\frac{\partial}{\partial n^-}$$ are the directional derivatives normal to the contour (away from the contour). $$\triangle$$

**Theorem ([Silva-Kuijlaars, 2015](https://arxiv.org/abs/1311.7026)):** Let $$V$$ be a polynomial (with complex coefficients) and $$\Gamma_{\mathrm{eq}}$$ be an associated solution to the max-min energy problem. Then we have the following.

1) $$\mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}}$$ is a finite collection of (bounded) analytic arcs.

2) $$\mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}}$$ is an $$S$$-curve in the external field $$\varphi = \Re V$$.

3) There exists a polynomial $$R$$ such that

$$R(z) = ( V^\prime(z) + 2 \pi i \, C[\mu^{\Gamma_{\mathrm{eq}}}](z))^2 \, .$$

4) $$\mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}}$$ consists of critical trajectories of the quadratic differential $$ - R(z)\, \mathrm{d}z^2 $$. $$\triangle$$

**Remark:** 1) Taking the square root we find $$\sqrt{R(z)} = V^\prime(z) + 2\pi i \, C[\mu^{\Gamma_{\mathrm{eq}}}](z)$$. The left hand side is analytic everwhere except on its branch cuts. The right hand side is analytic everywhere except on the support of $$\mu^{\Gamma_{\mathrm{eq}}}$$. Hence the branch cuts of $$\sqrt{R(z)}$$ are the support of the equilibrium measure.

2) By Plemelj's formula we have $$\frac{1}{i \pi } \sqrt{R(z)} \, \mathrm{d}z = \mathrm{d}\mu^{\Gamma_{\mathrm{eq}}}(z) > 0$$ whenever the density is positive. Squaring both sides (which amounts to ignoring the orientation of the contour) we find $$\mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}}$$ is a critical trajectory of $$ - R(z)\, \mathrm{d}z^2 $$.

3) From the fact that $$\sqrt{R(z)}$$ changes sign across the branch cut we have $$\sqrt{R(z)}_+ + \sqrt{R(z)}_- = 0$$. This gives

$$V^\prime(z) + \mathrm{p.v.}\int \frac{1}{w-z} \, \mathrm{d}\mu^{\Gamma_{\mathrm{eq}}}(w) = 0, \quad \quad \quad z \in \mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}} \, .$$

This equation can be regarded as a kind of "complexified" Euler-Lagrange equation, since taking its perpendicular to the contour one obtains the Euler-Lagrange and $$S$$-curve conditions. This again supports the interpretation of an ∞-dimensional saddle point, since it is analogous to the saddle point equation $$\varphi^\prime(z) = 0$$, which can also be broken down into two components representing the maximum along the curve and the minimum perpendicular to the curve.

The Euler-Lagrange condition for the minimisation of the energy on the curve states that there is a constant $$C_\varphi$$ such that $$\varphi(z) + U[\mu^{\Gamma_{\mathrm{eq}}}](z) \geq C_\varphi$$ throughout the curve, with equality when $$z \in \mathrm{supp} \, \mu^{\Gamma_{\mathrm{eq}}}$$. Hence let us define the *effective potential* as

$$\varphi_{\mathrm{eff}}(z)= \varphi(z) + U[\mu^{\Gamma_{\mathrm{eq}}}](z) - C_\varphi \, .$$

**Definition:** We say that the potential $$V$$ is *regular* if it is possible to choose $$\Gamma_{\mathrm{eq}}$$ such that

1) The polynomial $$R$$ vanishes nowhere on $$\Gamma_{\mathrm{eq}}$$ except at the endpoints of $$\mathrm{supp}\, \mu_{\mathrm{eq}}$$, and at these endpoints $$R$$ has only simple zeros. (Note this is equivalent to saying that the density of the equilibrium measure as square-root type vanishing.)

2) $$\varphi_{\mathrm{eff}}(z) > 0$$ for all $$z \in \Gamma_{\mathrm{eq}} \setminus \mathrm{supp}\, \mu_{\mathrm{eq}}$$.

$$V$$ is *one-cut regular* if it is regular and $$\mathrm{supp}\, \mu_{\mathrm{eq}}$$ is connected. $$\triangle$$

It is believed that the regular potentials are "generic", in the sense of being an open set of full measure within the space of all potentials. However for this max-min energy problem this is unproven (the claim that the regular potentials form an open set is a recent result of [Bertola, Bleher, Gharakhloo, McLaughlin and Tovbis, 2022](https://link.springer.com/chapter/10.1007/978-3-031-13851-5_10)).

**Definition:** Define the "$$g$$-function" as

$$g[\mu^{\Gamma_{\mathrm{eq}}}](z) \,\, \text{"="} \,\,   \int_{\Gamma_{\mathrm{eq}}} \ln(z-w) \, \mathrm{d}\mu_{\mathrm{eq}}(w) \, . $$

We put quotation marks because some care must be taken with the branch cut. To be more precise, we want to choose the branch cut of $$g[\mu^{\Gamma_{\mathrm{eq}}}]$$ along $$\Gamma_{\mathrm{eq}}$$ and require that

$$\frac{\mathrm{d}}{\mathrm{d}z} g[\mu^{\Gamma_{\mathrm{eq}}}](z) =   \int_{\Gamma_{\mathrm{eq}}} \frac{1}{z-w} \, \mathrm{d}\mu_{\mathrm{eq}}(w) \, . $$

Next, define the "complex energy" as

$$\mathcal{I}_{\Gamma_{\mathrm{eq}}}[\mu^\Gamma_{\mathrm{eq}}] \overset{\mathrm{def}}{=} \int_{\Gamma_{\mathrm{eq}}}\Big( \frac{1}{2} g^+[\mu^{\Gamma_{\mathrm{eq}}}] + \frac{1}{2} g^-[\mu^{\Gamma_{\mathrm{eq}}}] + 2 V\Big) \, \mathrm{d} \mu^{\Gamma_{\mathrm{eq}}} $$

where $$g^\pm[\mu^{\Gamma_{\mathrm{eq}}}]$$ are the left and right boundary values of the $$g$$-function up to the curve. $$\triangle$$

**Remark:** Note that the real part of the complex energy functional is the real energy functional, which justifies the name, $$\Re \mathcal{I}_{\Gamma_{\mathrm{eq}}}[\mu^\Gamma_{\mathrm{eq}}] = \mathtt{I}_V[\mu^\Gamma_{\mathrm{eq}}] \,$$ $$\triangle$$.

We now have enough tools to state the main theorem

**Theorem:** (Guionnet, Kozlowski, L., 2025) Let $$V(z) = \frac{z^\kappa}{\kappa} + \mathcal{O}(z^{\kappa-1})$$ be a polynomial, and one-cut regular as a potential. Then there exists a sequence of complex numbers $$(F_k(\beta,V))_{k=-2}^{\infty}$$ such that

$$\boxed{\ln \mathcal{Z}_{N , \Gamma}[V] \sim  \frac{\beta}{2} N \ln N + \frac{3 + \frac{\beta}{2}+ \frac{2}{\beta}}{12} \ln N + \sum_{k=-2}^{+\infty} \frac{1}{N^k } F_k(\beta,V) }$$

where this is an asymptotic series at scale $$\frac{1}{N}$$ in the sense of Poincaré. $$F_{-2}$$ and $$F_{-1}$$ have explicit formulae in terms of the equilibrium measure.

$$\begin{align*}
F_{-2}(\beta,V) &= - \frac{\beta}{2} \mathcal{I}_{\Gamma_{\mathrm{eq}}}[\mu^\Gamma_{\mathrm{eq}}] \\
F_{-1}(\beta,V) &= \Big( \frac{\beta}{2}-1 \Big) \Big[ \int_{\Gamma_{\mathrm{eq}}} \ln \frac{\mathrm{d}\mu^\Gamma_{\mathrm{eq}}}{\mathrm{d}z} \, \mathrm{d}\mu^\Gamma_{\mathrm{eq}} + \ln \frac{\beta}{2}  \Big]  + \frac{\beta}{2} \ln \frac{2\pi}{\mathrm{e}}\, . \,  \triangle
\end{align*}$$

Thus, the leading term at scale $$N^2$$ has the form of a complex *energy*, whereas the term at scale $$N$$ involves a complex *entropy* (note that $$\frac{\mathrm{d}\mu^\Gamma_{\mathrm{eq}}}{\mathrm{d}z}$$ is not real valued).

Let us sketch the proof of this result. The proof involves a combination of pre-existing theory and novel tricks. The pre-existing theory consists primarily in the method developed by Guionnet and Borot to analyse the partition function of β-ensembles on the real line. However the fact that the integrand is complex-valued introduces new obstacles which must be overcome.

<h2>A differential identity</h2>

Suppose we have a family of potentials $$\{ V_t \}_{t \in [0,1]}$$ smoothly parametrised by $$t$$, such that the potential we are interested in $$V = V_1$$. Then a short calculation shows that

$$\frac{\partial}{\partial t} \ln \mathcal{Z}_{N , \Gamma}[V_t] = - \beta N^2 \left\langle L_N^{(\mathbf{z})}\Big( \frac{\partial V_t}{\partial t} \Big) \right\rangle_{N, V_t}  $$

where

$$\left\langle \, \cdot \,  \right\rangle_{N, V} = \frac{1}{\mathcal{Z}_{N , \Gamma}[V]} \int_{\Gamma^N} \left( \, \cdot \,  \right) \,  \prod_{1 \leq i < j \leq N}(z_i - z_j)^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta V(z_k)} \, \mathrm{d}\mathbf{z}$$

is the expectation with respect to the "complex measure." If we then integrate up we have

$$\ln \mathcal{Z}_{N , \Gamma}[V] = \ln \mathcal{Z}_{N , \Gamma}[V_0] - \beta N^2 \int_0^1 \left\langle L_N^{(\mathbf{z})}\Big( \frac{\partial V_t}{\partial t} \Big) \right\rangle_{N, V_t} \, \mathrm{d}t \, . $$

If we take $$V_0$$ to be quadratic then $$\ln \mathcal{Z}_{N , \Gamma}[V_0]$$ can be computed by a Selberg integral and its asymptotics is then known via the asymptotics of the Barnes G-function. Thus our problem reduces to ab asymptotic expansion of $$\left\langle L_N^{(\mathbf{z})}\Big( \frac{\partial V_t}{\partial t} \Big) \right\rangle_{N, V_t}$$, ensuring that our error terms are sufficiently uniform in $$t \in [0,1]$$ that we can integrate them. $$L_N^{(\mathbf{z})}(f) = \frac{1}{N} \sum_{j=1}^N f(z_j)$$ is a linear statistic, hence our problem reduces to the asymptotics of moments (in this case the first moment) of a linear statistic.

<h2>Dyson-Schwinger equations</h2>

The method we use to study the asymptotics of moments of linear statistics is the method of *Dyson-Schwinger equations*. Define the fluctuation measure as

$$\mathrm{Fluct}_N = N \left( L_N^{(\mathbf{z})} - \mu^{\Gamma_{\mathrm{eq}}} \right) \, .$$

and let $$f_0, \dots, f_k$$ be a collection of analytic test functions, growing sufficiently slowly on $$\Gamma$$ such that the associated integrals converge. Then the DS equations read as follows.

$$\begin{align*}
\left\langle \mathrm{Fluct}_N( \Xi f_0) \prod_{p=1}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V}
&= \left( \frac{1}{\beta} - \frac{1}{2}\right) \mu^{\Gamma_{\mathrm{eq}}} (f_0^\prime) \left\langle \prod_{p=1}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V} \\
&+ \frac{1}{N}\left( \frac{1}{\beta} - \frac{1}{2}\right) \left\langle \mathrm{Fluct}_N(f_0^\prime) \prod_{p=1}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V}\\
&+ \frac{1}{\beta} \sum_{q=1}^k \mu^{\Gamma_{\mathrm{eq}}} (f_0 f_q^\prime) \left\langle \prod_{\substack{p=1\\p \neq q}}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V}  \\
&+ \frac{1}{\beta N} \sum_{q=1}^k  \left\langle \mathrm{Fluct}_N(f_0 f_q^\prime) \prod_{\substack{p=1\\p \neq q}}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V} \\
&+ \frac{1}{2 N} \left\langle \mathrm{Fluct}_N^{\otimes 2}\Big(\frac{f_0(z) - f_0(w)}{z-w}\Big) \prod_{\substack{p=1\\p \neq q}}^k \mathrm{Fluct}_N(f_p)  \right\rangle_{N,V}
\end{align*} $$

where $$\Xi$$ is the "master operator" which acts by

$$\Xi f (z) \overset{\mathrm{def}}{=} V^\prime(z) f(z) + \int_{\Gamma_{\mathrm{eq}}} \frac{f(z) - f(w)}{z-w} \, \mathrm{d}\mu^{\Gamma_{\mathrm{eq}}}(w)\, .$$

These equations can be derived by integration by parts. To make these equations useful we must do two things. Firstly, we must invert the master operator $$\Xi$$. When we are in the one-cut regular regime this is possible (up to constants, but a constant gives zero under the fluctuation measure). Hence we make this hypothesis. Replacing $$f_0 \to \Xi^{-1}(f_0)$$ throughout, our DS equations relate moments of $$\mathrm{Fluct}_N$$ at different orders. As it stands, these equations are not closed, however they are *asymptotically* closed. By this I mean that $$\mathrm{Fluct}_N(f) = \mathcal{O}(1)$$; more precisely, there are constants $$C_k$$ such that

$$ \left|\left\langle \left| \mathrm{Fluct}_N(f) \right|^k \right\rangle_{N,V} \right| \leq C_k \, .$$

This motivates the original factor of $$N$$ in the definition of $$\mathrm{Fluct}_N$$. Substituting this bound into the Dyson-Schwinger equations in we can recursively derive a $$\frac{1}{N}$$ expansion. For example

$$\left\langle L_N^{(\mathbf{z})}(f) \right\rangle_{N,V} = \mu^{\Gamma_{\mathrm{eq}}}(f) + \frac{1}{N} \left( \frac{1}{\beta} - \frac{1}{2} \right) \mu^{\Gamma_{\mathrm{eq}}} (\Xi^{-1}(f)^\prime) + \mathcal{O}(N^{-2}) \, . $$

**Remark:** Note that in the case of $$\beta = 2$$, it becomes a $$\frac{1}{N^2}$$ expansion.

**Remark:** Note that if we discard the subleading terms in the DS equations, we essentially have the formula from "Wick's theorem" which recursively relates Gaussian moments. This means that the linear statistics obey a kind of complex CLT, such that their cumulants of degree 3 and higher all tend to zero.

I am simplifying somewhat since in the paper it is useful to restrict to a compact subset of the contour, but doing so introduces boundary terms which must be dealt with. Moving forward, we immediately can obtain the $$\frac{1}{N}$$ expansion of the partition function. However I have not explained how to obtain the bound $$ \left\lvert\left\langle \left\lvert \mathrm{Fluct}_N(f) \right\rvert^k \right\rangle_{N,V} \right\rvert \leq C_k$$ and indeed this is the chief difficulty of the paper. The challenge to bounding $$\left\lvert\left\langle \cdot \right\rangle_{N,V} \right\rvert$$ is that it involves dividing by the complex partition function, and this is a complex (and so oscillatory) integral, so in principle it could be very small.

Let $$\gamma : \mathbb{R} \longrightarrow \mathbb{C}$$ be a parametrisation of the contour $$\Gamma_{\mathrm{eq}}$$. Define the "real model" as 

$$\begin{align*}
\mathbb{E}_{N,V}[\cdot] &= \frac{1}{\mathsf{Z}_{N,\gamma}[V]} \int_{\mathbb{R}^N} (\cdot) \prod_{1 \leq i < j \leq N}\lvert \gamma(x_i) - \gamma(x_j)\rvert^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta \tilde{\varphi}(x_j)} \, \mathrm{d}x  \\
\mathsf{Z}_{N,\gamma}[V] &= \int_{\mathbb{R}^N} \prod_{1 \leq i < j \leq N}\lvert \gamma(x_i) - \gamma(x_j)\rvert^\beta \prod_{k=1}^N \mathrm{e}^{-N \beta \tilde{\varphi}(x_k)} \, \mathrm{d}x
\end{align*}$$

where $$\tilde{\varphi} = \Re V \circ \gamma$$. Hence by the triangle inequality, for any non-negative function $$F(z_1, \dots, z_N) \geq 0$$ we have

$$\left\lvert\left\langle F \right\rangle_{N,V} \right\rvert \leq \frac{\mathsf{Z}_{N,\gamma}[V]}{|\mathcal{Z}_{N,\Gamma}[V]|}\mathbb{E}_{N,V}[F\circ \gamma] \, .$$

This "real model" $$\mathbb{E}_{N,V}$$ has the advantage that it can be studied by existing techniques. In particular 


$$\begin{align*}
\prod_{1 \leq i < j \leq N}\lvert \gamma(x_i) - \gamma(x_j)\rvert^\beta &= \prod_{1 \leq i < j \leq N}|x_i - x_j|^\beta \prod_{i,j=1}^N u(x_i, x_j) \prod_{k=1}^N v(x_k) \\
u(x, y) &= \left| \frac{\gamma(x) - \gamma(y)}{x-y}\right|^{\frac{\beta}{2}} \\
v(x) &= |\gamma^\prime (x)|^{-\frac{\beta}{2}}
\end{align*}$$

Thus this reduces to the usual β-ensemble but with some smooth 2-body and 1-body interactions, which has been studied by [Borot, Guionnet and Kozlowski, 2015](https://arxiv.org/abs/1312.6664). This is done by the method of Dyson-Schwinger equations, where now the "master operator" is more complicated and so more involved to invert (we remark that this is not an "integrable" model, which shows that DS equations are quite robust to the form of the interaction). The key conclusion is that one can compute asymptotically the moments of linear statistics like

$$\sum_{k=1}^N f(x_k) $$

and show that they are bounded in $$N$$. These techniques also allow one to prove a CLT for such linear statistics under the real model. This means that our problem reduces to bounding the ratio of the partition functions.

<h2>The ratio of the partition functions</h2>

In this final section we sketch how to bound $$\frac{ \lvert\mathcal{Z}_{N,\Gamma}[V]\rvert}{\mathsf{Z}_{N,\gamma}[V]}$$ from below. The idea is to separate the integrand of $$\mathcal{Z}_{N,\Gamma}[V]$$ into a modulus and a phase, and then treat the integral as an expectation of the oscillatory part with respect to the modulus. Given $$z \in \mathbb{C}$$ define the phase $$\mathrm{ph}(z) = \frac{z}{\lvert z \rvert}$$. Notice that $$\mathrm{ph}(rz) = \mathrm{ph}(z)$$ for any $$r > 0$$ and $$\mathrm{ph}(zw) = \mathrm{ph}(z) \mathrm{ph}(w)$$. Hence, using that $$\prod_{i < j }(x_i - x_j)^\beta > 0$$ for $$\beta > 0$$ even,

$$\begin{align*}&\mathrm{ph}\left[ \prod_{1 \leq i < j \leq N}(\gamma(x_i) -\gamma(x_j) )^\beta \prod_{k=1}^N \mathrm{e}^{- \beta N V(\gamma(x_k))} \gamma^\prime(x_k)\right] \\
&= \mathrm{ph}\left[ \prod_{1 \leq i < j \leq N}\left(\frac{\gamma(x_i) -\gamma(x_j)}{x_i - x_j} \right)^\beta \prod_{k=1}^N \mathrm{e}^{- \beta N V(\gamma(x_k))} \gamma^\prime(x_k)\right]\\
&= \mathrm{ph}\left[ \prod_{i,j=1}^N \left(\frac{\gamma(x_i) -\gamma(x_j)}{x_i - x_j} \right)^\frac{\beta}{2} \prod_{k=1}^{N} \mathrm{e}^{- \beta N  V(\gamma(x_k))} \gamma^{\prime}(x_k)^{1-\frac{\beta}{2}}\right] \\
&= \exp\left(\frac{i \beta N^2}{2}\int_{\mathbb{R}^2} a \, \mathrm{d}L_N^{(\mathbf{x})}\otimes \mathrm{d}L_N^{(\mathbf{x})} - i \beta N^2 \int_{\mathbb{R}} \Im V \, \mathrm{d}L_N^{(\mathbf{x})} + i N \left( 1-  \frac{\beta}{2}\right) \int_{\mathbb{R}}\arg \, \gamma^\prime \, \mathrm{d}L_N^{(\mathbf{x})} \right)
\end{align*}$$

where $$a(x,y) = \arg \, \frac{\gamma(x)-\gamma(y)}{x-y}$$. We observe that 

$$ \frac{ \mathcal{Z}_{N,\Gamma}[V]}{\mathsf{Z}_{N,\gamma}[V]}  = \mathbb{E}_{N,V}\left[ \mathrm{e}^{\frac{i \beta N^2}{2}\int_{\mathbb{R}^2} a \, \mathrm{d}L_N^{(\mathbf{x})}\otimes \mathrm{d}L_N^{(\mathbf{x})} - i \beta N^2 \int_{\mathbb{R}} \Im V \, \mathrm{d}L_N^{(\mathbf{x})} + i N \left( 1-  \frac{\beta}{2}\right) \int_{\mathbb{R}}\arg \, \gamma^\prime \, \mathrm{d}L_N^{(\mathbf{x})} }\right] \, .$$

Now, let us recall that the equilibrium measure for the real model is just the pullback of the equilibrium measure $$\mu_{\mathrm{eq}}$$ under $$\gamma$$. Let us call this equilibrium measure $$\nu_{\mathrm{eq}}$$. Let us now centre our empirical measure, $$L_N^{(\mathbf{x})} = \nu_{\mathrm{eq}} +  (L_N^{(\mathbf{x})} - \nu_{\mathrm{eq}}) $$. Then 

$$\begin{align*} \frac{ \mathcal{Z}_{N,\Gamma}[V]}{\mathsf{Z}_{N,\gamma}[V]}  &= \mathrm{e}^{\frac{i \beta N^2}{2}\int_{\mathbb{R}^2} a \, \mathrm{d}\nu_{\mathrm{eq}}\otimes \mathrm{d}\nu_{\mathrm{eq}} - i \beta N^2 \int_{\mathbb{R}} \Im V \circ \gamma \, \mathrm{d}\nu_{\mathrm{eq}} + i N \left( 1-  \frac{\beta}{2}\right) \int_{\mathbb{R}}\arg \, \gamma^\prime \, \mathrm{d}\nu_{\mathrm{eq}} }  \\
&\times \mathbb{E}_{N,V}\left[ \mathrm{e}^{\frac{i \beta N^2}{2}\int_{\mathbb{R}^2} a \, \mathrm{d}(L_N^{(\mathbf{x})}-\nu_{\mathrm{eq}})\otimes \mathrm{d}(L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}}) + i \beta N^2 \int_{\mathbb{R}} \mathcal{Q} \, \mathrm{d}(L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}}) + i N \left( 1-  \frac{\beta}{2}\right) \int_{\mathbb{R}}\arg \, \gamma^\prime \, \mathrm{d}(L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}}) }\right] \, .\end{align*}$$

where

$$\mathcal{Q}(x) = \int_{\mathbb{R}}a(x,y) \, \mathrm{d}\nu_{\mathrm{eq}}(y) - \Im V(\gamma(x)) \, .$$

We now notice something remarkable about $$\mathcal{Q}$$, which is that it is *constant* on $$\mathrm{supp}\, \nu_{\mathrm{eq}}$$. In particular

$$\mathcal{Q}^\prime(x) = \Im  \left[ \left( \mathrm{p.v.} \int_{\mathrm{R}} \frac{1}{\gamma(x) - \gamma(y)} \, \mathrm{d}\nu_{\mathrm{eq}}(y) - V^\prime(\gamma(x)) \right) \gamma^\prime(x)\right] \, .$$

This vanishes by the $$S$$-curve condition! By itself, this isn't quite enough. We need to do two things. Firstly, our particles range over all of $$\mathbb{R}$$, but it is convenient to condition them to lie in some bounded open set $$K \supset \mathrm{supp}\, \nu_{\mathrm{eq}}$$, which can be done up to exponentially small error. Secondly, $$\gamma$$ can be taken to be an analytic function on $$\mathrm{supp}\, \nu_{\mathrm{eq}}$$, and it natural to choose our contour to be the analytic extension of this to a neighbourhood $$K$$ of $$\mathrm{supp}\, \nu_{\mathrm{eq}}$$. When you do this, you find the $$S$$-curve condition also extends. Thus $$\mathcal{Q}$$ will be constant on $$K$$ also. Since $$L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}}$$ has zero net mass, 

$$\int_K \mathcal{Q} \, \mathrm{d}(L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}}) = 0 \, .$$

Thus

$$\begin{align*} \frac{ |\mathcal{Z}_{N,\Gamma}[V]|}{\mathsf{Z}_{N,\gamma}[V]}  &=  (1+ \mathcal{O}(\mathrm{e}^{-NC})) \mathbb{E}_{N,V}^K \left[ \mathrm{e}^{\frac{i \beta}{2}\int_{K^2} a \, \mathrm{d}\mathrm{Fluct}_N \otimes \mathrm{d}\mathrm{Fluct}_N  + i \left( 1-  \frac{\beta}{2}\right) \int_{K}\arg \, \gamma^\prime \, \mathrm{d}\mathrm{Fluct}_N }\right] \, .\end{align*}$$

where $$\mathrm{Fluct}_N = N (L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}})$$. We know from an analysis of the real model that $$\mathrm{Fluct}_N(f)$$ is asymptotically Gaussian for a smooth bounded function $$f$$. Thus, morally, what we have is 

$$\mathbb{E} \left[ \mathrm{e}^{i x^\mathsf{T} A x + i y^\mathsf{T} x }\right]$$

where $$A$$ is a deterministic matrix and $$x$$ is normally distributed vector in $$\mathbb{R}^n$$ and $$y$$ is a deterministic vector in $$\mathbb{R}^n$$ (in reality we are dealing with an infinite dimensional case). There is an explicit formula (which can be found in the paper) for the above in terms of $$A$$, $$y$$, and the mean and covariance of $$x$$, and crucially this formula is *non-zero*. This formula permits an infinite dimensional extension to when $$A$$ is an operator on a Hilbert space. In this way $$\frac{ \lvert \mathcal{Z}_{N,\Gamma}[V] \rvert }{\mathsf{Z}_{N,\gamma}[V]}$$ may be bounded from below by a constant. This illustrates what was meant at the beginning when it was said that the $$S$$-curve minimises the oscillations and hence optimises the triangle inequality. This essentially concludes the proof.

**Remark:** This choosing of the contour to be an analytic extension of $$\mathrm{supp}\, \mu_{\mathrm{eq}}$$ should be regarded as the infinite dimensional analogue of the steepest descent contour, since it stabilises the oscillations by eliminating the potentially dangerous term $$\mathrm{e}^{i \beta N^2 \int_{K} \mathcal{Q} \, \mathrm{d}(L_N^{(\mathbf{x})}- \nu_{\mathrm{eq}})}$$.

**Remark:** Before finishing something must be said about the interpolation. Everything said above is for a fixed potential, but in fact we need to find a smooth *family* of potentials for which all our hypotheses (in particular, one cut regular) are true. This is not entirely trivial, unlike for the case of real potentials where there is a nice interpolation. The solution we found is the following. Let $$\gamma : [0,1] \to \mathbb{C}$$ be an analytic curve such that $$\gamma(0) = 0$$ and $$\gamma(1) = 1$$. Then let

$$\gamma_t(x) = \frac{\gamma(tx)}{\gamma(t)}\, .$$

Clearly $$\gamma_1(x) = \gamma(x)$$ and $$\lim_{t \downarrow 0}\frac{\gamma(tx)}{\gamma(t)} = x$$. Furthermore $$\gamma_t(0) = 0$$ and $$\gamma_t(1) = 1$$, and $$\gamma_t$$ is injective since $$\gamma$$ is. Hence $$\gamma_t$$ is an interpolating family of curves from $$\gamma$$ to a straight line, with the same endpoints. The key observation, now, is that for a quadratic potential the $$S$$-curve is a straight line. Hence if we find the potential $$V_t$$ associated to the $$S$$-curve $$\gamma_t$$, we will have found our interpolating family. This construction easily generalises to curves where the endpoints are not $$0$$ and $$1$$. $$\triangle$$
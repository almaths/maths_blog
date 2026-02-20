---
layout: post
title:  "Skew-orthogonal polynomials as multiple-orthogonal polynomials"
date:   2024-08-19 12:00:00 +0100
permalink: /symplectic-skew-orthogonal/
katex: True
---

In classical random matrix theory one is interested in studying real-symmetric, complex-Hermitian or quaternion self-dual matrices with density given by

$$\frac{1}{Z_n} \mathrm{e}^{-\mathrm{tr} \, V(M)} \, \mathrm{d}M$$

where $$V: \mathbb{R} \longrightarrow \mathbb{R}$$ is a "potential" such that $$V(x)$$ tends to $$+\infty$$ as $$x \to \pm \infty$$ sufficiently fast such that the measure is normalisable. Such ensembles are often called "orthogonal", "unitary" or "symplectic" respectively because the measure is invariant under conjugations from $$\mathrm{O}(n)$$, $$\mathrm{U}(n)$$ or $$\mathrm{Sp}(n)$$. Under this distribution the eigenvalues of $$M$$ have the joint distribution

$$\rho_n(x_1, \dots, x_n) = \frac{1}{\tilde{Z}_n} \prod_{1 \leq i < j \leq n} |x_i - x_j|^\beta \prod_{k=1}^n \mathrm{e}^{- V(x_k)}$$

where $$\beta = 1, 2, 4$$ corresponds respectively to the real symmetric, complex-Hermitian and quaternion self-dual matrices. Of the three ensembles $$\beta = 2$$ has been by far the most studied, for reasons I shall mention in a moment; though from the point of view of quantum chaos $$\beta = 1,4$$ actually tends to be more relevant since $$\beta = 1,4$$ corresponds to quantum mechanical systems with time-reversal symmetry, while $$\beta = 2$$ corresponds to quantum mechanical systems without time-reversal symmetry, and the latter tends to be relatively rare in nature.

The greater attention given to $$\beta = 2$$ comes from its greater mathematical elegance. Such ensembles are related in a beautiful way to the theory of orthogonal polynomials and determinantal point processes, which in turn connects it to various "integrable" models in probability theory and the theory of integrable dynamical systems, e.g. the Toda lattice, Painlevé equations. For $$\beta = 1,4$$ instead of orthogonal polynomials, one has *skew-orthogonal polynomials*; that is, the inner product is replaced with a skew-symmetric bilinear form. Skew-orthogonal polynomials are the basis of polynomials in which this bilinear form is put into $$2 \times 2$$ skew-symmetric blocks. Instead of a determinantal point process, eigenvalues of $$\beta = 1,4$$ ensembles obey a Pfaffian point process whose kernel is expressed as a sum over skew-orthogonal polynomials (see this  [article](https://www.brandeis.edu/mathematics/people/adlerpublications/chapter5.pdf) for an accessible introduction). Unlike orthogonal polynomials, however, the theory of skew-orthogonal polynomials is less well-developed. Until recently, there was no analogue of the 3 term recurrence relation or Christoffel-Darboux formula, the latter of which is essential for the asymptotic analysis of correlation functions for $$\beta = 2$$. Furthermore, orthogonal polynomials possess a representation in terms of a Riemann-Hilbert problem ([Fokas, Its & Kitaev, 1992](https://link.springer.com/article/10.1007/BF02096594)), which has allowed a rigorous asymptotic analysis by the Deift-McLaughlin-Kriecherbauer-Venakides-Zhou (DMKVZ) scheme.

To study $$\beta = 1,4$$ ensembles asymptotically ($$n \to \infty$$) past work has involved making some reduction to the $$\beta = 2$$ case, for example, either by expressing the skew-orthogonal polynomials in terms of corresponding orthogonal polynomials, or by expressing the "pre-kernel" of the Pfaffian point process in terms of the "kernel" of the corresponding determinantal point process plus a finite rank correction. This approach works for potentials $$V$$ such that $$V^\prime$$ is a rational function. This latter approach due to Widom has allowed universality theorems to be proven for $$\beta = 1,4$$ by taking advantage of known results for $$\beta = 2$$. My impression is that Widom developed this method out of a dissatisfaction with the skew-orthogonal polynomial method. Ghosh in his [book](https://bookstore.ams.org/view?ProductCode=CRMM/28) on skew-orthogonal polynomials described Widom's method as "rigorous but cumbersome." Although this method has been successful in rigorously proving universality theorems, this method has struck me as somewhat "artificial" in that such $$\beta = 1,4$$ ensembles are exactly solvable models, but yet the asymptotic analysis does not treat them as such, instead treating them as a perturbation of a solvable model. A nice introduction to this method is given in the [book](https://www.ams.org/books/cln/018/) of Deift and Gioev.

This motivated me to investigate if these models could be approached in a way that respects their integrability. In this way we would have a completely "integrable" picture from start to finish. The result of my investigations is that for $$\beta = 4$$ and polynomial potentials, the answer is yes. This is the content of [my recent paper](https://arxiv.org/abs/2306.14107), [published](http://sigma-journal.com/2024/076/) in the journal SIGMA. In this work I found a Riemann-Hilbert representation for $$\beta = 4$$ skew-orthogonal polynomials. This was done by showing that $$\beta = 4$$ skew-orthogonal polynomials, for polynomial potentials $$V$$, can be viewed as a species of *multiple-orthogonal polynomials*. This work bears some similarity to work done by [Pierce (2006)](https://arxiv.org/abs/math/0610592) in which he showed that, again for polynomial potentials, $$\beta = 1$$ skew-orthogonal polynomials can be viewed as multiple-orthogonal polynomials and so admit a Riemann-Hilbert representation. In my work I also prove a Christoffel-Darboux-type formula, and so if a nonlinear steepest descent analysis could be carried out, we would have a new (and perhaps more elegant) proof of universality for $$\beta = 4$$.

In this blog post I want to sketch how this connection to multiple-orthogonality is achieved. To do this we need to recall the definition of skew-orthogonal polynomials. We start with the bilinear form

$$\langle P,Q \rangle_4 = \frac{1}{2} \int_{\mathbb{R}} (P(x)Q^\prime(x)-P^\prime(x)Q(x))\mathrm{e}^{-2V(x)}\, \mathrm{d}x $$

where $$P$$ and $$Q$$ are polynomials. $$\langle \cdot , \cdot \rangle_4$$ is clearly skew-symmetric. Somewhat less obviously $$\langle \cdot , \cdot \rangle_4$$ is non-degenerate in the following sense. If we define $$X^i(x) = x^i$$ then the matrix $$\left(\langle X^{i-1} , X^{j-1} \rangle_4 \right)_{i,j=1}^{2k} $$ is invertible for all $$k \in \mathbb{N}$$. This is proven in Proposition 2.1 of the paper. Moving forward, we say the sequence of monic polynomials $$\{ P_k \}_{k \in \mathbb{N}}$$, where $$P_k$$ has degree exactly $$k$$, is *skew-orthogonal* (of symplectic type) if

(1) $$\langle P_{2n}, P_{2m} \rangle_4 = \langle P_{2n+1}, P_{2m+1} \rangle_4 = 0$$ for all $$n,m \in \mathbb{N}$$.

(2) $$\langle P_{2n}, P_{2m+1} \rangle_4 = 0$$ for all $$n \neq m$$, for $$n,m \in \mathbb{N}$$.

This amounts effectively to finding a change of basis which puts the matrix $$\left(\langle X^{i-1} , X^{j-1} \rangle_4 \right)_{i,j=1}^{2k} $$ in block-diagonal form, with each block being a $$2 \times 2$$ skew-symmetric matrix. This sequence is not quite unique because we can replace $$P_{2n+1} \to P_{2n+1} + \alpha P_{2n}$$ for any $$\alpha \in \mathbb{R}$$. This degree of freedom may be fixed by requiring that the next-to-leading coefficient of $$P_{2n+1}$$ have a specified value. Given this extra constraint, the sequence is unique.

To make the connection to multiple-orthogonal polynomials, first observe that we can write

$$\langle P,Q \rangle_4 = \frac{1}{2} \int_{\mathbb{R}} P(x)\mathrm{e}^{-V(x)} \frac{\mathrm{d}}{\mathrm{d}x}\left( Q(x) \mathrm{e}^{-V(x)}\right) \, \mathrm{d}x .$$

Let us consider the even case only. Let $$V(x) = \gamma x^{D}+ \mathcal{O}(x^{D-1})$$ as $$x \to \infty$$ be our polynomial potential. Define the "dual" skew-orthogonal

$$\Psi_{2n}(x) = -\frac{1}{\gamma D} \mathrm{e}^{V(x)} \frac{\mathrm{d}}{\mathrm{d}x}\left( P_{2n}(x)\mathrm{e}^{-V(x)}\right) .$$

Note that $$\Psi_{2n}$$ is monic of degree $$2n+D-1$$ and $$P_{2n}$$ can be reconstructed from $$\Psi_{2n}$$ by integration. Then the requirement that $$\langle X^k,P_{2n} \rangle_4 = 0$$ for $$k = 0, \dots, 2n-1$$ is equivalent to

$$\int_{\mathbb{R}} x^k \Psi_{2n}(x)\mathrm{e}^{-2V(x)}\, \mathrm{d}x  = 0$$

for $$k=0, \dots, 2n-1$$. Note that $$\langle X^{2n},P_{2n} \rangle_4 = 0$$ however we get this "for free" from the other constraints by the skew-symmetry of the inner product. Thus $$\Psi_{2n}$$ is rather close to being an orthogonal polynomial. However notice that these conditions are not enough to determine $$\Psi_{2n}$$ since $$\Psi_{2n}$$ has $$2n+D-1$$ parameters while the above yields only $$2n$$ constraints. We thus need to find extra constraints. These extra constraints come from the fact that the map $$P \mapsto \mathrm{e}^{V} \frac{\mathrm{d}}{\mathrm{d}x}(P \mathrm{e}^{-V})$$ is not a surjective map within the space of polynomials. We need to characterise the image of this map. To do this, observe that

$$P_{2n}(x) = \gamma D \mathrm{e}^{V(x)} \int_{x}^{+\infty} \Psi_{2n}(y) \mathrm{e}^{-V(y)}\, \mathrm{d}y .$$

For an arbitrary monic polynomial $$Q$$ of degree $$2n+D-1$$,

$$x \mapsto \mathrm{e}^{V(x)} \int_{x}^{+\infty} Q(y) \mathrm{e}^{-V(y)}\, \mathrm{d}y$$

is not typically a polynomial. The above function is clearly entire, and so to determine if it is a polynomial we need only look at its $$x \to \infty$$ asymptotics in all possible directions in $$\mathbb{C}$$. Essentially, $$x \mapsto \mathrm{e}^{V(x)} \int_{x}^{+\infty} Q(y) \mathrm{e}^{-V(y)}\, \mathrm{d}y$$ will have different asymptotics in different sectors in $$\mathbb{C}$$ with Stokes jumps between these different sectors. The result is a polynomial if all these Stokes jumps are zero, which yields the "missing" constraints. The precise statement is the following.

**Proposition:** Let $$C_k = \mathrm{e}^{\frac{2\pi \mathrm{i} k}{D} }[0,+\infty)$$ for $$k = 0, \dots, D-1$$. $$C_0$$ has outgoing orientation (away from the origin) and $$C_1, \dots, C_{D-1}$$ have incoming orientation (towards the origin). Then define $$\Gamma_k = C_0 \cup C_k$$ for $$k = 1, \dots, D-1$$, with orientation carried over.

Then for any polynomial $$Q$$

$$x \mapsto \mathrm{e}^{V(x)} \int_{x}^{+\infty} Q(y) \mathrm{e}^{-V(y)}\, \mathrm{d}y$$

is a polynomial if and only if

$$\int_{\Gamma_k} Q(x) \mathrm{e}^{-V(x)}\, \mathrm{d}x = 0$$

for all $$k = 1, \dots, D-1$$. $$\triangle$$

The conditions then uniquely determine $$\Psi_{2n}$$.

**Corollary:** Let $$Q$$ be a monic polynomial of degree $$2n+D-1$$. Then $$Q = \Psi_{2n}$$ if and only if

$$\int_{\mathbb{R}}Q(x) x^j \mathrm{e}^{-2V(x)}\, \mathrm{d}x = 0 $$

for all $$j = 0,\dots, 2n-1$$, and

$$\int_{\Gamma_k} Q(x) \mathrm{e}^{-V(x)}\, \mathrm{d}x = 0$$

for all $$k = 1, \dots, D-1$$. $$\triangle$$

The above corollary effectively shows that $$\Psi_{2n}$$ can be thought of as a multiple-orthogonal polynomial of Type II. From these results the whole machinery of multiple-orthogonal polynomials and their associated Riemann-Hilbert problems can be brought to bear.

An additional surprising fact is that our original skew-orthogonal polynomial $$P_{2n}$$ can be viewed as multiple-orthogonal of Type I. Consider the problem of finding a monic polynomial $$P$$ of degree $$2n$$ and a collection of constants $$c_1, \dots, c_{D-1}$$ such that

$$\int_{\mathbb{R}} P(x) x^k \mathrm{e}^{-2V(x)}\, \mathrm{d}x + \sum_{j=1}^{D-1} c_j \int_{\Gamma_j} x^k \mathrm{e}^{-V(x)}\, \mathrm{d}x = 0$$

for all $$k = 0 , \dots, 2n+D-2$$. By taking linear combinations of these equations, we can of course replace $$x^k$$ by any polynomial of degree $$\leq 2n+D-2$$, in particular we can replace it with $$\Psi_0, \dots, \Psi_{2n-1}$$. This gives

$$\int_{\mathbb{R}} P(x) \Psi_{k}(x) \mathrm{e}^{-2V(x)}\, \mathrm{d}x = 0$$

for $$k = 0, \dots, 2n-1$$, and hence $$P=P_{2n}$$. The constants $$c_1, \dots, c_{D-1}$$ can also be computed. From this we see that $$P_{2n}$$ is multiple-orthogonal of Type I.

The situation for the odd polynomials $$P_{2n+1}$$ and their duals $$\Psi_{2n+1}$$ is a bit more complicated. The latter permit an interpretation in terms of Type II multiple orthogonality, but I could not find a way to represent $$P_{2n+1}$$ in terms of Type I multiple-orthogonality (though a mixture of both Type I and Type II conditions works). However a surprising result is that, from the point of view of random matrix theory, it is completely sufficient to study the even case. By which I mean, the even case yields a Riemann-Hilbert problem from which a Christoffel-Darboux-type formula can be derived, giving the pre-kernel of the Pfaffian point process, and so the RHP for the odd problem turns out to be unnecessary.
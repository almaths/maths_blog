---
layout: post
title:  "Making sense of eigenvalues of self-dual quaternion matrices"
date:   2023-07-24 12:00:00 +0000
permalink: /quaternion-eigenvalues/
katex: True
---

In my [recent work](https://arxiv.org/abs/2306.14107) I made a connection between the theory of self-dual quaternion random matrices and Riemann-Hilbert problems. As part of the background of this research, I needed to revisit the theory of self-dual quaternion random matrices, in particular the question how to make sense of the *eigenvalues* of such matrices. This is not entirely self-explanatory given quaternions do not commute. In this post I hope to give an accessible explanation of this.

First let us recall basic facts about quaternions. The algebra of quaternions $$\mathbb{H}$$ is the real span of 4 linearly independent elements $$1, e_1, e_2, e_3$$ with the relations

$$e_1^2 = e_2^2 = e_3^2 = -1 $$

$$e_1 e_2 = e_3 \quad \text{ etc. by cyclic permutations}$$

$$e_i e_j = -e_j e_i \quad \text{ for } i \neq j$$

It is convenient to identify these with $$2 \times 2$$ matrices

$$\begin{align*}1 \simeq \mathbb{I} = \left( \begin{matrix} 1 & 0 \\ 0 & 1 \end{matrix} \right), & & e_1  \simeq \left( \begin{matrix} i & 0 \\ 0 & -i \end{matrix} \right) \\
e_2  \simeq \left( \begin{matrix} 0 & 1 \\ -1 & 0 \end{matrix} \right), & & e_3  \simeq \left( \begin{matrix} 0 & i \\ i & 0 \end{matrix} \right)\end{align*}$$


In what follows it will be useful to complexify the quaternions $$\mathbb{H}_{\mathbb{C}}$$ so that for $$Q \in \mathbb{H}_{\mathbb{C}}$$

$$Q = q_0 \mathbb{I} + q_1 e_1 + q_2 e_2 + q_3 e_3 \quad (q_i \in \mathbb{C})$$

**Definition:** The *dual* of a quaternion $$Q= q_0 \mathbb{I} + q_1 e_1 + q_2 e_2 + q_3 e_3 \in \mathbb{H}_\mathbb{C}$$ is


$$Q^\mathsf{D} = q_0 \mathbb{I} - q_1 e_1 - q_2 e_2 - q_3 e_3$$

Note that $$Q \mapsto Q^\mathsf{D}$$ is a $$\mathbb{C}$$-linear (and not conjugate linear) operation.

**Lemma:** Using our $$2 \times 2$$ matrix representation of a quaternion $$Q \in \mathbb{H}_{\mathbb{C}}$$ we may write the dual

$$Q^\mathsf{D} = -e_2 Q^\mathsf{T} e_2$$

**Proof:** Straightfoward calculation. $$\square$$

**Definition:** The *adjoint* of a quaternion $$Q= q_0 \mathbb{I} + q_1 e_1 + q_2 e_2 + q_3 e_3 \in \mathbb{H}_\mathbb{C}$$ is

$$Q^\dagger = \overline{q_0} \mathbb{I} - \overline{q_1} e_1 - \overline{q_2} e_2 - \overline{q_3} e_3$$

Note that $$Q \mapsto Q^\dagger$$ is a conjugate-linear operation and given our matrix representation it is exactly the conjugate transpose of the matrix $$Q$$.

**Corollary:** A quaternion $$Q$$ is *real* (has real coefficients) if and only if $$Q^\dagger = Q^\mathsf{D} $$, i.e. $$Q^\dagger = -e_2 Q^\mathsf{T} e_2$$. Equivalently, a $$2 \times 2$$ matrix $$Q$$ is in the *real* span of $$\mathbb{I}, e_1, e_2, e_3$$ if and only if $$Q^\dagger = -e_2 Q^\mathsf{T} e_2$$. $$\triangle$$

We can now see the advantage of introducing $$\mathbb{H}_\mathbb{C}$$ even though we are really only interested in $$\mathbb{H}$$. Given an $$n \times n$$ (real) quaternion matrix $$\mathcal{M}$$ we identify this with a $$2n \times 2n$$ matrix $$M$$, and the condition that $$\mathcal{M}_{ij} = \mathcal{M}_{ji}^\mathsf{D}$$ becomes the requirement that

$$M = M^\mathsf{D} = M^\dagger$$

where $$M^\mathsf{D} = - J M^\mathsf{T}J$$ for $$J = \underbrace{e_2 \oplus \dots \oplus e_2}_{n \text{ times}}$$.

**Remark:** Define the non-degenerate skew-symmetric bilinear form $$\Omega : \mathbb{C}^{2n} \times \mathbb{C}^{2n} \to \mathbb{C}$$ by

$$\Omega(x,y) = x^\mathsf{T} J y$$

Then $$M = M^\mathsf{D}$$ is equivalent to $$\Omega(Mx,y) = \Omega(x,My)$$ for all $$x,y \in \mathbb{C}^{2n}$$. $$\triangle$$

**Definition:** The (non-compact) symplectic group $$\mathrm{Sp}(n)$$ is the group of $$2n \times 2n$$ matrices $$U$$ for which $$\Omega(Ux,Uy) = \Omega(x,y)$$ for all $$x,y \in \mathbb{C}^{2n}$$. The (compact) symplectic group is $$\mathrm{USp}(n)=\mathrm{Sp}(n) \cap \mathrm{U}(2n)$$. $$\triangle$$

It is easily seen that for $$U \in \mathrm{USp}(n)$$, $$U^{-1}= U^\mathsf{D} = U^\dagger$$, so that $$U$$ may be thought of as an $$n \times n$$ matrix with real quaternion entries whose dual is its inverse. Note that $$\mathrm{USp}(n)$$ is exactly the group which, acting by conjugation, preserves (real) quaternion self-duality.

**Proposition (Kramers' degeneracy):** Let $$M = M^\mathsf{D}$$ be a $$2n \times 2n$$ matrix. Then the characteristic polynomial of $$M$$ is an exact square. In particular, $$M$$ has generically $$n$$ eigenvalues each of multiplicity $$2$$.

**Proof:** Because $$M = M^\mathsf{D}$$ we have that $$(JM)^\mathsf{T} = - JM$$ and so

$$\det(\zeta \mathbb{I}- M) = \det(\zeta J- JM) = \left( \mathrm{pf}(\zeta J- JM)\right)^2 $$

for $$\zeta \in \mathbb{C}$$ and $$\mathrm{pf}$$ being the Pfaffian. Here we have used that $$\det J = 1$$. $$\square$$

**Remark:** Many works, including e.g. the textbooks of M. L. Mehta (*Random Matrices*) and P. Forrester (*Log-Gases and Random Matrices*), prefer to work with a so-called "quaternion determinant." Given a self-dual $$n \times n$$ quaternion matrix $$\mathcal{M}$$ with $$2n \times 2n$$ representative $$M$$, we define the *quaternion determinant*

$$\mathrm{Qdet}(\mathcal{M}) = \mathrm{pf}(JM)$$

Surprisingly, there is a theorem due Dyson (see Theorem 5.1.2 of Mehta's textbook) that shows that $$\mathrm{Qdet}$$ admits a Laplace-type formula in terms of a sum over permutations (ibid, Equation 5.1.5). All of this presumes that the matrix $$M$$ is self-dual, as far as I understand $$\mathrm{Qdet}$$ is not defined for non-self-dual matrices. $$\triangle$$

Finally, to conclude our discussion, we must give meaning to the notion of *diagonalising* quaternion self-dual matrices. Let $$\mathcal{M}$$ be an $$n \times n$$ self-dual quaternion matrix and $$M = M^\dagger = M^\mathsf{D}$$ be its $$2n \times 2n$$ representative. We aim to show that $$M$$ may be diagonalised by an element of $$\mathrm{USp}(n)$$. Let us assume for simplicity of exposition that $$M$$ has exactly $$n$$ (distinct) eigenvalues $$\lambda_1, \dots, \lambda_n \in \mathbb{R}$$ each of multiplicity $$2$$. Let $$v_k \in \mathbb{C}^{2n}$$ be an eigenvector, $$\| v_k \| = 1$$, with eigenvalue $$\lambda_k$$.

$$M v_k = \lambda_k v_k$$

By self-duality, $$w_k := J \overline{v_k}$$ is also an eigenvector with $$\lambda_k$$. $$w_k$$ and $$v_k$$ are distinct eigenvectors since $$\| w_k \| = 1$$ and $$\langle w_k , v_k \rangle = w_k^\dagger v_k = v_k^\mathsf{T} J v_k = 0$$. Then define the matrix

$$U = \left( \begin{matrix} \vert & \vert &\dots &\vert & \vert \\
v_1 & w_1 &\dots & v_n & w_n \\
\vert & \vert &\dots &\vert & \vert \end{matrix} \right)$$

From the construction it is clear that

$$U^{-1} M U = \mathrm{diag}(\lambda_1, \lambda_1, \dots, \lambda_n, \lambda_n)$$

so $$U$$ diagonalises $$M$$. Furthermore we claim $$U \in \mathrm{USp}(n)$$. This can be seen from the following. Firstly, since the columns of $$U$$ are orthonormal with respect to the standard Hermitian inner product on $$\mathbb{C}^{2n}$$, $$U$$ must be unitary ($$U^{-1} = U^\dagger$$). Secondly, again by construction $$J U J = - \overline{U} = -U^\mathsf{-T}$$, and hence $$U^\mathsf{D} = U^{-1}$$. This completes the proof that $$U \in \mathrm{USp}(n)$$.
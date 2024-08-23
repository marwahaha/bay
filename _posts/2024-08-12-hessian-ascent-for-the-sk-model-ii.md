---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model II"
date:   2024-08-12 11:15:40
blurb: "Optimizing the SK model via free probability and convex duality"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)

This is the second blog past in a **3-part** series on a recent result [[JSS24]](https://arxiv.org/abs/2408.02360) by [David](https://davidjekel.com/), [Jonathan](https://jshi.science/) and me. The first blog post can be found [here]().

A quick recap of what we accomplished in the first blog post:
1. We briefly went over the various representations of the Parisi formula, and then derived the Hessian ascent algorithm from the generalized TAP free energy.
2. We then motivated the two main statements we need to prove to demonstrate the success of our algorithm. We concluded by developing a _primal_ theory for the Parisi PDE and AC SDE, which will be quintessential in the proofs of these statements.
<br>

#### Table of Contents
1. [Proof Sketch](#proof-sketch)
   * [Analyzing the TAP-corrected Hessian via free probability](#analyzing-the-tap-corrected-hessian-via-free-probability)
   * [Empirical distribution of the coordinates of the iterates](#empirical-distribution-of-the-coordinates-of-the-iterates)
   * [Fluctuations of the generalized TAP free energy under fRSB](#fluctuations-of-the-generalized-tap-free-energy-under-frsb)
2. [Footnotes](#footnotes)
<br>

## Proof Sketch
### Analyzing the TAP-corrected Hessian via free probability
We now begin sketching the proof overview for the first main statement. As stated in the previous post, the goal is to create a matrix $$Q(t_j, \sigma_j) $$ for every step $$j \in [K] $$ with time $$t_j = j\eta $$ that, with high probability,
1. Projects into the top-eigenspace of $$\nabla^2 \mathsf{TAP} :=\sqrt{2}\beta\,\mathsf{GOE}(n) - D(t_j,\sigma_j) $$, and
2. Has "well-behaved" diagonal entries, so that we can assert convergence of the empirical distribution of the coordinates of every iterate of the algorithm to the primal version of the AC SDE in [[2.2]](#empirical-distribution-of-the-coordinates-of-the-iterates).

We proceed in chronological fashion, first giving a construction of $$Q $$ and showing it projects into the top-eigenspace, and then proving that this construction leads to diagonal entries with the desired behavior.

#### Projecting into the top-eigenspace
To project into the top-part of the spectrum of $$\nabla^2 \mathsf{TAP} $$, we will create an operator whose eigenspectrum is subtracted from the (idealized) operator norm of $$\nabla^2 \mathsf{TAP} $$ and then inverted. Intuitively, this should make the large eigenvalues more important than the small eigenvalues, while preserving the eigenvectors. Sans normalization and a projection away from the current iterate $$\sigma_j $$, therefore, a candidate choice is,

$$ \begin{equation} P(D)^2 := b\left(b^2\mathsf{Id} + (a(D) - (2\beta A_{sym} - D))\right)^{-1}\,, \end{equation} $$

where $$a(D) $$ is (approximately) the maximum value of the spectrum of $$\sqrt{2}\beta S - D(j\eta,\sigma_j) $$, and $$S $$ is an operator whose empirical eigenspectrum is the semi-circular law. As we shall see, $$b $$ will be chosen to be sufficiently small (but non-zero). It is not hard to see that as $$b \to 0 $$, the eigenvalues farthest from $$a(D) $$ will receive the smallest weight, and the ones closest to $$a(D) $$ will receive the largest. Two quick points:
* This choice of $$a(D) $$ implicitly suggests that we will compare $$2\beta A_{sym} $$ with $$\sqrt{2}\beta S $$; indeed, this is one of the key parts in the proof strategy.
  * We will use Gaussian concentration to compare $$P^2(D) $$ with its expected value (under the Gaussian randomness), and show that this deviation is desirably small.
  * We will then show that the expected value $$\mathbb{E}[P(D)^2] $$ is close to the "idealized" value of
  $$ \begin{equation} P^2_{S}(D) := b\left(b^2\mathsf{Id} + (a(D) - (\sqrt{2}\beta S - D))\right)^{-1}\,. \end{equation} $$

* The scalar $$b $$ will be chosen to be sufficiently small so that it does not cause $$P^2(D) $$ to deviate too much from projecting into the desired top-eigenspace. However, it will be non-zero all the same, so that the operator $$P^2 (D) $$ continues to be well-defined.

Let us now write a brief and informal statement that summarizes two (of the three) key qualities the operator $$P(D)^2 $$[^1] will have. In the statement below, we will assume that $$\frac{2\beta^2}{n}\mathsf{Tr}[D^{-2}] = 1$$.

$$ \begin{align*} \text{[Idealized operator norm]:}& \quad\quad \| \sqrt{2}\beta S - D\|_{\mathsf{op}} = \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}]\,. \\
\text{[Large-overlap with top-eigenspace]:} & \quad\quad \frac{1}{n}\left\langle P(D)^2,(\widetilde{a(D)}-(2\beta A_{sym} - D))^2\right\rangle \le \frac{O_\beta(\mathsf{Tr}[D^{-4}])}{n^{1.03}}\,. \end{align*} $$

The notation $$O_\beta(\cdot) $$ suppresses the dependencies on $$\beta = O(1/\epsilon)$$ for terms that decay sufficiently fast in $$n $$. For technical reasons, $$\widetilde{a(D)} $$ is not _exactly_ $$\|\sqrt{2}\beta S - D\|_{\mathsf{op}} $$, but a _small_ perturbation around it whereby,

$$ \begin{equation} \left| \widetilde{a(D)} - \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}] \right| \le \frac{O_\beta\left(\mathsf{Tr}[D^{-3}]\right)}{n^{1.02}}\,. \end{equation} $$

By putting the two statements above together with the observation that $$P(D)^2 $$ projects into the top-eigenspace of $$\sqrt{2}\beta S - D $$, one can reasonably believe that

$$ \begin{equation} \frac{1}{n} \|P(D) -  \Pi_{\delta'}\left(2\beta A_{sym} - D\right)\|_2^2 \le \mathsf{small} \,, \end{equation} $$

where $$\Pi_{\delta'}\left(\sqrt{2}\beta S - D\right) $$ is a projector into the top-$$\delta' n$$ dimensional eigenspace of the spectrum of $$2\beta A_{sym} - D $$.

We will now give overviews of the proofs for the two statements above. Giving the actual detailed proofs would require introducing various technical facts and auxiliary lemmata about non-commutative $$\mathcal{L}^p $$ spaces, the Cauchy-Stieljtes transform, the $$R $$-transform and analytic subordination. Consequently, to ensure that the broader TCS audience can follow along, we will try to sketch the main ideas with minimal possible reference to the underlying formalism[^2].


#### Approximating the diagonal of the projector
We now formally choose the distortion parameter $$b = \beta n^{-.01} $$. With this choice, we infer that the diagonal entries of $$P(D)^ 2$$ behave desirably, with high probability.

$$ \begin{equation} \text{[Diagonal entries]} \end{equation} $$

### Empirical distribution of the coordinates of the iterates

### Fluctuations of the generalized TAP free energy under fRSB

#### FOOTNOTES

[^1]: Those familiar with resolvent formalism will immediately notice that $$P^2(D) $$ is nothing but the (negative of the) imaginary part of the _resolvent_ $$g_{2\beta A_{sym} - D} $$, evaluated at a purely imaginary point $$ib $$, which is a slight perturbation from $$0 $$. Since, for a compactly supported measure, which we expect the empirical eigenspectrum of $$\sqrt{2}\beta A_{sym} - D $$ to converge to, the Cauchy transform (which is equivalent to the resolvent) is invertible everywhere except at poles, this is a loss-less way of studying analytic properties of the spectrum.

[^2]: Shortly, I will introduce the relevant complex-analytic basics and the relevant free-probability notions in a separate blog-post. In the very same blog-post, I will provide sufficient background that somebody can actually follow most of the proofs in detail (_except_ for the proof which reasons about the projection of $$P(D)^2 $$ onto the diagonal sub-algebra, which I will introduce in yet another post). For the impatient reader, [[Section 4.1, JSS24](https://arxiv.org/pdf/2408.02360)] should also act as an intense, but reasonably accessible proxy.

[^3]: When David, Jonathan and I were working on the project in the early days, we were thinking of this convex duality via the lens of mirror maps. It turns out that, because of the underlying convex duality, there is a connection to be made here through the information geometry of the underlying Hessian being akin to a (flat) Bregmannian manifold. However, this is beyond the scope of this post and something we are investigating in ongoing work.

[^4]: In a previous [post](https://juspreetsandhu.me/2022/02/08/the-sk-model-i#ruelle-probability-cascades), I described a _part_ of the RPC construction and the details afforded there are gentle and sufficient enough to understand the main ingredients in their construction, as well as their purpose. A slightly more detailed overview of the construction, along with how exactly the representation gets used in the [Hopf-Cole transform]() to solve the Parisi PDE for _atomic_ measures is provided in [[Appendix C, JSS24]](https://arxiv.org/pdf/2408.02360) and can be read by the interested reader. The estimates for the derivatives are proved between [[Lemma 2.12 & Lemma 2.13, JSS24]](https://arxiv.org/pdf/2408.02360) and stated in [[Proposition 2.11, JSS24]](https://arxiv.org/pdf/2408.02360).  

[^5]: For converting Lipschitz bounds into estimates of how much error propagates over a period of time, [Gronwall's inequality](https://en.wikipedia.org/wiki/Gr%C3%B6nwall%27s_inequality) is an indispensable tool that we use in Sections 2 & 5 of the paper. The combination of Ito's lemma & Gronwall's inequality allows us to, in fact, more or less have a mechanistic procedure for converting Lipschitz estimates into error bounds of various sorts.

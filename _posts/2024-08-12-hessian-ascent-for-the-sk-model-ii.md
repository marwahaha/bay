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
2. We then motivated the two main statements we need to prove to demonstrate the success of our algorithm. We concluded by developing a _primal_ theory for the Parisi PDE and AC SDE, which will be quintessential at various points in the main arguments.
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

$$ \begin{equation} P(D)^2 := b\left(b^2\mathsf{Id} + (\widetilde{a(D)} - (2\beta A_{sym} - D))^2\right)^{-1}\,, \end{equation} $$

where $$\widetilde{a(D)} $$ is (approximately) the maximum value of the spectrum of $$\sqrt{2}\beta S - D(j\eta,\sigma_j) $$, and $$S $$ is an operator whose empirical eigenspectrum is the semi-circular law. As we shall see, $$b $$ will be chosen to be sufficiently small (but non-zero). It is not hard to see that as $$b \to 0 $$, the eigenvalues farthest from $$a(D) $$ will receive the smallest weight, and the ones closest to $$\widetilde{a(D)} $$ will receive the largest. Two quick points:
* This choice of $$\widetilde{a(D)} $$ implicitly suggests that we will compare $$2\beta A_{sym} $$ with $$\sqrt{2}\beta S $$; indeed, this is one of the key parts in the proof strategy.
  * We will use Gaussian concentration to compare $$P^2(D) $$ with its expected value (under the Gaussian randomness), and show that this deviation is desirably small.
  * We will then show that the expected value $$\mathbb{E}[P(D)^2] $$ is close to the "idealized" value of

  $$ \begin{equation} P_{S}(D)^2 := b\left(b^2\mathsf{Id} + (\widetilde{a(D)} - (\sqrt{2}\beta S - D))^2\right)^{-1}\,. \end{equation} $$

* The scalar $$b $$ will be chosen to be sufficiently small so that it does not cause $$P^2(D) $$ to deviate too much from projecting into the desired top-eigenspace. However, it will be non-zero all the same, so that the operator $$P^2 (D) $$ continues to be well-defined.

Let us now write a brief and informal statement that summarizes two (of the three) key qualities the operator $$P(D)^2 $$ will have[^1]. In the statement below, we will assume that $$\frac{2\beta^2}{n}\mathsf{Tr}[D^{-2}] = 1$$.

$$ \begin{align*} \text{[Idealized operator norm]:}& \quad\quad \| \sqrt{2}\beta S - D\|_{\mathsf{op}} = \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}]\,. \\
\text{[Large-overlap with top-eigenspace]:} & \quad\quad \frac{1}{n}\left\langle P(D)^2,(\widetilde{a(D)}-(2\beta A_{sym} - D))^2\right\rangle \le \frac{O_\beta(\mathsf{Tr}[D^{-4}])}{n^{1.03}}\,. \end{align*} $$

The notation $$O_\beta(\cdot) $$ suppresses the dependencies on $$\beta = O(1/\epsilon)$$ for terms that decay sufficiently fast in $$n $$. For technical reasons, $$\widetilde{a(D)} $$ is not _exactly_ $$\|\sqrt{2}\beta S - D\|_{\mathsf{op}} $$, but a _small_ perturbation around it whereby,

$$ \begin{equation} \left| \widetilde{a(D)} - \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}] \right| \le \frac{O_\beta\left(\mathsf{Tr}[D^{-3}]\right)}{n^{1.02}}\,. \end{equation} $$

By putting the two statements above together with the observation that $$P(D)^2 $$ projects into the top-eigenspace of $$\sqrt{2}\beta S - D $$, one can reasonably believe that

$$ \begin{equation} \frac{1}{n} \|P(D) -  \Pi_{\delta'}\left(2\beta A_{sym} - D\right)\|_F^2 \le \mathsf{small} \,, \end{equation} $$

where $$\Pi_{\delta'}\left(\sqrt{2}\beta S - D\right) $$ is a projector into the top-$$\delta' n$$ dimensional eigenspace of the spectrum of $$2\beta A_{sym} - D $$.

We will now give overviews of the proofs for the two statements above. Giving the actual detailed proofs would require introducing various definitions, technical facts and auxiliary lemmata about non-commutative $$\mathcal{L}^p $$ spaces, the Cauchy-Stieljtes transform, the $$R $$-transform and analytic subordination. Consequently, we will try to sketch the main ideas with minimal possible reference to the underlying formalism. That being said, it is simply impossible to explain even the ideas without invoking some basic terminology[^2]. Therefore, those readers that are satisfied with the two statements above may freely skip the rest of this subsection without any loss of ability to following the remaining two subsections in this post.

_Idealized Operator Norm_: We locate the maximum of the bulk spectrum of $$\sqrt{2}\beta S - D $$. Let us lay down the strategy to do so:
1. First, we compute the _explicit_ form of a function[^3] which reduces the resolvent[^4] of the sum of two (free) matrices $$X $$ and $$Y $$ to the resolvent of just one matrix $$X $$ evaluated at a shifted argument. This is a standard (and first-half of a remarkable) fact borrowed from [[Bia98]](), and says that:

    $$ \begin{equation} \exists! f: \mathbb{H} \to \mathbb{H}, \text{s.t.},\,\text{i) }f(z) = z + O(1)\text{ for z large enough, and ii) }g_{X+Y}(z) = g_X(f(z))\,.   \end{equation} $$

2. Using the explicit form of $$f^{-1} $$ for our specific case with $$X = -D $$ and $$Y = \sqrt{2}\beta S $$ we will conclude that, under a normalization assumption about the trace of $$D^{-2} $$, the function $$f^{-1} $$ is strictly monotone in a certain region of the complex plane and, therefore, invertible.

3. At this point, we will reason about the analytic behavior of the resolvent of $$\sqrt{2}\beta S-D $$ at the shifted argument. In particular, we will demonstrate that the resolvent is analytic in $$\left(2\beta^2 g_{-D}(0), \infty\right) $$ to obtain an upper bound and then do some approximation to obtain a matching lower bound.


First, we use the fact that the $$R $$-transform is defined as $$r_X(z) = g^{-1}_X(z) - 1/z $$, and that for of a sum of two (asymptotically) free matrices, it simplifies to $$r_{X+Y}(z) = r_X(z) + r_Y(z) $$, to quickly obtain that

$$ \begin{equation} z = g_{\sqrt{2}\beta S - D}\left(g^{-1}_{-D}(z) + r_{\sqrt{2}\beta S}(z)\right)\,, \end{equation} $$

for $$z $$ in a neighborhood around $$0 $$. <br>
Then, combining the above equation with the fact that the $$R $$-transform of a semi-circular operator with variance $$2\beta^2 $$ is $$2\beta^2 z $$, and that $$g^{-1}_{-D}(z) $$ will be in a neighborhood of $$\infty $$, one can deduce that,

$$ \begin{equation} g_{-D}(w) = g_{\sqrt{2}\beta S - D}(w + 2\beta^2 g_{-D}(w))\, ,\end{equation} $$

where we use the fact that $$w = g^{-1}_{-D}(z) $$, for $$w $$ in a neighborhood of $$\infty $$. At this point, since $$w $$ is sufficiently large[^5], one can apply the magical fact of [[Bia98]]() and conclude that, for $$w = f(z) $$,

$$ \begin{equation}  g_{\sqrt{2}\beta S - D}(w + 2\beta^2g_{-D}(w)) = g_{- D}(w) = g_{-D}(f(z)) = g_{\sqrt{2}\beta S - D}(z)\, .\end{equation} $$

Then, using the injectivity of $$g $$ itself, it straightforwardly follows that:
* $$f^{-1}(z) = z + 2\beta^2 g_{-D}(z) $$, and
* $$f(\cdot) $$ is injective.

With some more effort, one can conclude that,

$$ \begin{equation} z \in \mathsf{dom}(f^{-1}) \iff  \mathsf{Im}\left(z + 2\beta^2 g_{-D}(z)\right) > 0\,.\end{equation} $$

Consequently, the domain of $$f^{-1} $$ is a subset of the upper-half (complex) plane where the imaginary part of $$2\beta^2 g_{-D}(z) $$ is "countered" by the imaginary part of $$z $$ itself. It will be possible to strengthen this observation, under the assumption that $$\frac{2\beta^2}{n}\mathsf{Tr}[D^{-2}] = 1 $$, and conclude that  $$\mathsf{sgn}\left(\mathsf{Im}(z)\right) = \mathsf{sgn}\left(\mathsf{Im}(f^{-1})\right) $$. Combining this fact with basic properties about the resolvent allows one to conclude that $$f^{-1} $$ is injective on $$\{z = x + iy \mid x \in (0,\infty), y \in \mathbb{R} \} $$.

_Large-overlap with top-eigenspace_:

#### Approximating the diagonal of the projector
We now formally choose the distortion parameter $$b = \beta n^{-.01} $$. With this choice, we infer that the diagonal entries of $$P(D)^ 2$$ behave desirably, with high probability.

$$ \begin{equation} \text{[Diagonal entries]} \end{equation} $$

### Empirical distribution of the coordinates of the iterates

### Fluctuations of the generalized TAP free energy under fRSB

#### FOOTNOTES

[^1]: Those familiar with resolvent formalism will immediately notice that $$P^2(D) $$ is nothing but the (negative of the) imaginary part of the _resolvent_ $$g_{2\beta A_{sym} - D} $$, evaluated at a purely imaginary point $$ib $$, which is a slight perturbation from $$0 $$. Since, for a compactly supported measure, which we expect the empirical eigenspectrum of $$2\beta A_{sym} - D $$ to converge to, the Cauchy transform (which is equivalent to the resolvent) is invertible everywhere except at poles, this is a loss-less way of studying analytic properties of the spectrum.

[^2]: I will introduce the relevant complex-analytic basics and the relevant free-probability notions in a separate blog-post. In the very same blog-post, I will provide sufficient background so that a motivated reader can actually follow most of the proofs in detail (_except_ for the proof which reasons about the projection of $$P(D)^2 $$ onto the diagonal sub-algebra, which I will introduce in yet another post). For the impatient reader, [[Section 4.1, JSS24](https://arxiv.org/pdf/2408.02360)] should also act as an intense, but reasonably accessible proxy.

[^3]: A function such as this is known in the literature as a **subordination** function. In fact, this function can be extended and shown to be analytic on $$\bar{\mathbb{H}} $$. It allows one to, in some sense, "separate" the support of the spectrum of the two freely independent operators by pushing the domain of one inside the resolvent to be significantly farther than the other.

[^4]: See [here](https://en.wikipedia.org/wiki/Resolvent_formalism) for a definition of the resolvent, and a statement of the famous resolvent identity.

[^5]: It is a classical fact known in free probability theory that a semi-circular operator $$S $$ and an element $$D \in \text{M}_n(\mathbb{C}) $$ with $$\| D\|_{\mathsf{op}} \le O(1)$$ are asymptotically free.

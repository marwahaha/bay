---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model II"
date:   2024-08-12 11:15:40
blurb: "Optimizing the SK model via free probability and convex duality"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)

This is the second blog past in a **3-part** series on a recent result [[JSS24]](https://arxiv.org/abs/2408.02360) by [David](https://davidjekel.com/), [Jonathan](https://jshi.science/) and me. The first blog post can be found [here](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-i).

A quick recap of what we accomplished in the first blog post:
1. We briefly went over the various representations of the Parisi formula, and then derived the Hessian ascent algorithm from the generalized TAP free energy.
2. We then motivated the two main statements we need to prove to demonstrate the success of our algorithm. We concluded by developing a _primal_ theory for the Parisi PDE and AC SDE, which will be quintessential at various points in the main arguments.
<br>

#### Table of Contents
1. [Proof Sketch](#proof-sketch)
   * [Analyzing the TAP-corrected Hessian via free probability](#analyzing-the-tap-corrected-hessian-via-free-probability)
   * [Empirical distribution of the coordinates of the iterates](#empirical-distribution-of-the-coordinates-of-the-iterates)
2. [Footnotes](#footnotes)
<br>
<br>
<br>

## Proof Sketch
### Analyzing the TAP-corrected Hessian via free probability
We now begin sketching the proof overview for the first main statement. As stated in the previous post, the goal is to create a matrix $$Q(t_j, \sigma_j) $$ for every step $$j \in [K] $$ with time $$t_j = j\eta $$ that, with high probability,
1. Projects into the top-eigenspace of $$\nabla^2 \mathsf{TAP} :=2\beta A_{sym} - D(t_j,\sigma_j) $$, and
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

* The scalar $$b $$ will be chosen to be sufficiently small so that it does not cause $$P(D)^2 $$ to deviate too much from projecting into the desired top-eigenspace. However, it will be non-zero all the same, so that the operator $$P(D)^2 $$ continues to be well-defined.

Let us now write a brief and informal statement that summarizes two (of the three) key qualities the operator $$P(D)^2 $$[^1] will have. In the statement below, we will assume that $$\frac{2\beta^2}{n}\mathsf{Tr}[D^{-2}] = 1$$.

$$ \begin{align*} \text{[Idealized operator norm]:}& \quad\quad \mathsf{maxSpec}\left(\sqrt{2}\beta S - D\right) = \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}]\,. \\
\text{[Large-overlap with top-eigenspace]:} & \quad\quad \frac{1}{n}\left\langle P(D)^2,(\widetilde{a(D)}-(2\beta A_{sym} - D))^2\right\rangle \le \frac{O_\beta(\mathsf{Tr}[D^{-4}])}{n^{1.03}}\,. \end{align*} $$

The notation $$O_\beta(\cdot) $$ suppresses the dependencies on $$\beta = O(1/\epsilon)$$ for terms that decay sufficiently fast in $$n $$. For technical reasons, $$\widetilde{a(D)} $$ is not _exactly_ $$\mathsf{maxSpec}\left(\sqrt{2}\beta S - D\right) $$, but a _small_ perturbation around it whereby,

$$ \begin{equation} \left| \widetilde{a(D)} - \frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}] \right| \le \frac{O_\beta\left(\mathsf{Tr}[D^{-3}]\right)}{n^{1.02}}\,. \end{equation} $$

By putting the two statements above together with the observation that $$P(D)^2 $$ projects into the top-eigenspace of $$\sqrt{2}\beta S - D $$, one can reasonably believe that

$$ \begin{equation} \frac{1}{n} \|P(D) -  \Pi_{\delta'}\left(2\beta A_{sym} - D\right)\|_F^2 \le \mathsf{small} \,, \end{equation} $$

where $$\Pi_{\delta'}\left(\sqrt{2}\beta S - D\right) $$ is a projector into the top-$$\delta' n$$ dimensional eigenspace of the spectrum of $$2\beta A_{sym} - D $$.

We will now give overviews of the proofs for the two statements above. Giving the actual detailed proofs would require introducing various definitions, technical facts and auxiliary lemmata about non-commutative $$\mathcal{L}^p $$ spaces, the Cauchy-Stieljtes transform (denoted for the spectral measure of an operator $$A $$ as $$g_A(z) $$), the $$R $$-transform and analytic subordination. Consequently, we will try to sketch the main ideas with minimal possible reference to the underlying formalism. That being said, it is simply impossible to explain even the ideas without invoking some basic terminology[^1]. Therefore, those readers that are satisfied with the two statements above may freely skip the rest of this subsection without any loss of ability to following the remaining two subsections in this post.

_Idealized Operator Norm_: We locate the maximum of the bulk spectrum of $$\sqrt{2}\beta S - D $$. Let us lay down the strategy to do so:
1. Those familiar with resolvent formalism will immediately notice that $$P(D)^2 $$ is nothing but the (negative of the) imaginary part of the _resolvent_[^2] $$g_{2\beta A_{sym} - D} $$, evaluated at a purely imaginary point $$ib $$, which is a slight perturbation[^3] from $$0 $$. First, we compute the _explicit_ form of a function[^4] which reduces the resolvent of the sum of two (free) matrices $$X $$ and $$Y $$ to the resolvent of just one matrix $$X $$ evaluated at a shifted argument. This is a standard (and first-half of a remarkable) fact borrowed from [[Bia98]](), and says that:

    $$ \begin{equation} \exists! f: \mathbb{H} \to \mathbb{H}, \text{s.t.},\,\text{i) }f(z) = z + O(1)\text{ for z large enough, and ii) }g_{X+Y}(z) = g_X(f(z))\,,   \end{equation} $$

where $$\mathbb{H} $$ denotes the upper-half plane.

2. Using the explicit form of $$f^{-1} $$ for our specific case with $$X = -D $$ and $$Y = \sqrt{2}\beta S $$ we will conclude that, under a normalization assumption about the trace of $$D^{-2} $$, the function $$f^{-1} $$ is strictly monotone in a certain region of the complex plane and, therefore, invertible.

3. At this point, we will reason about the analytic behavior of the resolvent of $$\sqrt{2}\beta S-D $$ at the shifted argument. In particular, we will demonstrate that the resolvent is analytic in $$\left(2\beta^2 g_{-D}(0), \infty\right) $$ to obtain an upper bound and then do some approximation to obtain a matching lower bound.


We begin by obtaining an explicit form of $$f^{-1} $$ and inferring basic analytic properties about (its injectivity and domain). <br>
Note that the $$R $$-transform is defined as $$r_X(z) = g^{-1}_X(z) - 1/z $$, and that for of a sum of two free operators, it simplifies to $$r_{X+Y}(z) = r_X(z) + r_Y(z) $$. With this, one quickly obtains that

$$ \begin{equation} z = g_{\sqrt{2}\beta S - D}\left(g^{-1}_{-D}(z) + r_{\sqrt{2}\beta S}(z)\right)\,, \end{equation} $$

for $$z $$ in a neighborhood around $$0 $$. <br>
Combining the above functional equation with the fact that the $$R $$-transform of a semi-circular operator with variance $$2\beta^2 $$ is $$2\beta^2 z $$, and that $$g^{-1}_{-D}(z) $$ will be in a neighborhood of $$\infty $$, one can deduce that,

$$ \begin{equation} g_{-D}(w) = g_{\sqrt{2}\beta S - D}(w + 2\beta^2 g_{-D}(w))\, ,\end{equation} $$

where we define $$w = g^{-1}_{-D}(z) $$. Since $$w $$ is sufficiently large[^5], one can apply the magical fact of [[Bia98]]() and conclude that, for $$w = f(z) $$,

$$ \begin{equation}  g_{\sqrt{2}\beta S - D}(w + 2\beta^2g_{-D}(w)) = g_{- D}(w) = g_{-D}(f(z)) = g_{\sqrt{2}\beta S - D}(z)\, .\end{equation} $$

Then, using the injectivity of $$g_{\sqrt{2}\beta S - D} $$ itself, it straightforwardly follows that:
* $$f^{-1}(z) = z + 2\beta^2 g_{-D}(z) $$, and
* $$f(\cdot) $$ is injective.
* In fact, with some more effort, one can conclude that,

    $$ \begin{equation} z \in \mathsf{dom}(f^{-1}) \iff  \mathsf{Im}\left(z + 2\beta^2 g_{-D}(z)\right) > 0\,.\end{equation} $$

Consequently, the domain of $$f^{-1} $$ is a subset of the upper-half (complex) plane where the imaginary part of $$2\beta^2 g_{-D}(z) $$ is "countered" by the imaginary part of $$z $$ itself. Under the assumption that $$\frac{2\beta^2}{n}\mathsf{Tr}[D^{-2}] = 1 $$, one may strengthen this observation via some simple algebra and conclude that  $$\mathsf{sgn}\left(\mathsf{Im}(z)\right) = \mathsf{sgn}\left(\mathsf{Im}(f^{-1}(z))\right) $$. Combining this fact with basic properties about the resolvent allows one to conclude that $$f^{-1} $$ is injective on $$A = \{z = x + iy \mid x \in (0,\infty), y \in \mathbb{R} \} $$.

At this point, the upper bound follows by using the injectivity of $$f^{-1} $$ in the desired region to conlcude that $$g_{-D}\circ f^{-1} $$ agree on the region of injectivity, and that $$g_{-D}\circ f^{-1} = g_{\sqrt{2}\beta S - D} $$ is analytic on the region $$f^{-1}(A) $$. It is a standard property of the Cauchy-Stieljtes transform that $$f^{-1}(\mathsf{Re}(A)) = \left(\frac{2\beta^2}{n}\mathsf{Tr}[D^{-1}], \infty\right)$$ cannot, then, be in the support of the spectral measure of $$\sqrt{2}\beta S - D $$.

To obtain the lower bound, we do a Taylor expansion of $$f^{-1} $$ near $$0 $$, observing that the first derivative is $$0 $$ and the second derivative is _strictly_ greater than $$0 $$. Using these observations, in conjunction with an approximation argument then demonstrates that the upper and lower values of $$g_{\sqrt{2}\beta S - D} $$ do not agree on $$(f^{-1}(0)-\delta, f^{-1}(0)) $$ for a sufficiently small $$\delta $$. This immediately yields that it is not analytic in the neighborhood. This closes the argument.

_Large-overlap with top-eigenspace_: We now focus on formally showing that the matrix $$P(D)^2 $$ will have large overlap with the top-eigenspace of the shifted Hessian.

We will formally choose the distortion parameter $$b $$, and explicitly compute the overlap between $$P(D)^2 $$ and the subspace spanned by the eigenvectors of $$2\beta A_{sym} - D $$ that are small. Our goal will be to show that there is a valid choice of $$b $$ where this overlap is desirably minimal, which immediately implies that $$P^2(D) $$ must have large overlap with the top-eigenspace. It is easy to see that, by the definition of $$P(D)^2 $$,

$$ \begin{equation} \frac{1}{n}\left\langle P(D)^2, (\widetilde{a(D)}-(2\beta A_{sym} - D))^2\right\rangle = \frac{1}{n}\mathsf{Tr}\left[b(b^2 + (\widetilde{a(D)}-(2\beta A_{sym} - D))^2)^{-1}(\widetilde{a(D)}-(2\beta A_{sym}-D))^{-1}\right]\,. \end{equation} $$

At this point one can use the fact that we will choose $$b > 0 $$ to assert that the above can be bounded as

$$ \begin{equation} \frac{1}{n}\left\langle P(D)^2, (\widetilde{a(D)}-(2\beta A_{sym} - D))^2\right\rangle \le \frac{b}{n}\mathsf{Tr}\left[\mathsf{Id}_n\right] = b\,, \end{equation} $$

and, therefore, it now remains to define and bound $$b $$. Recall from the calculation for the _idealized operator norm_ that $$z \in \mathsf{dom}(f^{-1}) \iff \mathsf{Im}(z + 2\beta^2 g_{-D}(z)) > 0 $$. We will choose $$b = \mathsf{Im}(z + 2\beta^2 g_{-D}(z)) $$ with $$z = ic = i\,\text{small}_n$$. Then,

$$ \begin{equation} \mathsf{Im}(ic + 2\beta^2 g_{-D}(ic)) = c + 2\beta^2\mathsf{Im}(g_{-D}(ic)) = c + \frac{2\beta^2}{n} \left(-b\mathsf{Tr\left[(b^2\mathsf{Id}_n + D^2)^{-1}\right]}\right)\, , \end{equation} $$

which follows with some straightforward algebra and the definition of the resolvent. At this point, adding and subtracting $$D^{-2} $$ inside the trace and invoking the resolvent identity, one obtains,

$$ \begin{equation} \mathsf{Im}(ic + 2\beta^2 g_{-D}(ic)) = c\left(1-\frac{2\beta^2}{n}\mathsf{Tr}\left[D^{-2}\right] + \frac{2\beta^2}{n}\mathsf{Tr}\left[(b^2\mathsf{Id}_n + D^2)^{-1}D^{-2}\right]\right)\,. \end{equation} $$

Using the assumption of the trace of $$\frac{2\beta^2}{n}D^{-2} $$ being normalized to $$1 $$ and the fact that $$b > 0$$ immediately yields,

$$ \begin{equation} b = \mathsf{Im}(ic + 2\beta^2 g_{-D}(ic)) \le \frac{2\beta^2 c^{3}}{n}\mathsf{Tr}\left[D^{-4}\right]\,.\end{equation} $$

A choice of $$c = \beta n^{-.01} $$, along with the fact that our diagonal matrix coming from the entropic correction will be positive-definite and satisfying a bound $$D^{-1} \preceq (1 + O(\epsilon))\mathsf{Id}_n $$ uniformly, will then yield the desired result.

#### Approximating the diagonal of the projector
We have shown that $$P(D)^2 $$ has desirably large overlap with the top-eigenspace of the shifted Hessian, and that the idealized operator norm is $$\approx \frac{2\beta^2}{n}\mathsf{Tr}\left[D^{-1}\right] $$ uniformly over the choices of $$D $$. We now give an overview of how to infer that the diagonal entries of $$P(D)^2 $$ behave desirably. This is essential to arranging the fact that the dynamics of the algorithm converge to the primal Auffinger-Chen SDE; without asserting this convergence, it is unclear how to prove that the energy achieved by the final iterate output by the algorithm is sufficiently large.

We do this via two steps:
* Using the second-half of the remarkable result of [[Bia98]]() to show that the conditional expectation of a free sum can be studied by looking at the algebra of one of the components. Specifically, the result tells us that the same subordination function $$f $$ used above **also** preserves the projection of the free sum $$X + Y $$ onto the diagonal algebra when $$X $$ is an element of the same algebra. Namely,

    $$ \begin{equation} \mathsf{diag}[(z - X + Y)^{-1}] = (f(z) - X)^{-1} \,.\end{equation} $$

    The equation above will allow us to reduce studying $$\mathsf{diag}\left(\tilde{z} - (\sqrt{2}\beta S - D)\right)^{-1} $$ to simply studying $$(f(\tilde{z}) + D)^{-1} $$.

* Comparing the "idealized" free-probability analog $$\mathsf{diag}(\tilde{z} - \sqrt{2}\beta S - D)^{-1} = \left(f(\tilde{z})+ D\right)^{-1} $$ with the "true" random-matrix analog $$\mathsf{diag}\left(\tilde{z} - \left(2\beta A_{sym} - D\right)\right)^{-1} $$. This is done by using concentration of measure for Lipschitz functions of Gaussians to show that,

    $$ \begin{equation} \left\| \mathsf{diag}\left(\tilde{z} - \left(2\beta A_{sym} - D\right)\right)^{-1} - \mathbb{E}\left[\mathsf{diag}\left(\tilde{z} - \left(2\beta A_{sym} - D\right)\right)^{-1}\right] \right\|_2 \le \text{small}_1\,, \end{equation} $$

    and then using a free interpolation, which may be regarded as a non-commutative version of Gaussian interpolation that is ubiquituous in spin-glass theory, to obtain the fact that

    $$ \begin{equation} \left\| \left(f(\tilde{z})+ D\right)^{-1} - \mathbb{E}\left[\mathsf{diag}\left(\tilde{z} - \left(2\beta A_{sym} - D\right)\right)^{-1}\right] \right\|_2 \le \text{small}_2\,. \end{equation} $$

    The free interpolation used to prove the above statement is the most technically involved part of section-4 (and likely the paper). It extends an idea of Collins, Guionnet and Parraud [[CGP22]()] which attempts to understand the quantitative strength of fluctuations between the expected value of $$\frac{1}{n}\mathsf{Tr}\left[f(X_n)\right] $$ and $$\tau\left(f(S)\right) $$; $$X_n $$ is a finite random matrix, $$f $$ is a reasonably "nice" function (not necessarily a polynomial), $$S $$ is an idealized operator in a non-commutative probability space, and $$\tau $$ is the tracial operator in this space. The free interpolation looks at

    $$ \begin{equation} h(t) := \mathbb{E}\left[\frac{1}{n}\mathsf{Tr}\left[(zI_{nk} - 2\beta[\sqrt{1-t}A_{sym} + \sqrt{t}B_{sym}]+ D\otimes I_k)^{-1}(D' \otimes I_k)\right]\right]\,,  \end{equation} $$

    with

    $$ \begin{equation} h(0) = \mathbb{E}\left[\frac{1}{n}\mathsf{Tr}\left[(zI_n -2\beta A_{sym} + D)^{-1}D'\right]\right]\,, \end{equation} $$

    and

    $$ \begin{equation} h(1) = \mathbb{E}\left[\frac{1}{n}\mathsf{Tr}\left[(zI_{nk} - \beta())^{-1}\right]\right]\,. \end{equation} $$

    Here, $$h(1) $$ becomes a semi-circular element as the limit $$k \to \infty $$ is taken at the end of the argument and almost-surely convergence is invoked. The element $$h(0) $$ is the quantity we want to bound, and since we know the limiting spectrum of $$h(1) $$, it remains to show the time-derivative of $$h(t) $$ is sufficiently small and apply the fundamental theorem of calculus. We apply a lemma about the resolvent to expand the derivative $$\frac{d}{dt} h(t) $$; we then evaluate the "easy" term using Gaussian integration-by-parts, and use various Lipschitz estimates on the resolvent $$G(z,t) $$ and the conjugated-resolvent $$h(t) $$ in conjunction with the Poincare inequality to estimate the "hard" term.

The above line of reasoning, with a choice of $$\tilde{z} = f^{-1}(z) = z + 2\beta^2 g_{-D}(z) $$ where $$z = a + ib = \widetilde{a(D)} + ic$$, leads to the conclusion that

$$ \begin{equation} \| \mathsf{diag}\left(\tilde{z} - \left(2\beta A_{sym} - D\right)\right)^{-1} - (z + D)^{-1}\|_2 \le \text{small}_1 + \text{small}_2 \,,\end{equation} $$

at which point taking the imaginary component of the operators immediately yields the desired result,

$$ \begin{equation} \text{[Diagonal entries]} \quad\quad \left\| \mathsf{diag}(P(D)^2) - c(c^2 + D^2)^{-1}\right\|_2 \le \text{small}_1 + \text{small}_2\,.\end{equation} $$

### Empirical distribution of the coordinates of the iterates
With a clear construction of the projection operator $$P(D)^2 $$, and proofs that it projects into the top-eigenspace with desirable diagonal behavior, we are ready to demonstrate that the empirical distribution of the coordinates of the sequence of iterates

$$ \begin{equation} \left(\sigma_1,\dots,\sigma_K\right) \sim \left(\mathcal{N}\left(0, P(D_1)^2\right), \dots,  \mathcal{N}\left(0,P(D_k)^2\right)\right)\, , \end{equation} $$

will approximate, at every step (with high probability), the distribution of the primal Auffinger-Chen SDE $$dY^\gamma_t = \frac{\sqrt{2}\beta}{\partial_{2,2}\Lambda_\gamma(t, Y^\gamma_t)}dW_t $$. We will set an external "clock" where $$t = j\eta $$ for some sufficiently small $$\eta(\epsilon) $$ and $$j \in [K] $$.

To demonstrate that the empirical distribution converges to the desired SDE, we will use an inductive argument and follow a strategy of, once again, using the Lipschitz-ness of the underlying SDE to show concentration of the empirical distribution around its expectation. We will then compare the expected empirical distribution of the coordinates under the inductive hypothesis with a discretization of a step of the primal AC SDE.

To start, we renormalize $$P(D)^2 $$ so that its (normalized) trace is $$1 \pm o_n(1) $$; we also project $$P(D)^2 $$ orthogonal to the current iterate $$\sigma_i $$. This normalization will be helpful in the convergence computations, and the projection will be particularly useful in the energy analysis that follows. The final covariance matrices $$\{Q_i(D)^2\}_{i \in [K]} $$ are then given as,

$$ \begin{equation} Q_i(D)^2 := 2\beta n^{.01}\Pi_{\sigma_i^{\perp}}P(D)^2\Pi_{\sigma_i^{\perp}}\,, \end{equation} $$

where $$D := \left(\frac{2\beta^2}{n}\sum_{j=1}^n \partial_{2,2}\Lambda_\gamma\left(i\eta,(\sigma_i)_j\right)^{-2}\right)^{1/2}\mathsf{diag}\left[\partial_{2,2}\Lambda_\gamma(i\eta,,\sigma_i)\right] $$. It is straightforward to verify with the new definition and using the three conclusions of the prior section that:
1. The normalized trace is,

    $$ \begin{equation} \frac{1}{n}\mathsf{Tr}\left[Q_i(D)^2\right] = 1 \pm O_\beta(n^{-.01})\,.\end{equation} $$

2. The diagonal of $$Q_i(D)^2 $$ is close to $$D^{-2} $$,

    $$ \begin{equation} \left\| \mathsf{diag}\left[Q_i(D)^2\right] - 2\beta^2 D^{-2}\right\|_2 \le  O_\beta\left(n^{.49}\right)\,.\end{equation} $$


In fact, we can also assert that the diagonal of $$Q_i(D)^2 $$ has $$2 $$-norm that is of $$O_\beta(\sqrt{n}) $$ and that the operator norm of the covariance matrix is $$O(n^{.04}) $$. As stated above we will proceed in two steps:
* Consider a _fixed_ iterate $$\sigma_i $$ with desirably small Wasserstein distance to the distribution of $$Y_{i\eta} $$. Then, we first show that,

    $$ \begin{equation} d_{W,2}\left(\mathbb{E}\left[\mathsf{emp}(\sigma_{i+1})\right], Y^\gamma_{(i+1)\eta}\right)^2 \le d_{W,2}\left(\mathsf{emp}(\sigma_{i}), Y^\gamma_{i\eta}\right)^2 + \text{small}_{\beta}\left(\eta,\epsilon\right)\,,\end{equation} $$

    To show this, we first combine a standard estimate using Gronwall's inequality, the Lipschitz-ness of the driving function $$\frac{\sqrt{2}\beta}{\partial_{2,2}\Lambda_\gamma(t,y)} $$, and an application of Ito-isometry to conclude that,

    $$ \begin{equation} \left\|Y^\gamma_{(i+1)\eta} - \left(Y^\gamma_{i\eta} + \frac{\sqrt{2}\beta}{\partial_{2,2}\Lambda_\gamma(t,Y^\gamma_t)\vert_{t=i\eta}}\left(W_{(i+1)\eta} - W_{i\eta}\right)\right)\right\|^2_2 \le O_\beta(\eta^2)\,,\end{equation} $$

    allowing us to approximate the primal AC SDE in small, discretized time increments with sufficiently small quantitative error. At this point, the inductive hypothesis kicks in to accumulate the errors from the previous step, provided we can bound the action of (Lipschitz) test functions on the iterate $$\sigma_{i+1} = \sigma_i + \sqrt{\eta}Q_i(D)z_i $$ where $$z_i \sim \mathcal{N}(0,\mathsf{Id}_n) $$.

    They key point here is to notice that $$\sqrt{\eta}(Q_i(D)z_i)_j $$ can be replaced by $$\left(Q^2_i(D)\right)_{j,j}^{1/2} W$$ with $$W \sim \mathcal{N}(0,\eta) $$ while retaining equality in distribution. This allows us to study the normalized counting measure on the coordinates $$[n] $$. Letting $$\Sigma $$ denote the random variable that takes value $$(\sigma_i)_j $$ on $$j \in [n] $$, and $$\Pi $$ be one that takes value $$\left(Q^2_i(D)\right)_{j,j}^{1/2} W$$ on $$j \in [n] $$, a simple rewrite then implies that,

    $$ \begin{equation} \mathbb{E}[\mathsf{emp}(\sigma_{i+1})] \overset{d}{=} \Sigma + \Pi W\, , \end{equation} $$

    and consequently, optimally coupling $$\Sigma $$ and $$Y^{\gamma}_{i\eta} $$,

    $$ \begin{equation} \|\Sigma -  Y^\gamma_{i\eta}\|_2 = d_{W,2}\left(\mathsf{emp}(\sigma_i), Y^\gamma_i\right)\,. \end{equation} $$

    At this point, we use the definition of the $$W_2$$ distance in conjunction with a triangle inequality to conclude

    $$ \begin{align*} d_{W,2}\left(\mathbb{E}\left[\mathsf{emp}(\sigma_{i+1})\right], Y^\gamma_{(i+1)\eta}\right) &\le \|(\Sigma + \Pi W)- Y^\gamma_{(i+1)\eta}\|_2 \\
    &\le \| \Sigma - Y^\gamma_{i\eta}\|_2 + \| \Pi W - (Y^\gamma_{(i+1)\eta} - Y^\gamma_{i\eta})\|_2\,. \end{align*} $$   

     Now, we apply the inductive hypothesis, and use the discretized primal AC SDE estimate above in conjunction with the properties of $$Q_i(D)^2 $$ mentioned above to conclude that

     $$ \begin{equation} d_{W,2}\left(\mathbb{E}\left[\mathsf{emp}(\sigma_{i+1})\right], Y^\gamma_{(i+1)\eta}\right) \le \text{small}_{\beta}(\eta,\gamma) + o_n(1)\,. \end{equation} $$

* Having shown this, we use concentration of measure arguments to argue that the empirical distribution doesn't fluctuate a lot. Namely,

    $$ \begin{equation} d_{W,2}\left(\mathsf{emp}\left(\sigma_{i+1}\right), \mathbb{E}\left[\mathsf{emp}\left(\sigma_{i+1}\right)\right]\right) \le o_n(1)\,.\end{equation} $$

    To do this we use the fact that we can upper bound, with high probability over the randomness of the iterate, the maximum value of every coordinate. Then, since the operator norm of the covariance is bounded, we can apply Lipschitz test functions on a single updated iterate and project the underlying space into a large enough cube which contains the current iterate. On doing this, one can use concentration for Lipschitz functions to bound the behavior of $$1 $$-Lipschitz functions, and make the strength of the concentration fight the size of a net over all Lipschitz functions.

#### FOOTNOTES

[^1]: I will introduce the relevant complex-analytic basics and the relevant free-probability notions in a separate blog-post. In the very same blog-post, I will provide sufficient background so that a motivated reader can actually follow most of the proofs in detail (_except_ for the proof which reasons about the projection of $$P(D)^2 $$ onto the diagonal sub-algebra, which I will introduce in yet another post). For the impatient reader, [[Section 4.1, JSS24](https://arxiv.org/pdf/2408.02360)] should also act as an intense, but reasonably accessible proxy.

[^2]: See [here](https://en.wikipedia.org/wiki/Resolvent_formalism) for a definition of the resolvent, and a statement of the famous resolvent identity.

[^3]: Since, for a compactly supported measure, which we expect the empirical eigenspectrum of $$2\beta A_{sym} - D $$ to converge to, the Cauchy transform (which is equivalent to the resolvent) is invertible everywhere except at poles, this is a loss-less way of studying analytic properties of the spectrum.

[^4]: A function such as this is known in the literature as a **subordination** function. In fact, this function can be extended and shown to be analytic on $$\bar{\mathbb{H}} $$. It allows one to, in some sense, "separate" the support of the spectrum of the two freely independent operators by pushing the domain of one inside the resolvent to be significantly farther than the other.

[^5]: We work in the free product of the algebras generated by a semi-circular operator $$S $$ and an element $$D \in \text{M}_n(\mathbb{C}) $$ with $$0 < \| D\|_{\mathsf{op}} \le O(1)$$.

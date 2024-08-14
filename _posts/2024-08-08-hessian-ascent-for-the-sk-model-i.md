---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model I"
date:   2024-08-08 11:15:40
blurb: "Optimizing the SK model via free probability and convex duality"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)


Recently, [David](https://davidjekel.com/), [Jonathan](https://www.jshi.science/) and I put out a preprint on the arXiv [(Potential Hessian Ascent: The Sherrington-Kirkpatrick Model, JSS24)](https://arxiv.org/abs/2408.02360).

This work might seem somewhat opaque and technical to the broader TCS (and even probabilist) audience, so I have decided to write a **3-part** blog post explaining the background, motivation, proof structure and main technical innovations of this result, as well as its relationships to other works in the literature and the open questions it naturally induces.

We begin by introducing relevant past work: the Sherrington-Kirkpatrick model, the Parisi formula, and alternative representations of the latter via the Auffinger-Chen stochastic differential equation (SDE) and the generalized TAP free energy. We then sketch our development of the _primal_ theory for the Parisi PDE (and Auffinger-Chen SDE) which is central to the proofs of the algorithm’s key theorems. In later posts we will show how the Hessian ascent algorithm maximizes the generalized TAP free energy by giving a proof sketch for these key statements [(Part-2)](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-ii), and discuss connections to other work in the literature, ongoing work and open questions [(Part-3)](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-iii).  

**<u>The paper in a nutshell</u>**: [[JSS24]](https://arxiv.org/abs/2408.02360) introduces and analyzes a Hessian ascent algorithm for the SK model, the update rules of which are motivated (in part) to resolve a conjecture of Eliran Subag [[Sub18, Pg. 8, Ising Spins]](https://arxiv.org/abs/1812.04588). Subag gave an (essentially) equivalent algorithm for the same models on the sphere. Due to various technical reasons stemming from geometry, the analysis and conceptual understanding for the Hessian ascent algorithm on the cube is significantly more demanding than it is on the sphere.

The model is given by the following cost function,

$$ \begin{equation} H(\sigma) = \langle \sigma, A\sigma\rangle\,,\end{equation}$$

where $$A $$ is a random matrix with i.i.d. $$\mathcal{N}(0,1/n) $$ entries. It can be shown by standard concentration arguments for Lipschitz functions of Gaussians that the maximum value of this cost function is concentrated around $$\mathbb{E}[\max_{\sigma \in \{-1,1\}^n}H(\sigma)] $$.

At a high-level, the algorithm is fairly simple:
* Set the starting point $$\sigma_0 = (0,\dots,0), $$ and step-size $$\eta = \text{small} $$.
* For $$i \in [K = O(1/\eta)] $$:
    * Set $$Q_i = \text{smooth projector on to top-eigenspace of TAP-corrected Hessian orthogonal to }\sigma_{i-1} $$.
    * Set $$\sigma_i = \sigma_{i-1} + \sqrt{\eta}z $$, where $$z \sim \mathcal{N}(0,Q^2_i) $$.
* Return $$\sigma_K $$ after truncating it and rounding it to $$\{-1,1\}^n $$.

The details of what the "TAP corrected Hessian" of the SK model are will be introduced later, and do not matter for now. The main points are the simplicity of the update rules, and the fact that we are following the top-eigenspace of a certain Hessian matrix induced by the SK model. With this high-level description in hand, we now briefly explain what the main motivations for rigorously developing and analyzing this algorithm are.

**Primary motivation**: <u>Give a conceptually simple and purely spectral algorithm for these models on the cube</u>.

  Ideally, the algorithm should implement the same principle on the cube that Subag's algorithm implements on the sphere, leading to a unified principle to optimize (to the conjectured limit for $$\mathsf{poly}(n) $$-time algorithms) random polynomials over _any_ "reasonable" underlying domain $$D $$.

  The Hessian ascent algorithm designed in this paper accomplishes this on the cube, in the same way it was done on the sphere by Subag, and **resolves Subag's conjecture affirmatively**. The universal principle is to:
  > Maximize the objective plus the generalized TAP correction term in small, orthogonal increments starting from the center and eventually ending up at the boundary.

  As it turns out, this generalized TAP correction is a precise measure of entropy, and its derivatives are understood by a particular partial differential equation (which is closely related to the famous Parisi PDE). Deriving and understanding the properties of this PDE, which has an initial condition given by the entropy of the Bernoulli-$$1/2 $$ random variable, is critical in a part of the energy analysis for the algorithm.

There are also a few **secondary motivations**:
1. Understand what properties about random polynomials can be efficiently certified using a sum-of-squares hierarchy that Jonathan and I introduced in a previous work [[SS24]](https://arxiv.org/abs/2401.14383).
2. Ultimately, give a principled, mathematically rigorous and domain-generalized derivation of the Parisi formula (in certain regimes).

The two points above have connections to different areas of the literature, and are motivated by a desire to understand algorithms such as Hessian ascent, AMP and low-degree polynomials as part of a larger family of algorithms that can (potentially) be unified and generalized in various ways -- more on this later in [Part-3](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-iii).
<br>

#### Table of Contents
1. [The Sherrington-Kirkpatrick Model](#the-sherrington-kirkpatrick-model)
   * [The Parisi formula and Auffinger-Chen Representation](#the-parisi-formula-and-auffinger-chen-representation)
   * [The generalized TAP free energy](#the-generalized-tap-free-energy)
   * [A primal theory for the Parisi PDE via convex duality](#a-primal-theory-for-the-parisi-pde-via-convex-duality)
2. [Footnotes](#footnotes)
<br>

## The Sherrington-Kirkpatrick Model
We briefly introduce the Sherrington-Kirkpatrick model as an optimization problem. We describe the expected limit of the optimal value, as a standard Gaussian concentration inequality implies that the value is extremely concentrated around its expectation. This limit exists almost-surely, and is captured by the famous Parisi formula at zero temperature. This formula has a long history in the statistical physics and probability theory literature (see [Bolthausen's overview of the proof of the Parisi formula](http://archive.numdam.org/article/SB_2004-2005__47__349_0.pdf)) that we will not spend much time on in this post. For us, the critical facts are that:
* The Parisi formula is a variational formula over certain measures, and that it is strictly convex over this space having a unique minimizer. Furthermore, one can efficiently find this minimizer since the variational formula is independent of $$n $$ (though, still infinite-dimensional, and requiring a clever technique of [[JT16]](https://arxiv.org/abs/1502.04398)).
* The Parisi formula involves two terms, one of which is a solution to the Parisi PDE. This solution can be rewritten in terms of an optimal stochastic control problem, and this rewrite is called the Auffinger-Chen representation.
* There is (yet) another alternative representation for the Parisi formula which _extends_ into the interior of the solution domain. This quantity, developed rigorously and systematically by Subag [[Sub18]](https://arxiv.org/abs/1804.10576) for the sphere, and then again by Chen, Panchenko and Subag [[CPS18]](https://arxiv.org/abs/1812.05066v2) for the cube, is the heart of the algorithmic design principle mentioned above. In fact, Subag is able to use the ideas related to this representation to give a derivation of the Parisi formula on the sphere (known as the Crisanti-Sommers formula) from _first principles_. Unfortunately, this is (seemingly) not quite yet in reach for the cube. Nonetheless, there is remarkable conceptual similarity between the generalized TAP equation on the sphere and the cube in the way of its component terms, and this is enough to demonstrate how Subag's conjecture about step-wise optimization on the cube is resolved via our approach.

We begin by defining the Parisi formula and the Auffinger-Chen representation, and mention a few reasons why these representations are insufficient to analyze our algorithm's performance. Then, we introduce the generalized TAP free energy and show how, by following Subag's principle, we can derive the correct form of our algorithm (a Hessian ascent) through a simple line of reasoning. This also naturally leads to the framing of the two main statements we wish to prove about our algorithm. We conclude by developing alternative representations of the Parisi PDE and Auffinger-Chen SDE in _primal_ space and analyzing various properties about these, as they form the backbone of our algorithm's analysis.   
<br>


### The Parisi formula and Auffinger-Chen Representation
The Parisi formula was proposed by Giorgio Parisi in 1979 to give the asymptotic free energy density of the Sherrington-Kirkpatrick (SK) model. It was first rigorously proved in a series of works by Guerra and Talagrand, the culmination of which was in the hard direction of the proof by Talagrand.  

The expected value of the maximum of the SK model hamiltonian is given by the zero-temperature limit of the famous Parisi formula. The Parisi formula for the SK model at inverse temperature $$\beta = 1/T $$ is given by,

$$ \begin{equation} P_\beta(\mu) = \Phi_\mu(0,0) - \beta^2\int_0^1t\mu(t)dt\,,\end{equation} $$

where $$\mu(.) $$ is a cumulative-distribution function for a probability measure with support in $$[0,1] $$ and $$\Phi_\mu $$ is the solution to the following parabolic PDE,

$$ \begin{equation} \partial_t \Phi_\mu(t,x) + \beta^2\left(\partial_{xx} \Phi(t,x) + \mu(t)(\partial_x \Phi(t,x))^2\right)= 0\,,\end{equation}$$

with initial condition $$\Phi(1,x) = \log(2\cosh(x)) $$. Various properties about the solutions of this equation, and a zero temperature variant of it, are known. The main two points for a reader of this post about this PDE are:
1. Its solution exists, is well-posed, and has various pleasant regularity properties: namely, the spatial derivatives are Lipschitz, and the mixed derivatives have regularity depending on the properties of $$\mu $$.
2. $$\Phi $$ is actually the Fenchel-Legendre (FL) dual (in $$x $$) of an entropic function, and the initial condition $$\Phi(1,x) $$ in particular is the FL dual of the Bernoulli-$$1/2 $$ entropy.

#### The primal picture

This last point is important: a large part of the previous work has focused on the "dual" picture with $$\Phi $$ and its argument $$x \in (-\infty, \infty) $$. We find it more enlightening to work in the "primal picture" with $$\Lambda(t,y) $$ as the FL dual to $$\Phi(t,x) $$. Then the argument $$y \in [-1,+1] $$ of $$\Lambda $$ can be interpreted as the mean of a binary distribution over the two possible output values $$\{-1,+1\} $$ for each coordinate $$\sigma_i $$, and when $$t $$ approaches 1 near the end of the algorithm, $$-\Lambda $$ itself becomes the Shannon entropy of that binary distribution.

We can then apply this coordinate-wise to $$\sigma \in [-1,+1]^n $$ to obtain a product distribution over $$\{-1,+1\} $$ whose mean is $$\sigma $$ and whose entropy is $$-\sum_i\Lambda(t, \sigma_i) $$ when $$t $$ is near 1.
Since $$H(\sigma) $$ is (ignoring the negligible diagonal entries of $$A $$) equal to the expected value of $$H(\sigma^*) $$ when $$\sigma^* $$ is sampled from that product distribution, the Gibbs variational principle tells us that the free energy of this distribution is $$\beta H(\sigma) - \sum_{i}\Lambda(t, \sigma_i) $$.

Thus, if we maximize the objective function $$\beta H(\sigma) - \sum_{i}\Lambda(t, \sigma_i) $$, we are simply maximizing the free energy of the product distribution centered on $$\sigma $$.

The differential equation defining $$\Phi(t,x) $$ in the Parisi formula can be recast as a differential equation defining $$\Lambda(t,y) $$ for $$t < 1 $$ (and this is discussed further in [(1.3)](#a-primal-theory-for-the-parisi-pde-via-convex-duality)), and so we can interpret this PDE as characterizing a "modified entropy" which tells us about the exact explore-exploit trade-off needed at different stages of the algorithm in order to maximize the free energy at the end.

This treatment of $$\Lambda$$ was initiated by Chen, Panchenko, and Subag in work on the generalized TAP representation [[CPS18]](https://arxiv.org/pdf/1812.05066v3). We expand on these ideas by developing the PDE and SDE theory for the primal representation, and for certian technical reasons, we must work with a "smoothened" version of $$\Lambda $$, written as $$\Lambda_\gamma $$: the development of this primal theory is discussed further in [(1.3)](#a-primal-theory-for-the-parisi-pde-via-convex-duality).

There is some recursive dependence between how $$\Lambda$$ affects the evolution of $$\sigma$$ in the algorithm and how the future evolution of $$\sigma$$ affects the desired explore-exploit trade-off that $$\Lambda$$ should encode. To understand (and explicitly compute quantities related to) how $$\sigma$$ and $$\sum_i\Lambda(t,\sigma_i)$$ co-evolve in the algorithm, we need to utilize a reformulation of $$\Lambda $$ (and $$\Lambda_\gamma $$) as an _optimal stochastic control_ problem. Using the SDE underlying this representation, the computation of the quantities such as the value of $$\Lambda $$ at the final point of our (randomized) algorithm become tractable.

For $$\Phi $$ itself, this reformulation as a _stochastic optimal control_ problem is known as the **Auffinger-Chen representation** [[AC15]](). It gives a rewrite for the Parisi PDE. More specifically, if $$\Phi $$ is a solution to the Parisi PDE, then,

$$ \begin{equation} \Phi(0,x) = \mathbb{E}\left[\Phi\left(1, x + X_1\right) - \beta^2\int_0^1\mu(s)(\partial_x \Phi(s,X_s))^2 ds\right]\, ,\end{equation} $$

where $$X_s $$ is given by the strong solution to the following Ito drift-diffusion process,

$$ \begin{equation} dX_s = 2\beta^2 \mu(t)\partial_x \Phi(s,X_s)ds + \sqrt{2}\beta dW_s\,, \,\,\, X_0 = 0\,. \end{equation}$$

As hinted at above, we need to write down a new SDE corresponding to the Auffinger-Chen (AC) SDE above that runs in the primal space, where the iterates of our Hessian ascent algorithm reside. This is born directly out of the fact that the quantities in the generalized TAP correction include terms dependent on $$\Lambda $$. With the primal SDE (introduced with the primal Parisi-like PDE in [(1.3)](#a-primal-theory-for-the-parisi-pde-via-convex-duality)) we can demonstrate that the coordinates of every iterate of the Hessian ascent algorithm (whp) converge, in Wasserstein distance, to the primal version of the AC SDE. Then, we can use Ito calculus to (approximately) compute various important quantities under the dynamics of the algorithm.

This is perhaps one reason why the approximate-message passing (AMP) algorithms previously used to optimize these problems are somewhat opqaue to a particularly neat conceptual interpretation: the iterates $$\{u_k\}_{k} $$ of the AMP algorithm [[Mon19]](https://arxiv.org/abs/1812.10897) are updated by discretizing the AC SDE, taking a single step of the discretized SDE, then updating the current iterate by multiplying the Hessian of the objective with a non-linearity evaluated at this new step rescaled by the current iterate value coordinatewise, and finally correcting the whole computation by subtracting a term known as the _Onsager correction_. Upon running these complicated updates for a sufficient number of iterations, these iterates are then again rescaled (by a small parameter $$\sqrt{\delta} $$) and transformed (yet) again (non-linearly) into the _primal_ space, where the iterates are truncated and rounded onto the hypercube. <br>
Au contraire, as we saw, the Hessian ascent algorithm remains _squarely_ in the primal space, leading to a very simple update rule as well as a very clean conceptual interpretation a-la Subag's principle mentioned in the introduction.
<br>


### The generalized TAP free energy
In the previous subsection, we introduced the Parisi formula, the Auffinger-Chen representation, and mentioned how they work in a dual space. We then stated that, given the fact that the TAP correction involves a FL dual to the solution of the Parisi PDE, for various reasons in the analysis of the Hessian ascent algorithm, we will need to develop a primal version of the Parisi PDE and AC SDE.

We now table the development of this _primal_ PDE and SDE to [(1.3)](a-primal-theory-for-the-parisi-pde-via-convex-duality), and first write down the form of the generalized TAP free energy on the cube introduced in [[CPS18]](https://arxiv.org/abs/1812.05066v2). We will briefly interpret the generalized TAP correction and write down its gradient at a critical point, yielding the generalized TAP equation. Then, observe that jumping from one critical point (where the gradient is $$0 $$) to the next entails moving along the _kernel_ of the _Hessian_ of the generalized TAP free equation. It will turn out that, in order to move along a path of critical points of the generalized TAP equation by moving into this kernel, one will automatically be forced to climb the top-eigenspace of,

$$ \begin{equation} \nabla^2\bigg(\langle \sigma, A \sigma\rangle - \text{FL dual to }\Phi(t,x)\bigg) = A - \nabla^2\left(\text{FL dual to }\Phi(t,x)\right)\, .\end{equation} $$

For largely historical reasons[^1] we will use a _convex_ dual, while [[CPS18]](https://arxiv.org/abs/1812.05066v2) use a _concave_ dual to define the generalized TAP free energy correction term. Let us begin by defining the FL dual to the solution $$\Phi $$ of the Parisi PDE,

$$\begin{equation} \Lambda(t,y) = \sup_{x \in \mathbb{R}}\left(xy - \Phi(t,x)\right) \end{equation}\,, $$

for every $$t \in [0, 1] $$. It is well known that $$\Phi $$ is strictly convex in $$x $$ [[Section 2, CPS18]](https://arxiv.org/abs/1812.05066v2), and using this (along with the regularity properties of $$\Phi(t, x) $$ mentioned prior) it is not hard to ascertain the following facts:
* The maximizer of the FL dual is unique for every $$y \in [-1,1] $$ and obeys the relationship $$y = \partial_x \Phi(t,x) $$.
* In fact, the maps $$y = \partial_x \Phi(t,x) $$ and $$x = \partial_y \Lambda(t,y) $$ define a change of coordinates from $$\mathbb{R} \to \{-1,1\} $$ and vice-versa. There are two points to notice here:
    * The derivatives of $$\Phi $$ and $$\Lambda $$ are maps that take us between (convexly) dual spaces. While the generalized TAP energy, involving $$\Lambda $$ stays in the primal space, the standard Parisi machinery is in the dual space. This implies that, beneath the entire framework, there is some mysterious convex geometry at work.
    * They are inverse functions of each other, except for some minor issues that arise on the corners where $$\partial_{y}\Lambda = \partial_x\Phi^{-1} $$ blows up. To overcome this, a regularization to the FL dual is introduced (see [[Section 2, JSS24]](https://arxiv.org/abs/2408.02360) and [1.3]()). In fact, the algorithm actually follows this _regularized_ FL dual that has desirable smoothness properties around the corners of the cube. As we shall see in [1.3](), this regualized FL dual stays uniformly close to $$\Lambda $$ inside the cube, and is negligible outside, meaning the errors introduced are manageable.

The use of the convex dual actually ends up being convenient, because at a critical point $$\sigma $$,

$$ \begin{equation} \frac{1}{n}\nabla H(\sigma) = -\nabla\mathsf{TAP}(\sigma)\,, \end{equation} $$

and our choice of FL dual gives the following form for the TAP correction term,

$$ \begin{equation} \mathsf{TAP}_{\text{emp}}(\sigma) = -\int \Lambda(t,\sigma) d(\text{emp}(\sigma)) - \beta^2\int_{\frac{1}{n}\|\sigma\|^2_2}^1 t\mu(t)dt\, , \end{equation} $$

analogous to what is given in [[Eq. 1.27, CPS18]](https://arxiv.org/abs/1812.05066v2). The generalized TAP equation for mean-field spin glasses on the cube is introduced in [[Theorem 2, CPS18]](https://arxiv.org/abs/1812.05066v2). Using the explicit representation of the generalized TAP equation along with the facts about the FL duals above, some calculus yields that the critical points $$\sigma $$ satisfy the following eigenvector equation,

$$ \begin{equation}  \left(\beta A - \left(2\beta^2\int_q^1 \mu(t)dt\right)\mathsf{Id}\right)\sigma = \left(\partial_{\sigma_1}\Lambda(q,\sigma_1),\dots,\partial_{\sigma_n}\Lambda(q,\sigma_n)\right)\,, \end{equation} $$

where $$q = \frac{1}{n}\|\sigma\|^2_2 $$. The representation above is equivalent to [[Remark 6, CPS18]](https://arxiv.org/abs/1812.05066v2). Now comes the crucial part: We would like to have an algorithm that follows small (orthogonal) updates, such that, _every_ point is a critical point along the path. This means that we _must_ actually proceed in a direction where the Hessian of the TAP equation (projected orthogonal to the current location) is zero. Equivalently, we must stay in the kernel of the Hessian projected orthogonal to the current iterate.

After taking the gradient of the above equation, applying chain rule, using the FL duality rewrites, and discarding terms that are rank-$$1 $$ or along the current iterate ($$\sigma $$), we arrive at the fact that the update $$\Delta \sigma $$ must satisfy the following condition vis-a-vis the Hessian of the TAP equation,

$$ \begin{equation} \Delta\sigma\left(2\beta A_{\text{sym}} -\underbrace{\sum_{i} \partial_{\sigma_i\sigma_i}\Lambda(q,\sigma_i) e_ie_i^{\mathsf{T}}}_{D'(t,\sigma)} - \frac{2\beta^2}{n}\sum_{i \in [n]}\partial_{x_ix_i}\Phi(q,x_i)\mathsf{Id}\right)\Delta\sigma \approx 0\end{equation}\,, $$

where $$A_{\text{sym}} = (A + A^{\mathsf{T}})/2 $$ is distributed as $$\sqrt{2}\,\mathsf{GOE}(n) $$. Therefore, _if_ we are at a critical point $$\sigma $$ and we wish to make a small $$\approx \eta $$ sized increment that jumps to the next TAP state (critical point), it is non-negotiable that the quadratic form with the matrix in the Hessian term above be (approximately) $$0 $$. This basically implies that we want to take (small) steps in the eigenspace of

$$ \begin{equation} 2\beta A_{\text{sym}} - D'(q,\sigma) \end{equation} $$

with value

$$ \begin{equation} \approx \frac{2\beta^2}{n}\sum_i \partial_{x_ix_i}\Phi(q,x_i) = \frac{2\beta^2}{n}\mathsf{Tr}[D'^{-1}(q,\sigma)]\, , \end{equation}$$

where we use the fact that whenever $$\partial_{xx}\Phi(t,x) > 0 $$, its reciprocal is well-defined and equal to $$\partial_{yy} \Lambda(t,y) $$ when $$x $$ and $$y $$ satisfy the change of coordinates implied by FL duality (that is, they are critical points in their respective bases). This identity is called the [Crouzeix identity in convex analysis](https://link.springer.com/article/10.1007/BF01584350), and is an important observation in working out the details of the primal Parisi theory ([1.3]())  as well as understanding the conceptual basis on which the free-probabilistic analysis of the TAP-corrected Hessian proceeds.

As it turns out, the desired value will be achieved in the _top-eigenspace_ of the TAP corrected Hessian (see [(2.1)](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-ii)) and, therefore, we will need an iterative argument where we can construct a covariance matrix $$Q^2(\sigma) $$, such that $$Q(\sigma) $$ smoothly projects into the top eigenspace of,

$$ \begin{equation} \text{TAP-corrected Hessian} = 2\beta A_{\text{sym}} - D'(q,\sigma) \overset{d}{=} \sqrt{2}\beta\,\mathsf{GOE}(n) - D'(q,\sigma)\,. \end{equation} $$

We _will_ be able to accomplish this, _conditioned_ on the fact that,

$$ \begin{equation} \frac{2\beta^2}{n}\mathsf{Tr}[D'^{-2}(q,\sigma)] = 1\,, \end{equation} $$

which will in-turn require a re-normalization of $$D'(t,\sigma) $$. To have this be the case, we will choose the final diagonal correction matrix as,

$$ \begin{equation} D(t,\sigma) = \left(\frac{2\beta^2}{n}\sum_{i \in [n]}\partial_{\sigma_i,\sigma_i}\Lambda(t,\sigma_i)^{-2}\right)^{1/2}D'(t,\sigma)\,. \, \end{equation} $$

It will be shown in [(2.1)](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-ii) that this is the right scaling to obtain the desired eigenvalue in the top-eigenspace of $$2\beta A_{\text{sym}} - D'(t,\sigma) $$. Having achieved an understanding of what space we want to project to at every step, the remaining critical task is to arrange the (efficient) construction of a matrix $$Q(\sigma) $$ that projects into this top eigenspace _and_ has diagonal entries that behave roughly like $$D'(t,\sigma)^{-2} $$. You may (correctly) inquire:
> Why must the diagonal entries of square-root of the covariance matrix behave akin to $$D(t,\sigma)^{-2} $$?

The answer is that without this arrangement, we will not be able to demonstrate that the empirical distribution of the coordinates of

$$ \begin{equation} \sigma_{k+1} = \sigma_k + \eta^{1/2}Q(\sigma_k)w\,,\,\,\, w \sim \mathcal{N}(0,\mathsf{Id})\, , \end{equation} $$

converges to the primal version of the AC SDE that we desire, where $$\sigma_k $$ is a critical point (and having empirical coordinate distribution sufficiently close to the primal AC SDE itself). Without this, various points in the analysis remain intractable (including the two points mentioned in the [(1.1)]()).

Consequently, in [(2.2)](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-ii) we will show that given that the update comes from an (appropriately rescaled) eigenvector in the top-eigenspace of the TAP-corrected Hessian, the empirical distribution of the coordinates of its iterates will converge (in Wasserstein-$$2 $$ distance) to the _primal_ version of the AC SDE (with high probability).

At this point, it then remains for us to demonstrate the following inductive steps:
1. Under the condition that $$(2\beta^2)/n \mathsf{Tr}[D^{-2}(t,\sigma)] = 1 $$, we can (efficiently) compute a matrix $$Q(\sigma) $$ which projects into the top-eigenspace of $$2\beta\,A_{\text{sym}} - D(t,\sigma) $$ with diagonals that are $$\approx D(t,\sigma)^{-2} $$ (see [2.1]()) to choose the update $$\Delta\sigma := \sqrt{\eta}\,Q(\sigma)w$$ with $$w \sim \mathcal{N}(0, \mathsf{Id}_n) $$, and
2. Provided that $$\Delta\sigma $$ is chosen as above at every step, it is the case that $$\mathsf{Wass}_2(\mathsf{emp}(\sigma_j), Y_j) \le \mathsf{small}_j $$ for every $$j \in [K] $$, where $$Y_j $$ is generated by the primal AC SDE process.
<br>


### A primal theory for the Parisi PDE via convex duality

Having built up to the the two main properties that we will need to prove to conduct a successful analysis of the algorithm, we now focus on building the analytic tools that will be used (time and again) in the proofs of these two properties. These tools are:
1. The introduction of a $$\gamma $$-_regularized_ (or $$\gamma $$-smoothened) FL dual to $$\Phi $$, termed $$\Lambda_\gamma $$, along with its regularity properties. In particular, we would like that for (sufficiently) large $$\beta $$, the regularized FL dual is a uniformly good approximate for $$\Lambda $$ with "sane" derivatives, especially near the corners of the hypercube.
2. The introduction of the primal version of the Parisi PDE and the AC SDE, written for $$\Lambda_\gamma $$ and $$\Lambda $$, along with Wasserstein distance bounds between these.

We begin by stating the definition of $$\Lambda_\gamma $$,

$$ \begin{equation} \Lambda_\gamma(t,y) = \sup_{x \in \mathbb{R}}\left(xy - \Phi_\gamma(t,x)\right) = \sup_{x \in \mathbb{R}}\left(xy - \Phi(t,x) - \frac{\gamma}{2}x^2\right)\,. \end{equation} $$

By similar convex analytic reasons as the ones hinted at in [(1.2)](#the-generalized-tap-free-energy), we have that,

$$ \begin{equation} \partial_y\Lambda_\gamma = \partial_x\Phi_\gamma^{-1} = \left(\partial_x \Phi + \gamma x\right)^{-1}\,, \end{equation} $$

and

$$ \begin{equation} y = \partial_x \Phi(t,x) + \gamma x\,. \end{equation} $$

Unfortunately, for $$\Lambda $$ itself, one obtains that $$\partial_y \Lambda = (\partial_x \Phi)^{-1} $$ which goes to $$\infty $$ when $$y > 1 $$ or $$y < -1 $$. This is the main reason for introducing the regularization.

_<u>Estimates for </u>$$\Lambda $$_: We now focus on continuity estimates for $$\Lambda $$, which will be especially important in estimating how well $$\Lambda_\gamma $$ approximates the former in the solid cube (uniformly). To obtain these, we use a stochastic expression for $$\partial_x \Phi $$ (see [[Lemma 2.3, JSS24](https://arxiv.org/pdf/2408.02360)]) as an average over a function of the process $$X_t $$ that solves the AC SDE. This allows us to obtain regularity estimates for $$\partial_x \Phi(t,x) $$ and $$\partial_{xx} \Phi(t,x) $$ slightly more refined than those written down in the literature (see [[Proposition 2, AC15]](https://arxiv.org/pdf/1402.5132), [[JT16]](https://arxiv.org/abs/1502.04398), [[Chapter-14.7, Tal11]](https://link.springer.com/book/10.1007/978-3-642-22253-5)). For accomplishing this, the idea is simple:
> Using the stochastic expression for $$\partial_x \Phi(t,x) = \mathbb{E}[\tanh(X_1)]$$ provided by [JT'16], wield Ito calculus with an application of Gronwall's inequality to bound the MGF of $$X_t $$. Then, use bounds for hyperbolic functions to sharpen the estimates from the literature.

Using the coordinate change of $$y = \partial_x \Phi $$ and $$x = \partial_y \Lambda $$ one can transfer these bounds to the primal space and obtain Lipschitz estimates for $$\Lambda $$ itself, and these estimates are fairly tight around the corners $$-1 $$ and $$1 $$ [[Proposition 2.6, JSS24]](). For instance, using the coordinate transfer scheme in conjunction with the MGF bound strategy outlined above, we can show that

$$ \begin{equation} |\partial_y \Lambda(t,y)| \le \frac{1}{2}\log\left(\frac{2}{1-|y|}\right) + 4\beta^2(1-t)\, , \end{equation} $$

which tells us that the gradient of $$\Lambda $$ blows up logarithmically as $$y \to -1/1 $$. This immediately allows us to obtain a uniform continuity estimate on $$\Lambda $$ itself by some integration and the application of the fundamental theorem of calculus,

$$ \begin{equation} |\Lambda(t,y') - \Lambda(t,y)| \le \frac{1}{2}|y-y'|\left(\log\left(\frac{2}{|y - y'|}\right) + 1 + 8\beta^2(1-t)\right)\,. \end{equation} $$

_<u>Inf-convolution and convergence of </u>$$\Lambda_\gamma \to \Lambda $$_: By definition, it is clear that $$\Lambda_\gamma = \Lambda $$ when $$\gamma = 0 $$. However, for the analysis, we need a quantitative estimate on the uniform convergence of the former to the latter. In doing this, we use the fact that, by definition, $$\Lambda_\gamma $$ is a convolution of the concave function $$-\frac{\gamma}{2}x^2$$ with the function $$\Phi $$ and, therefore, will obey an **inf-convolution** rule for FL duals. This, in conjunction with the continuity estimates for $$\Lambda $$ above, will allow us to obtain the desired quantitative uniform convergence estimates.

The _inf-convolution_ formula tells us that,

$$ \begin{equation} \Lambda_\gamma(t,y) = \inf_{y' \in [-1,1]} \left(\Lambda(t,y') + \frac{1}{2\gamma}(y'-y)^2\right)\, .\end{equation} $$

This is proved by noticing that since $$\partial_x \Phi $$ is strictly increasing in $$x $$, so is $$\partial_y \Lambda = (\partial_x \Phi)^{-1} $$; this implies that $$(y + \gamma \partial_y\Lambda(t,y)) $$ is also strictly increasing and so there is a unique point $$y' $$, such that, $$y' + \partial_y\Lambda(t,y') = y $$, which is a critical point of the function inside the infimum to minimize.

At this point, substituting

$$ \begin{equation} y = y' + \gamma\partial_y\Lambda_\gamma(t,y)\,, \end{equation} $$    

into the expression for $$\Lambda_\gamma = xy - \Phi_\gamma(t,x)$$ and some algebra concludes the inf-convolution formula. We now combine the inf-convolution formula with the uniform continuity estimates obtained for $$\Lambda $$ above to immediately conclude (with some algebra) that

$$ \begin{equation} \Lambda - (1+4\beta^2)(2\beta^2\gamma)^{(4\beta^2)/(1+4\beta^2)}\le \Lambda_\gamma \le \Lambda\,, \end{equation}$$

which will be a perfectly fine estimate for a choice of $$\gamma $$ that is exponentially small in $$\beta $$, and $$\beta = O(1/\epsilon) $$, with $$\epsilon > 0 $$ coming from the desired $$(1-\epsilon) $$-approximation ratio.

At this point, we are ready to state the three main concluding results from this section:
* The (smoothed) primal Parisi PDE and AC SDE in terms of $$\Lambda_\gamma $$,
* The $$L^2 $$-distance between the processes generated by the AC SDE of $$\Lambda $$ and $$\Lambda_\gamma $$, and
* Various estimates for the derivatives of $$\Lambda_\gamma $$ that are proved using the RPC representation for $$\Phi $$ that I will not discuss in this post[^2].

_<u>The PDE, SDE & Wasserstein bounds between the primal SDEs</u>_: We begin with the PDE, and then move on the SDE(s) and the relevant Wasserstein distance bounds between them.

Using the Crouzeix identity, the coordinate transforms implied by the FL duality, the characterization of the critical points and a substitution of the Parisi PDE, we can obtain the following PDE for $$\Lambda_\gamma $$,

$$ \begin{equation} \partial_t \Lambda_\gamma(t,y) = \beta^2\left(\frac{1}{\partial_{y,y}\Lambda_\gamma(t,y)} - \gamma + \mu(t)\left(y - \gamma \partial_y \Lambda(t,y)\right)^2\right)\end{equation}\,. $$

Recall that in [(1.1)](#the-parisi-formula-and-auffinger-chen-representation) we mentioned that it would be critical to understand the rate at which $$\Lambda_\gamma $$ changes at various points in time: for this purpose, the PDE above becomes critical. Using the PDE above, we can express the time derivative above in terms of the spatial derivatives, and this will be important in the Taylor expansion analysis that bounds the fluctuations of the modified objective function at each step. In fact, in the Taylor expansion analysis, we will isolate the RHS of the PDE into three components:
* The first will consist of the $$\frac{1}{\partial_{y,y}\Lambda_\gamma(t,y)} $$ term, and this will be a critical factor in terms of its interaction with the Hessian term (as we also saw in [(1.1)](#the-parisi-formula-and-auffinger-chen-representation)).
* The second will consist of just $$\gamma $$, and if $$\Delta t$$ is sufficiently small (along with $$\gamma $$ itself being small), this will simply be a small error term at each time step.
* The third and final term will be handled, as also mentioned in [(1.1)](#the-parisi-formula-and-auffinger-chen-representation), using convergence of the iterates to the primal AC SDE and shown to be sufficiently small at every step.

Using a bit of Ito calculus, as is done in [[Proposition 2.10, JSS24]](https://arxiv.org/pdf/2408.02360), the following primal AC SDEs for $$\Lambda $$ and $$\Lambda_\gamma $$ are obtained,

$$ \begin{equation} dY_t = \frac{\sqrt{2}\beta}{\partial_{y,y}\Lambda(t,y)}dW_t\,, Y_0 = 0\,\,\,\text{and}\,\,\,\,\, dY^{\gamma}_t = \frac{\sqrt{2}\beta}{\partial_{y,y}\Lambda_\gamma(t,y)}dW_t\,,\,Y^{\gamma}_0 = 0\,.\end{equation} $$

Then, using the estimates for the derivatives of $$\Lambda $$ and $$\Lambda_\gamma $$ (as briefly mentioned below) as Lipschitz bounds in conjunction with Gronwall's inequality (again)[^3] and some Ito calculus, one obtains the following final bound [[Lemma 2.17, JSS24](https://arxiv.org/pdf/2408.02360)],

$$ \begin{equation} \| Y^{\gamma}_t - Y_t\|^2_{L^2} \le 2\gamma^2\left(e^{10\beta^2 t} - 1\right) \,.\end{equation} $$

This estimate will be crucially used to understand the accumulation of error we obtain from following the Hessian of a TAP corrected energy with $$\Lambda_\gamma $$, as opposed to the true TAP corrected energy which invokes $$\Lambda $$ directly.

_<u>Derivative estimates</u>_: We now conclude the first post with a small statement that gives various estimates for the derivatives of $$\Lambda_\gamma $$.

The spatial derivatives are bounded as,

$$ \begin{equation} \frac{1}{1+\gamma} \le \partial_{y,y}\Lambda_\gamma \le \frac{1}{\gamma}\,,\,|\partial_{y,y,y}\Lambda_\gamma| \le \frac{2}{\gamma^2} \,,\end{equation} $$

and the spatial and temporal Lipschitz estimates for the driver of the "smoothed" primal AC SDE are,

$$ \begin{equation}|\partial_y\left(\frac{1}{\partial_{y,y}\Lambda_\gamma}\right)| \le 2\,,\,|\partial_t\left(\frac{1}{\partial_{y,y}\Lambda_\gamma}\right)| \le 14\beta^2 \,.\end{equation} $$

These estimates are stated with careful formal precision in [[Proposition 2.11, JSS24](https://arxiv.org/pdf/2408.02360)] and rely primarily on the representation of $$\Phi $$ using the RPCs, and working with the polynomial representations induced by them.

#### FOOTNOTES

[^1]: When David, Jonathan and I were working on the project in the early days, we were thinking of this convex duality via the lens of mirror maps. It turns out that, because of the underlying convex duality, there is a connection to be made here through the information geometry of the underlying Hessian being akin to a (flat) Bregmannian manifold. However, this is beyond the scope of this post and something we are investigating in ongoing work.

[^2]: In a previous [post](https://juspreetsandhu.me/2022/02/08/the-sk-model-i#ruelle-probability-cascades), I described a _part_ of the RPC construction and the details afforded there are gentle and sufficient enough to understand the main ingredients in their construction, as well as their purpose. A slightly more detailed overview of the construction, along with how exactly the representation gets used in the [Hopf-Cole transform]() to solve the Parisi PDE for _atomic_ measures is provided in [[Appendix C, JSS24]](https://arxiv.org/pdf/2408.02360) and can be read by the interested reader. The estimates for the derivatives are proved between [[Lemma 2.12 & Lemma 2.13, JSS24]](https://arxiv.org/pdf/2408.02360) and stated in [[Proposition 2.11, JSS24]](https://arxiv.org/pdf/2408.02360).  

[^3]: For converting Lipschitz bounds into estimates of how much error propagates over a period of time, [Gronwall's inequality](https://en.wikipedia.org/wiki/Gr%C3%B6nwall%27s_inequality) is an indispensable tool that we use in Sections 2 & 5 of the paper. The combination of Ito's lemma & Gronwall's inequality allows us to, in fact, more or less have a mechanistic procedure for converting Lipschitz estimates into error bounds of various sorts.

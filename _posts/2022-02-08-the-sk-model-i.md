---
layout: post
title:  "The SK Model I"
date:   2022-02-08 11:15:40
blurb: "The Aizenman-Sims-Starr Representation & Ruelle Probability Cascades"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)
<br />
The goal of this post is to introduce the _Sherrington-Kirkpatrick_ model, which is the canonical starting point in Spin-Glass Theory. This model (along with a slightly generalized family) is one of the first models for which the infamous _Parisi-Variational Principle_ was (formally) proven to be true. It is also a model that instigated the introduction of many other concepts via the _Replica-Symmetry Breaking_ ansatzen such as _Ultrametricity_, the _TAP Equations_, and the _Ghirlanda-Guerra Identities_. It is the prototypical model of a <strong>Mean-Field Spin Glass</strong> and has recently had algorithmic implications, which in conjunction with the growing body of work on the _Overlap-Gap Property_ have shed light on the average-case complexity of a large family of optimization problems[^1].

<br />


#### Table of Contents
1. [The Sherrington-Kirkpatrick Model](#the-sherrington-kirkpatrick-model)
   * [An Average-Case Optimization Problem](#an-average-case-optimization-problem)
   * [Covariance and Overlaps](#covariance-and-overlaps)
   * [Guerra-Tonnineli Interpolation](#guerra-tonnineli-interpolation)
   * [Gaussian Concentration](#gaussian-concentration)
2. [Aizenman-Sims-Starr Scheme](#aizenman-sims-starr-scheme)
    * [The ASS Functional](#the-ass-functional)
    * [Invariance Symmetries](#invariance-symmetries)
3. [Ruelle Probability Cascades](#ruelle-probability-cascades)
4. [Footnotes](#footnotes)

<br />

## The Sherrington Kirkpatrick Model
We briefly introduce the Sherrington-Kirkpatrick model as an optimization problem. We will be interested in the almost-surely limit of the optimal value, and we will use an elementary estimate to establish the appropriate normalization in searching for this limit value. The covariance of the optimization problem under consideration will be computed and shown to be a function of the _overlap_ between two solutions - This is a fundamental and critical observation that will show up many times as a quantity of interest in the proof of the Parisi Variational Principle. We will introduce the infamous [Guerra-Tonnineli interpolation](#guerra-tonnineli-interpolation) and use it to demonstrate the existence of (but not compute exactly) the "smoothed" verison of the limit we are interested in _on average_. Using the interpolation again in conjunction with some convexity arguments, we will then establish a [gaussian concentration inequality](#gaussian-concentration) which will imply that the limit exists almost-surely.
<br />

### An Average-Case Optimization Problem
Given $$n^2$$ i.i.d. standard normal ($$\mathcal{N}(0, 1) $$) variables $$\{J_{ij}\}_{i, j \in [n]}$$, we are interested in the optimal value of the following optimization problem over the hypercube,

$$ \begin{equation} \max_{\sigma \in \{\pm 1\}^n} \frac{1}{\sqrt{n}}\sum_{i, j=1}^n J_{i,j}\sigma_i\sigma_j\, . \end{equation}$$

The problem above can be seen as asking for the $$\mathsf{MAX}$$-$$\mathsf{CUT}$$ of a complete graph with i.i.d. $$\mathcal{N}(0,1)$$ weights. As mentioned, we will be interested in the _typical_ optimal value of the above problem as $$n$$ gets large. The normalized quadratic form above will be expressed as $$H_n(\sigma)$$ which, in physics language, translates to asking for the value of the hamiltonian $$H_n$$ under the configuration of spins $$\sigma \in \{\pm 1\}^n$$. Therefore, we will ask the following question instead,

$$ \begin{equation} \lim_{n \to \infty} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma) \right]\, . \end{equation}$$

While this problem seems well defined, we need to nornalize it carefully so that the limit we wish to compute does not diverge.

<br />

### Covariance and Overlaps
A simple observation (formalized below) will reveal that the fluctuations (variance) of the term above lead to a divergent limit. Therefore, we must normalize it appropariately. To know what normalization is appropriate, we compute the covariance of the underlying gaussian process explicitly,

$$ \begin{align} \mathbb{E}[H_n(\sigma^1)H_n(\sigma^2)] &= \frac{1}{n}\mathbb{E}\left[ \sum_{i_1, j_1}J_{i_1, j_1}\sigma^1_{i_1}\sigma^2_{j_1}\sum_{i_2, j_2}J_{i_2, j_2}\sigma^1_{i_2}\sigma^2_{j_2}\right] \\ &= \frac{1}{n}\mathbb{E}\left[\sum_{i_1, j_1}J^2_{i_1, j_1}\sigma^1_{i_1}\sigma^2_{i_1}\sigma^1_{j_1}\sigma^2_{j_1}\right] \\ &= \frac{1}{n}\left(\sum_{i =1}^n \sigma^1_i\sigma^2_i\right)^2 = n\cdot \left(\frac{1}{n}\langle \sigma^1, \sigma^2\rangle\right)^2\, . \end{align} $$

The above calculation reveals that the covariance of the quantity we wish to maximize is a $$O(n)$$ scaling of the _normalized overlap_ between two configurations $$\sigma^1 $$ and $$\sigma^2 $$. Since the latter quantity is a constant, this tells us that the extra normalization term we wish to add divides the quantity to optimize by $$n$$. This will yield the final average quantity whose limit we wish to compute precisely,

$$ \begin{equation} \lim_{n \to \infty} \frac{1}{n} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma)\right]\, , \end{equation} $$

and this quantity will be termed the _Ground State Energy_ of the system. When studying such quantities, it is a standard trick in statistical physics to understand this quantity using its _smoothed_ approximation, known as the _Free Energy Density_ of the system. This quantity depends on the smoothing parameter $$\beta$$ and is defined as,

$$ \begin{equation} F_{n, \beta} := \frac{1}{n}\mathbb{E}\left[\log\left(\sum_{\sigma \in \{ \pm 1\}^n} e^{\beta H_n(\sigma)}\right)\right] \, , \end{equation} $$

where the exponential summation term,

$$ \begin{equation} Z_{n, \beta} := \sum_{\sigma \in \{\pm 1\}^n}e^{\beta H_n(\sigma)}\, , \end{equation} $$

is termed the _Partition function_ in Statistical Physics. As we shall see, this softmax trick induces a measure called the _Gibbs Measure_ which will turn out to be very closely related to the _Free Energy_ of the model and will appear as soon as a derivative of the latter is taken with respect to an appropriate parameter in the hamiltonian. Understanding the asymptotic geometric structure of the gibbs, therefore, will be of paramount importance in formally proving the Parisi Variational Principle.

The smoothed free energy can be related to the ground state energy in the regime that the smoothing parameter $$\beta \to \infty$$. This is made precise in the following,

<u><strong>(Proposition-1)</strong></u>: The following holds for all $$\beta > 0$$,

$$ \begin{equation} \lim_{n \to \infty} \frac{1}{n} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma)\right] \leq \lim_{n \to \infty} \frac{1}{\beta n}\mathbb{E}\left[\log\left(\sum_{\sigma \in \{\pm 1\}^n}e^{\beta H_n(\sigma)}\right)\right] \leq \lim_{n \to \infty} \frac{1}{n} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma)\right] + \frac{\log(2)}{\beta} \, . \end{equation} $$

Note that the above immediately implies (via an application of Holder's inequality) that,

$$ \begin{equation} \lim_{n \to \infty} \frac{1}{n} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma)\right] = \lim_{\beta \to \infty}\frac{1}{\beta}\left(\lim_{n \to \infty} F_{n, \beta}\right)\, , \end{equation} $$

if we assume that the $$n \to \infty$$ limit (also called the _thermodynamic limit_) of the free energy density exists.
<br />

### Guerra-Tonnineli Interpolation
We now introduce a smooth interpolation between three independent instances of the SK model, and in conjunction with the lemmas of _Fekete_ and _Stein_, use it show that the thermodynamic limit of the free energy density exists. <br />

Before introducing the interpolation formally, we briefly state _Fekete's_ and _Stein's_ lemmata. The proofs for these are elementary and require no more than basic algebraic manipulations and integration by parts, and are omitted in this blog post[^2].

<u><strong>(Stein's Lemma)</strong></u>: Given a differentiable function $f$ that doesn't grow too fast, and a jointly-gaussian process $$\{g(\sigma)\}_{\sigma \in \Sigma} $$, the following holds,

$$ \begin{equation} \mathbb{E}[g_{\sigma}f(g)] = \sum_{\sigma' \in \Sigma}\mathbb{E}\left[g_{\sigma}g_{\sigma'} \right]\mathbb{E}\left[\partial_{\sigma'}f(g)\right]\, . \end{equation}$$

<u><strong>(Fekete's Lemma)</strong></u>: If a sequence $$\{x_n\}_{n = 1}^{\infty} $$ is _superadditive_ ($$x_n + x_m \leq x_{n + m}\, ,\, \forall n, m \geq 1 $$), then,

$$ \begin{equation} \lim_{n \to \infty} \frac{x_n}{n} = \sup_{n \geq 1}\frac{x_n}{n}\, . \end{equation}$$

We now state the interpolation, and give some intuition about its structure before proving the main lemma which asserts the existence of the limit for the free energy density as $$n \to \infty $$. <br />

<u><strong>(Guerra-Tonnineli Interpolation)</strong></u>: Given three independent instances of the SK model on $$n, m$$ and $$n + m$$ particles, the interpolation is defined $$\forall t \in [0, 1]$$ as follows,

$$ \begin{equation} H^t(\sigma) = \sqrt{t} H_{n + m}(\rho\cdot\tau) + \sqrt{1-t}\left(H_n(\rho) + H_m(\tau)\right)\, . \end{equation} $$

Notice that the interpolation above is equivalent to two independent copies of the SK model of size $$n$$ and $$m$$ at $$t = 0 $$, and becomes a single copy of the SK model of size $$n + m$$ at $$t=1 $$. The square-roots are introduced to keep the variance of a single gaussian interaction in the interpolated gaussians within the graphs of size $$n$$ and $$m$$ to $$1 $$, while slowly increasing the variance of the gaussian interactions _between_ the two graphs from $$0$$ to $$1$$. As we shall see, the interpolation helps prove that the free energy is _superadditive_ and an application of _Fekete's_ lemma immediately implies the existence of the limit.

<u><strong>(Lemma-1)</strong></u>: The free energy exists in the thermodynamic limit: $$\lim_{n \to \infty} F_{n, \beta}$$ exists. <br />
<u><i>Proof</i></u>: We will analyze the free energy density $$\phi(t)$$ of the interpolated hamiltonian $$H^t$$ at every $$t \in [0, 1] $$. Note that since $$J_{i, j}$$ are continuously distributed and the expectation is a convex combination of continuous variables, $$\phi(t) $$ is continuous. The change of free energy density as a function of $$t $$ is then given as,

$$ \begin{align} \partial_t\phi(t) &= \frac{1}{n + m}\mathbb{E}\left[ \frac{1}{Z_t}\partial_t(Z_t)\right] = \frac{1}{n + m}\mathbb{E}\left[ \frac{1}{Z_t}\left(\sum_{\sigma \in \{\pm 1\}^n}\partial_te^{ H^t(\sigma)}\right)\right] \\ &= \frac{1}{n + m}\mathbb{E}\left[ \langle \partial_t H^t(\sigma)\rangle_t\right]\end{align} $$

where $$\langle . \rangle_t $$ denotes the average with respect to the Gibbs measure at $$t $$. The above relationship states that the free energy density changes proportional to the average rate of change of the energy of the interpolated hamiltonian. The above expression is interpolated by applying a wonderful lemma that allows us to rewrite the expected value of some jointly gaussian vector $$\{x_{\sigma}\} $$ in terms of a term that subtracts the "covariance" between $$\{x_{\sigma}\} $$ and another gaussian vector $$\{y_{\sigma}\} $$ from the "overlap" terms.

<u><strong>(Gaussian Covariance for Gibbs Average </strong></u>[\[Lemma 1.1, Pa14\]](https://link.springer.com/book/10.1007/978-1-4614-6289-7)<u><strong>)</strong></u>: Given two jointly gaussian vectors $$\{x_{\sigma}\} $$ and $$\{y_{\sigma}\} $$, the following can be said about the iterated average of $$x_{\sigma} $$ (where the average is with respect to the Gibbs for $$\{y_{\sigma}\} $$),

$$ \begin{equation} \mathbb{E}[\langle x_\sigma \rangle] = \mathbb{E}[\langle \mathbb{E}[x_{\sigma_1}y_{\sigma_1}]\rangle - \langle \mathbb{E}[x_{\sigma_1}y_{\sigma_2}]\rangle]\end{equation}\, . $$


The lemma above allows us to re-write the average rate of change in the ground state energy of the interpolated hamiltonian in terms of a covariance term with respect to the hamiltonian itself, which can be expanded and evaluated in terms of the overlaps. Setting $$x_{\sigma} = \partial_t H^t(\sigma) $$ and $$y_{\sigma} = H^t(\sigma) $$, we have

$$ \begin{align} \frac{1}{n + m}\mathbb{E}\left[ \langle \partial_t H^t(\sigma)\rangle_t\right] &= \frac{1}{n + m}\mathbb{E}\left[\langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right]\rangle - \langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle\right] \end{align}\, .$$

The terms inside the expectations can easily be evaluated by taking the derivatives and multiplying to be,

$$ \begin{align} &\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right] = \mathbb{E}\left[ \frac{1}{2}(H_{n + m}(\sigma_1))^2 - (H_n(\rho_1) + H_m(\tau_1))^2\right] \\
&= \mathbb{E}\left[(n + m)\left(\frac{\langle\sigma_1, \sigma_1\rangle}{n + m}\right)^2 - n\left(\frac{\langle \rho_1, \rho_1 \rangle}{n}\right)^2 - m\left(\frac{\langle \tau_1, \tau_1 \rangle^2}{m}\right)\right] \\ &= 0\, ,\end{align} $$

where we used the fact that $$H_n $$ and $$H_m $$ are independent to zero out the covariance term in-between them in conjunction with the fact that the [covariance of the underlying gaussian process is the square of its normalized overlaps](#covariance-and-overlaps).

A similar calculation immediately reveals that,
$$ \begin{align} &\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right] = \mathbb{E}\left[\frac{1}{2}(H_{n + m}(\sigma_1)H_{n + m}(\sigma_2)) - H_n(\rho_1)H_n(\rho_2) - H_m(\tau_1)H_m(\tau_2)\right] \\ &= \frac{1}{2}\mathbb{E}\left[ (n + m)\left(\frac{\langle \sigma_1, \sigma_2 \rangle}{n + m}\right)^2 - n\left(\frac{\langle \rho_1, \rho_2 \rangle}{n}\right)^2  - m\left(\frac{\langle\tau_1, \tau_2 \rangle}{m}\right)^2\right] \\ &= \frac{1}{2}\mathbb{E}\left[ \frac{n}{n+m}\left(\frac{n}{n + m} - 1\right)\left(\frac{\langle \rho_1, \rho_2 \rangle}{n}\right)^2 + \left(\frac{m}{n + m} - 1\right)\left(\frac{\langle \sigma_1, \sigma_2 \rangle}{m}\right)^2\right]\end{align} $$

Note that the term above is negative since it is equivalent to

$$ \begin{align} \frac{1}{m + n}\left(2AB - mA - nB\right)\, ,\forall n, m \geq 1 \end{align} $$

where $$A = n\left(\frac{\langle \rho_1, \rho_2\rangle}{n}\right)^2$$ and $$B = m\left(\frac{\langle \tau_1, \tau_2\rangle}{m}\right)^2 $$, and we used the fact that $$\|A\|, \|B\| \leq 1 $$ and 
$$ \begin{equation} \langle \sigma_1, \sigma_2 \rangle = \frac{n}{m + n}\left(\frac{\langle\rho_1, \rho_2\rangle}{n}\right) + \frac{m}{m + n}\left(\frac{\langle \tau_1, \tau_2 \rangle}{m}\right)\, . \end{equation}$$

Note that the above immediately implies that,

$$ \begin{equation} \partial_t\phi(t) = \frac{1}{n + m}\mathbb{E}\left[\langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right]\rangle - \langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle\right] = \langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle\right] \geq 0\, . \end{equation} $$

The fact that the square of the overlaps is a convex function in conjunction with the observations that,

$$ \begin{align} \phi(0) = \frac{n}{m + n}F_n + \frac{m}{m + n}F_m \, , \\ \phi(1) = F_{n + m} \, , \end{align} $$

yields that the free energy is superadditive. This yields that the limit of the free energy density is well-defined in the thermodynamic limit with an application of _Fekete's Lemma_. <br />

### Gaussian Concentration
The goal of this section is to "boost" the previous lemma showing that the free energy density of the SK model is well-defined on average to an _almost-surely_ statement. In order to do that it is crucial to prove that there is concentration of the <strong>log-partition</strong> function (under the gaussians).
<br />

## Aizenman-Sims-Starr Scheme
<br />

### The ASS functional
<br />

### Invariance Symmetries
<br />

## Ruelle Probability Cascades
<br />

#### FOOTNOTES

[^1]: There are many references to this body of work which will be given in due course. Nonetheless, the following two surveys are a nice start: [\[B05\]](http://www.numdam.org/item/SB_2004-2005__47__349_0.pdf), [\[G21\]](https://arxiv.org/pdf/2109.14409.pdf).

[^2]: Stein's Lemma can be proved by a simple integration-by-parts argument for every coordinate $$i $$ in conjunction with Fubini's theorem and the chain-rule. The lemma itself simply asserts that the correlation between a gaussian variable and some function of it is equivalent to a sum of scalings (by the correlations) of the average rate of change of the function. Likewise, Fekete's Lemma boils down to algebraic manipulation of the sequence, where we take the limit infimum of $$x_n/n $$ and compare it to limit infimum of some $$x_m $$ where $$m $$ is chosen to be the supremum as a divisor with a remainder term. The lemma itself merely asserts that for an appropriately growing sequence, the empirical average in the limit simply picks out the largest contributing term.

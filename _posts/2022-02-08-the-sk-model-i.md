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
    * [The Parisi Variational Principle](#parisi-variational-principle)
    * [The RPC Construction](#rpc-construction)
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

While this problem seems well defined, we need to nornalize it carefully so that the limit we wish to compute does not diverge. <br />

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

<u><strong>(Stein's Lemma)</strong></u>: Given a differentiable function $$f $$ that doesn't grow too fast, and a jointly-gaussian process $$\{g(\sigma)\}_{\sigma \in \Sigma} $$, the following holds,

$$ \begin{equation} \mathbb{E}[g_{\sigma}f(g)] = \sum_{\sigma' \in \Sigma}\mathbb{E}\left[g_{\sigma}g_{\sigma'} \right]\mathbb{E}\left[\partial_{\sigma'}f(g)\right]\, . \end{equation}$$

<u><strong>(Fekete's Lemma)</strong></u>: If a sequence $$\{x_n\}_{n = 1}^{\infty} $$ is _superadditive_ ($$x_n + x_m \leq x_{n + m}\, ,\, \forall n, m \geq 1 $$), then,

$$ \begin{equation} \lim_{n \to \infty} \frac{x_n}{n} = \sup_{n \geq 1}\frac{x_n}{n}\, . \end{equation}$$

We now state the interpolation, and give some intuition about its structure before proving the main lemma which asserts the existence of the limit for the free energy density as $$n \to \infty $$. <br />

<u><strong>(Guerra-Tonnineli Interpolation)</strong></u>: Given three independent instances of the SK model on $$n, m$$ and $$n + m$$ particles, the interpolation is defined $$\forall t \in [0, 1]$$ as follows,

$$ \begin{equation} H^t(\sigma) = \sqrt{t} H_{n + m}(\rho\cdot\tau) + \sqrt{1-t}\left(H_n(\rho) + H_m(\tau)\right)\, . \end{equation} $$

Notice that the interpolation above is equivalent to two independent copies of the SK model of size $$n$$ and $$m$$ at $$t = 0 $$, and becomes a single copy of the SK model of size $$n + m$$ at $$t=1 $$. The square-roots are introduced to keep the variance of a single gaussian interaction in the interpolated gaussians within the graphs of size $$n$$ and $$m$$ to $$1 $$, while slowly increasing the variance of the gaussian interactions _between_ the two graphs from $$0$$ to $$1$$. As we shall see, the interpolation helps prove that the free energy is _superadditive_ and an application of _Fekete's_ lemma immediately implies the existence of the limit.

<u><strong>(Lemma-1)</strong></u>: The free energy exists in the thermodynamic limit: $$\lim_{n \to \infty} F_{n, \beta}$$ exists. <br />
<u><i>Proof</i></u>: We will analyze the free energy density $$\phi(t)$$ of the interpolated hamiltonian $$H^t$$ at every $$t \in [0, 1] $$. Note that since $$J_{i, j}$$ are continuously distributed and the expectation is a convex combination of continuous variables, $$\phi(t) $$ is continuous. The change of free energy density as a function of $$t $$ is then given as,

$$ \begin{align} \partial_t\phi(t) &= \frac{1}{n + m}\mathbb{E}\left[ \frac{1}{Z_t}\partial_t(Z_t)\right] = \frac{1}{n + m}\mathbb{E}\left[ \frac{1}{Z_t}\left(\sum_{\sigma \in \{\pm 1\}^n}\partial_te^{ H^t(\sigma)}\right)\right] \\ &= \frac{1}{n + m}\mathbb{E}\left[ \langle \partial_t H^t(\sigma)\rangle_t\right]\, ,\end{align} $$

where $$\langle . \rangle_t $$ denotes the average with respect to the Gibbs measure at $$t $$. The above relationship states that the free energy density changes proportional to the average rate of change of the energy of the interpolated hamiltonian. The above expression is evaluated by applying a wonderful lemma that allows us to rewrite the expected value of some jointly gaussian vector $$\{x_{\sigma}\} $$ in terms of a term that subtracts the "covariance" between $$\{x_{\sigma}\} $$ and another gaussian vector $$\{y_{\sigma}\} $$ from the "overlap" terms[^3].

<u><strong>(Gaussian Covariance for Gibbs Average <a href="https://link.springer.com/book/10.1007/978-1-4614-6289-7">[Lemma 1.1, Pa14]</a>)</strong></u>: Given two jointly gaussian vectors $$\{x_{\sigma}\} $$ and $$\{y_{\sigma}\} $$, the following can be said about the iterated average of $$x_{\sigma} $$ (where the average is with respect to the Gibbs for $$\{y_{\sigma}\} $$),

$$ \begin{equation} \mathbb{E}[\langle x_\sigma \rangle] = \mathbb{E}[\langle \mathbb{E}[x_{\sigma_1}y_{\sigma_1}]\rangle - \langle \mathbb{E}[x_{\sigma_1}y_{\sigma_2}]\rangle]\end{equation}\, . $$


The lemma above allows us to re-write the average rate of change in the ground state energy of the interpolated hamiltonian in terms of a covariance term with respect to the hamiltonian itself, which can be expanded and evaluated in terms of the overlaps. Setting $$x_{\sigma} = \partial_t H^t(\sigma) $$ and $$y_{\sigma} = H^t(\sigma) $$, we have

$$ \begin{align} \frac{1}{n + m}\mathbb{E}\left[ \langle \partial_t H^t(\sigma)\rangle_t\right] &= \frac{1}{n + m}\mathbb{E}\left[\langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right]\rangle_t - \langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle_t\right] \end{align}\, .$$

The terms inside the expectations can easily be evaluated by taking the derivatives and multiplying to be,

$$ \begin{align} &\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right] = \mathbb{E}\left[ \frac{1}{2}(H_{n + m}(\sigma_1))^2 - (H_n(\rho_1) + H_m(\tau_1))^2\right] \\
&= (n + m)\left(\frac{\langle\sigma_1, \sigma_1\rangle}{n + m}\right)^2 - n\left(\frac{\langle \rho_1, \rho_1 \rangle}{n}\right)^2 - m\left(\frac{\langle \tau_1, \tau_1 \rangle^2}{m}\right) \\ &= 0\, ,\end{align} $$

where we used the fact that $$H_n $$ and $$H_m $$ are independent to zero out the covariance term in-between them in conjunction with the fact that the [covariance of the underlying gaussian process is the square of its normalized overlaps](#covariance-and-overlaps).

A similar calculation immediately reveals that,

$$ \begin{align} &\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right] = \mathbb{E}\left[\frac{1}{2}(H_{n + m}(\sigma_1)H_{n + m}(\sigma_2)) - H_n(\rho_1)H_n(\rho_2) - H_m(\tau_1)H_m(\tau_2)\right] \\ &= \frac{1}{2}(n + m)\left(\frac{\langle \sigma_1, \sigma_2 \rangle}{n + m}\right)^2 - n\left(\frac{\langle \rho_1, \rho_2 \rangle}{n}\right)^2  - m\left(\frac{\langle\tau_1, \tau_2 \rangle}{m}\right)^2 \\ &= \frac{1}{2}\frac{n}{n+m}\left(\frac{n}{n + m} - 1\right)\left(\frac{\langle \rho_1, \rho_2 \rangle}{n}\right)^2 + \left(\frac{m}{n + m} - 1\right)\left(\frac{\langle \sigma_1, \sigma_2 \rangle}{m}\right)^2 + \frac{2mn}{m + n}\left(\frac{\langle\rho_1,\rho_2\rangle\langle\tau_1,\tau_2\rangle}{mn}\right)^2 \end{align} $$

Note that the term above is negative since it is equivalent to

$$ \begin{align} \frac{1}{m + n}\left(2AB - mA - nB\right)\, ,\forall n, m \geq 1\, , \end{align} $$

where $$A = n\left(\frac{\langle \rho_1, \rho_2\rangle}{n}\right)^2$$ and $$B = m\left(\frac{\langle \tau_1, \tau_2\rangle}{m}\right)^2 $$, and to assert negativity we used the facts that $$|A|, |B| \leq 1 $$ and
$$ \begin{equation} \langle \sigma_1, \sigma_2 \rangle = \frac{n}{m + n}\left(\frac{\langle\rho_1, \rho_2\rangle}{n}\right) + \frac{m}{m + n}\left(\frac{\langle \tau_1, \tau_2 \rangle}{m}\right)\, . \end{equation}$$

This immediately implies that,

$$ \begin{align} \partial_t\phi(t) &= \frac{1}{n + m}\mathbb{E}\left[\langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_1) \right]\rangle_t - \langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle_t\right] \\ &= \frac{-1}{m+n}\mathbb{E}\left[\langle\mathbb{E}\left[\partial_tH^t(\sigma_1)H^t(\sigma_2)\right]\rangle_t\right] = -\mathbb{E}\left[\langle\textsf{g}(\text{overlaps}^2)\rangle_t\right] \geq 0\, , \end{align} $$

since the expectations of non-negative random variables are non-negative. Additionally, these are convex combinations of the squares of the overlap, which are convex functions themselves. Therefore, $$\partial_t\phi(t)$$ is a convex function. Then, the observations that,

$$ \begin{align} &\phi(0) = \frac{n}{m + n}F_n + \frac{m}{m + n}F_m \, , \\ &\phi(1) = F_{n + m} \, , \end{align} $$

yield that the free energy is superadditive (since $$\partial_t\phi(t) \geq 0 $$ and convexity implies $$\phi(0) \leq \phi(1) $$). This shows that the limit of the free energy density is well-defined in the thermodynamic limit with an application of _Fekete's Lemma_. <br />

### Gaussian Concentration
The goal of this section is to "boost" the previous lemma showing that the free energy density of the SK model is well-defined on average to an _almost-surely_ statement. In order to do that it is crucial to show that there is concentration of the <strong>log-partition</strong> function (under the gaussians). To do this, we will use a gaussian concentration inequality that is proved using a _gaussian interpolation_ technique (which is very similar to the [Guerra-Tonnineli interpolation](#guerra-tonnineli-interpolation)) in conjunction with some elementary convexity properties[^3].

<u><strong>(Gaussian Concentration Inequality)</strong></u>: Given a $$b $$-lipschitz function $$F:\mathbb{R}^n \to \mathbb{R} $$, and a jointly gaussian process $$\{g_i\}_{i=1}^n$$ with bounded covariance $$C $$, the following holds $$\forall \epsilon > 0 $$

$$ \begin{equation} \Pr_g\left[\left|F(g_1,\dots,g_n) - \mathbb{E}_g\left[F(g_1,\dots,g_n)\right]\right| \geq \epsilon\right] \leq 2e^{-\frac{\epsilon^2}{4b^2 C}}\, .  \end{equation}$$

Note that choosing $$F(\{g_{\sigma}\}) = \log\left(\sum_{\sigma \in \{\pm 1\}^n}e^{g_{\sigma}}\right) = \log(Z_n)$$ with $$C = \mathbb{E}[H_n(\sigma)^2] = 1 $$ and $$b = \sqrt{C} = 1$$, immediately yields concentration for the partition function/free energy density as,

$$ \begin{equation} \Pr_g\left[\left|\frac{1}{n}\log\left(Z_n\right) - \frac{1}{n}\mathbb{E}\left[\log(Z_n)\right]\right| \geq \epsilon \right] \leq 2e^{-\frac{\epsilon^2 n}{4}}\, .\end{equation} $$

With the result above, we can argue that the limit of the _free energy density_ exists _almost-surely_ for large enough system sizes. <br />

## Aizenman-Sims-Starr Scheme
We now introduce a scheme that will derive an expression for the change in the free energy of the system when a "cavity" is created, or equivalently, we decrease the size of the input instance by removing one vertex. Rewriting the free energy as a telescoping sum over the cavities then allows one to get an expression for the free energy - This term is known as the _Aizenman-Sims-Starr (ASS) functional_. Unfortunately, it is not possible to show that the limit exists for the ASS functional, but we can use it to _lower bound_ the free energy of the SK model by showing that its $$\lim\inf $$ exists. The close relationship of the Gibbs measure with the free energy is made clear in two ways:
1. The ASS functional depends on the Gibbs measure, and a telescoping sum shows that the free energy is a Gibbs-averaged quantity.
2. We will also state a result that reveals a close relationship between the structure of the _overlap distribution_ induced by the Gibbs measure and the _ASS functional_ - As will become clearer later, this is a deep result that hints at a tight connection between the two seemingly different approaches of the _Replica Method_ and _Cavity Method_ used extensively in Statistical Physics. <br />

### The ASS functional
We begin by introducing the quantity that measures the difference in the free energy on instances where the size differs by one,

$$ \begin{equation} A_j = \mathbb{E}[\log(Z_{j+1})] - \mathbb{E}[\log(Z_j)]\, ,\, \forall j \in \{0,\dots,n-1\}\, . \end{equation} $$

This immediately yields the following observation,

$$ \begin{equation} F_{n, \beta} = \frac{1}{n}\mathbb{E}[\log(Z_n)] = \frac{1}{n}\sum_{i=0}^{n-1}A_i\, .\end{equation} $$

To compute $$A_n $$ we begin by computing $$H_{n + 1} $$ and try to write it as $$H_n $$ plus some other term. This other term will then roughly correspond to the increase in the free energy from the additional of the cavity vertex.

$$ \begin{equation} H_{n+1}(\sigma.\tau) = H'_{n}(\sigma) + \tau\cdot y(\sigma) + \frac{1}{\sqrt{n+1}}g_{n+1, n+1}\end{equation} \, ,$$

where,

$$ \begin{align} H'_{n}(\sigma) &= \frac{1}{\sqrt{n+1}}\sum_{i, j = 1}^n g_{i, j}\sigma_i\sigma_j\, , \\ y(\sigma) &= \frac{1}{\sqrt{n+1}}\sum_{i=1}^n(g_{i, n+1} + g_{n+1, i})\sigma_i \, ,\\ &\tau \in \{\pm 1\} \, . \end{align} $$

Note that, upto a different normalization factor, the expression above for $$H_{n+1}(\sigma) $$ establishes the desired recursive formulation. In order to write _directly_ in terms of the hamiltonian $$H_n(\sigma) $$ (as opposed to $$H'_n(\sigma) $$) we use two observations:
  * The sum of two independent gaussians is another gaussian.
  * The covariance of a centered gaussian process can characterize it.

To begin, note that,
$$ \begin{equation}  \mathbb{E}[H'_n(\sigma_1)H'_n(\sigma_2)] = \frac{n^2}{n + 1}\left(\frac{\langle \sigma_1, \sigma_2 \rangle}{n}\right)^2\, . \end{equation} $$

We'd like it to be the case that when we add another centered gaussian process $$\{x_\sigma\} $$ that is normalized appropriately, this covariance will be become _exactly_ $$n $$. So, we need a process with covariance $$\frac{n}{n + 1} $$ times the overlap on $$n $$ vertices that is _independent_ of the gaussians in $H'_n$. Note that such a process can be given by,

$$\begin{equation} x(\sigma) = \frac{1}{\sqrt{n(n + 1)}}\sum_{i, j =1}^ng'_{i, j}\sigma_i\sigma_j\, ,\end{equation} $$

where $$\{g'_{i,j}\} $$ are independent of the gaussian interactions in $$H'_n(\sigma) $$. Therefore, in distribution,

$$ \begin{equation} H_n(\sigma) \overset{d}{=} H'_n(\sigma) + x(\sigma)\, , \end{equation} $$

which yields the final following cavity equation for the _Aizenman-Sims-Starr_ functional,

$$ \begin{align} A_n &= \mathbb{E}[\log(Z_{n+1})] - \mathbb{E}[\log(Z_n)] \\ &= \mathbb{E}\left[\log\left(\sum_{\sigma, \tau}e^{\beta (H'_n(\sigma) + \tau\cdot y(\sigma) + (n + 1)^{-1/2}g_{n+1, n+1})}\right)\right] - \mathbb{E}\left[\log\left(\sum_{\sigma}e^{\beta(H'_n(\sigma) + x(\sigma))}\right)\right] \\ &= \mathbb{E}\left[\log\left(\left\langle \sum_{\tau}e^{\beta \tau\cdot y(\sigma)}\right\rangle'\right)\right] - \mathbb{E}\left[\log\left(\langle e^{\beta x(\sigma)}\rangle'\right)\right] \\ &= \mathbb{E}\left[\log\left(\left\langle 2\cosh(\beta y(\sigma)) \right\rangle'\right)\right] - \mathbb{E}\left[\log\left(\langle e^{\beta x(\sigma)}\rangle'\right)\right]\, ,\end{align} $$

where we used the fact that the self-interacting gaussian term averages to $$0 $$ as the gaussians in consideration are centered, and $$\langle . \rangle' $$ denotes a Gibbs measure with respect to the re-scaled hamiltonian $$H'_n(\sigma) $$ on $$n $$ vertices. A few comments about the functional $$A_n $$ are in order:

* The functional depends on the average with respect to the Gibbs measure $$\langle . \rangle' $$ defined on the (slightly rescaled) SK hamiltonian on $$n $$ vertices.
* The quantity with respect to which the gibbs measure is made has gaussians that are completely _independent_ of the gaussians in $$x(\sigma) $$ and $$y(\sigma) $$.
* The gaussian processes $$x(\sigma) $$ and $$y(\sigma) $$ are defined on $$n $$ vertices.
* The decomposition works because a _common_ quantity on $$n $$ vertices ($$H'_n $$) is found, which allows us to write the free energy on $$n $$ and $$n + 1$$ vertices as the shared process _plus_ some _independent_ gaussian processes.

The last point is critical: Since we can decouple the hamiltonians on both vertices as a shared gaussian process plus two independent processes, we can actually _average_ (with respect to the Gibbs measure) over the contribution that comes from the cavity vertex, keeping the same iterated average over the Gibbs from both free energy terms that contribute to $$A_n $$.  <br />

### Invariance Symmetries
<br />

## Ruelle Probability Cascades
The Ruelle Probability Cascades are a random measure on a separable Hilbert Space. To motivate this structure, it is critical to introduce a statement of the _Parisi Variational Principle_. There are many equivalent statements of the principle, but we will focus on the one that helps motivate the [Rulle Probability Cascades](#ruelle-probability-cascades). The Parisi Variational Principle essentially gives a variational representation of the optimization problem which constitutes finding the free-energy density of the SK model at _any temperature_ $$\beta $$ - In this sense (among others) the Parisi Variational Principle is stronger than the problem we are interested in since we merely need the $$\beta \to \infty $$ limit, as that specifically that corresponds to the $$\mathsf{MAX} $$-$$\mathsf{CUT} $$ problem over the complete graph with i.i.d. standard gaussian weights. With this definition in mind, we will introduce the construction of the Ruelle Probability Cascades. In a future post, we will then introduce the **Guerra RSB bound** to show that one can upper bound the limit of $$H_n $$ as $$n \to \infty $$ with the _Parisi Variational Principle_[^4]. In doing this, the Ruelle Probability Cascades (RPCs) will play a critical role since they will provide an alternative representation of the Parisi Variational Principle.

### The Parisi Variational Principle
We specify the Parisi Variational Principle after fixing two sequences of parameters that generate the "overlaps" (a sequence $$\{q_i\}^r_{i=0} $$) in the support of a distribution $$\xi $$ along with the cumulative density associated with them given by (though an abuse of notation) by a sequence $$\{\xi_i\}_{i=-1}^r $$. Note that the overlaps are between 0 and 1, so this information can be uniquely by identified by a probability distribution $$\mathcal{D} $$ over the interval $$[0, 1] $$. So, we formally specify this distribution as,
$$ \begin{equation} 0 = \xi_{-1} < \xi_0 < \dots < \xi_{r-1} < \xi_r = 1\, , \end{equation} $$  
and
$$ \begin{equation} 0 = q_0 < q_1 < \dots < q_{r-1} < q_r = 1\, \end{equation} $$
with the constraint that $$\xi(\{q_a\}) = \xi_{a} - \xi_{a-1} $$ enforcing that we are thinking of $$\xi $$ as a CDF. <br />
We begin by introducing an interesting sequence of random variables where,
* A base random variable is defined by a "smoothed" absolute value given to the increase in overlap between two parameters $$q_a $$ and $$q_{a+1} $$ weighted by some standard gaussian variable $$z_{a} $$.
* The rest of the variables are defined recurively as the average (over a bunch of independent standard gaussians) of a $$\xi $$-weighted sum of prior random variables.

It is not easy to motivate why one would define such a sequence of random variables without some knowledge of the _Replica Method_ and the introduction of the _Replica-Symmetry Breaking_ ansatz that ensues to evaluate a certain modified form of a moment of the free energy of the model under consideration assuming some form of "weights" on overlaps between solutions independent copies of the model.

<u><strong>(Parisi Sequence)</strong></u>: Given i.i.d. standard gaussian random variables $$\{z_a\}_{a=1}^r $$, the _Parisi Sequence_ is a sequence of random variables $$\{X^{\xi}_a\}_{a=1}^r $$ defined recursively as,

$$ \begin{align} X^{\xi}_r &= \log\left(2\cosh\left(\sum_{1 \leq a \leq r}\sqrt{2}\beta(q_a - q_{a-1})^{1/2}z_a \right)\right)\, ,\\ X^{\xi}_b &= \frac{1}{\xi_b}\log\left(\mathbb{E}_{z_{b+1}}\left[\exp(\xi_{b+1}X^{\xi}_{b+1})\right]\right)\, ,\, 0 \leq b \leq r-1\, .\end{align} $$

The Parisi-Variational Principle is given by taking the last (or first, depending on how one views the sequence) term of the _Parisi Sequence_ $$X^{\xi}_0 $$ minus the temperature normalized variance of the overlap distribution. The _Parisi Sequence_ is obtained as a solution of the so-called Parisi PDEs using the infamous Hopf-Cole Transformation.

<u><strong>(Parisi Variational Principle)</strong></u>: The Parisi Variational Principle is an optimization problem over the space of distributions with support $$[0, 1] $$ and is given for all $$\beta > 0 $$ as,

$$ \begin{equation} \lim_{n \to \infty} F_{n, \beta} = \inf_{\xi \in \mathcal{D}[0, 1]}\left(X^{\xi}_0 - \beta^2\int_{0}^1t\xi(t)dt\right)\end{equation}\, , $$

where $$\xi(t) $$ is being interpreted as the CDf of the distribution over the overlaps $$\{q_a\} $$.

The Parisi Variational Principle is claiming that at any finite temperature, the free energy denesity of the SK model in the thermodynamic limit can be expressed as a variational optimization problem over the space of distributions with support $$[0, 1] $$ which minimize a so-called "free entropy" term $$X^{\xi}_0 $$ obtained by solving a PDE and the "variance" of the overlaps induced by the solutions. <br />

### The RPC Construction
<br />

#### FOOTNOTES

[^1]: There are many references to this body of work which will be given in due course. Nonetheless, the following two surveys are a nice start: [\[B05\]](http://www.numdam.org/item/SB_2004-2005__47__349_0.pdf), [\[G21\]](https://arxiv.org/pdf/2109.14409.pdf).

[^2]: Stein's Lemma can be proved by a simple integration-by-parts argument for every coordinate $$i $$ in conjunction with Fubini's theorem and the chain-rule. The lemma itself simply asserts that the correlation between a gaussian variable and some function of it is equivalent to a sum of scalings (by the correlations) of the average rate of change of the function. Likewise, Fekete's Lemma boils down to algebraic manipulation of the sequence, where we take the limit infimum of $$x_n/n $$ and compare it to limit infimum of some $$x_m $$ where $$m $$ is chosen to be the supremum as a divisor with a remainder term. The lemma itself merely asserts that for an appropriately growing sequence, the empirical average in the limit simply picks out the largest contributing term.

[^3]: In [this]() upcoming post, various results from the so-called "Gaussian Toolbox" will be stated and briefly proved. This is a very useful set of techniques to have command over in order to prove properties about mean-field spin glasses, and relate them to the behavior of random instances of sparse CSPs.

[^4]: In [this]() upcoming post, we introduce the Guerra-RSB bound and prove that it can be used to show that the Parisi Variaional Principle, represented as an appropriately parameterized RPC, can be used to upper bound the free energy density of the SK model. We will also introduce the famed _Ghirlanda-Guerra_ identities.

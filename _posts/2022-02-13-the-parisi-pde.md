---
layout: post
title:  "The Parisi PDE"
date:   2022-02-12 11:15:40
blurb: "Gaussian Interpolation, The Heat Equation & Hopf-Cole Transformation"
og_image: /assets/img/content/post-example/Banner.jpg
---

<br />
The goal of this post is to show how the Parisi PDE can be solved _completely explicitly_ using the Hopf-Cole transform, which linearizes the PDE. In doing this, many interesting insights are uncovered, including a form of the Parisi formula that hints at a connection to the _Ruelle Probability Cascades_ due to the representation seen previously via the Aizenman-Sims-Starr scheme. In fact, the solution of the Parisi PDE by using the Hopf-Cole transform in conjunction with Gaussian interpolation (that solves the heat equation) and the continuity argument of Guerra hints at a completely explicit representation that turns out to be equivalent to averaging certain gaussian processes over the Ruelle Probability Cascades. This equivalence is the backbone of, both, the proof of the upper bound by Guerra (using the RSB interpolation scheme) and the proof of the lower bound by Panchenko (using the Dovbysh-Sudakov exchangability theorem and the fact that the RPCs satisfy the Ghirlanda-Guerra identities).

<br />


#### Table of Contents
1. [Gaussian Distributions](#gaussian-distributions)
   * [Gaussian Integration by Parts](#gaussian-integration-by-parts)
   * [Gaussian Integration with Gibbs Averages](#gaussian-integration-with-gibbs-averages)
   * [Gaussian Interpolation](#gaussian-interpolation)
2. [The Heat Equation](#the-heat-equation)
   * [Separation of Variables](#separation-of-variables)
   * [A solution via Gaussian Interpolation](#a-solution-via-gaussian-interpolation)
3. [Hopf-Cole Transformation](#hopf-cole-transformation)
   * [The Parisi PDE](#the-parisi-pde)
   * [Piece-wise Hopf-Cole Variables](#piece-wise-hopf-cole-variables)
   * [The Parisi Random Variables](#the-parisi-random-variables)
4. [Footnotes](#footnotes)

<br />

## Gaussian Distributions
We are interested in various analytic properties of Gaussian distributions. While the multivariate gaussian distribution (written $$\mathcal{N}(\mu, \sigma^2)^{\otimes n}$$) is a typical instantiation from this family, it will not be general enough for the setting of Spin Glasses. The reason, as we saw in the [previous post](), is that the Gaussian process associated with the random optimization problems of interest has a subtle covariance structure which depends on the overlap between two configurations. We will, therefore, study the analytic properties of _jointly_ gaussian distributions. To do so, we first introduce some standard notation and important definitions.

**<u>[Jointly Gaussian Process]:</u>** A jointly gaussian process $$\{g_\sigma\}_{\sigma \in \Sigma}$$ will be thought of as a vector of gaussian random variables $$(g_{\sigma_1},\dots,g_{\sigma_{|\Sigma|}})$$ indexed by $$\Sigma$$ (an index set) specified by the mean vector and covariance matrix:
$$\begin{equation}
  \mu_g = (\mu(g_{\sigma_1}),\dots,\mu(g_{\sigma_{|\Sigma|}}))\, .
\end{equation}$$
$$\begin{equation}
  \textsf{Cov}(g_{\sigma}, g_{\sigma'}) = \mathbb{E}[g_{\sigma}g_{\sigma'}] - \mu(g_{\sigma})\mu(g_{\sigma'})\, .
\end{equation}$$

Often, for the purposes of Spin Glasses, the gaussians will be centered ($$\mu(g_{\sigma_i}) = 0 $$) and the index set will simply be the set of all valid spin configurations. Typically, this set will also happen to be convex ($$\{\pm 1\}^n$$ or $$\mathcal{S}^{n-1}(\sqrt{n})$$).

A jointly gaussian process is basically a collection of gaussian random variables with some covariance structure that captures their interdependencies. We typically want order to understand _functions_ of these jointly gaussian processes, since we are interested in optimizing the value of these functions over the underlying index set $$\Sigma $$. To approach this goal analytically, two techniques that are often helpful are _Gaussian Integration by Parts_ and _Gaussian Interpolation_.

### Gaussian Integration by Parts
Gaussian integration by parts allows one to understand the correlation between a function of gaussian variables and the gaussian itself. It turns out that if the function is smooth with reasonable growth properties, the correlation between the gaussian and a function that acts on it is simply related to the average rate of change of the function - This observation is codified as _Stein's Lemma_.

**<u>[Stein's Lemma]</u>:** Given a gaussian random variable $$g \sim \mathcal{N}(0, \sigma^2) $$ and a $$C^1 $$ function $$f: \mathbb{R} \to \mathbb{R} $$ that satisfies mild growth conditions, the following holds,
$$\begin{equation}
  \mathbb{E}[g f(g)] = \sigma^2\mathbb{E}[f'(g)]\, .
\end{equation} $$
_Proof:_ The proof for this fact as actually quite elementary and follows by an application of integration by parts coupled with the observation that,
$$\begin{equation}
  xe^{-x^2/2\sigma^2} = -\sigma^2\partial_x(e^{-x^2/2\sigma^2})\, .
\end{equation}$$
It is then elementary to compute,
$$\begin{align}
  \mathbb{E}[gf(g)] &= a\int x f(x)e^{-x^2/2\sigma^2}dx = -a\sigma^2 f(x)e^{-x^2/2\sigma^2}\vert^{+\infty}_{-\infty} + a\sigma^2\int f'(x)e^{-x^2/2\sigma^2}dx \\
  &= a\sigma^2\int f'(x)e^{-x^2/2\sigma^2} = \mathbb{E}[f'(g)]\, ,
\end{align} $$
where $$a = \frac{1}{\sqrt{2\pi\sigma^2}} $$ and we use the fact that $$f $$ has slower than exponential growth.

As it turns out, by an elementary use of Fubini's theorem and a simple orthogonalization argument (making $$g _i $$ independent of the other gaussian vectors in the process), one can extend Stein's Lemma to hold for an arbitrary jointly gaussian process.

**<u>[High-Dimensional Stein's Lemma]</u>:** Given a jointly gaussian process $$\{g _i\}_{i \in [n]} $$ with mean 0 and a covariance matrix $$\textsf{Cov}_{ij} = \mathbb{E}[g_{i}g_j] $$, and a smooth function $$f: \mathbb{R}^{n} \to \mathbb{R} $$ that doesn't grow too fast, the following holds:
$$\begin{equation}
  \mathbb{E}[g_if(g_1,\dots,g_n)] = \sum_{j=1}^n\textsf{Cov}_{ij}\mathbb{E}[\partial_jf(g_1,\dots,g_n)]\, .
\end{equation}$$
_Proof Sketch_: We provide here only a sketch, and the details are left to the reader as an exercise. We first orthogonalize with respect to $$g_i $$ and consider a "shifted" gaussian process $$\{g'_i\} $$:
$$\begin{equation}
    g'_j = g_j - \frac{\textsf{Cov}_{ij}}{\textsf{Var}_i}g_i\, .
\end{equation}$$
Note, then, that holding $$\{g'_i\} $$ fixed we can compute the expected correlation between $$g_i $$ and $$f(g_1,\dots,g_n) $$ which yields,
$$
\begin{equation}
  \mathbb{E}_{g_i}[g_if(g_1,\dots,g_n)] = \textsf{Var}_i\cdot\mathbb{E}\left[\partial_if\left(g_1 + \frac{\textsf{Cov}_{i1}}{\textsf{Var}_i}x,\dots,g_n + \frac{\textsf{Cov}_{in}}{\textsf{Var}_i}x\right)\right]\, ,
\end{equation} $$
and this follows by an application of the standard Stein's Lemma. To complete the argument, one also computes the expectation over $$\{g'_i\} $$ and uses Fubini's theorem followed by the chain rule, which complete the argument.

### Gaussian Integration with Gibbs Averages

### Gaussian Interpolation

## The Heat Equation

### Separation of Variables

### A solution via Gaussian Interpolation

## Hopf-Cole Transformation
We will now intrdouce a simple change of variables that will allow us to go between non-linear Hamilton-Jacobi equations (such as the Parisi PDE) and the standard heat equation - This change of variables is termed as the **Hopf-Cole Transform**. As it will turn out, it is not sufficient to do this just _once_, but rather this needs to be done on a set of intervals on which the support of the Parisi measure is defined. Nonetheless, courtesy a continuity argument by Guerra, it is possible to approximate the Parisi formula to arbitary precision (in Total Variation Distance) using step-functions on a sufficiently large number of intervals as the supoprt of the CDF, thereby allowing us to rewrite the Free-Entropy term $$\phi_\zeta(x, t) $$ in the Parisi formula to explicitly depend on the typical sequence of parameters:
$$
\begin{equation}
m_0 = 0 < m_1 < \dots < m_{k-1} < m_k = 1\, ,
\end{equation} $$  
$$
\begin{equation}
q_0 = 0 < q_1 < \dots < q_{k} < q_{k+1} = 1\, ,
\end{equation}
$$
which will explicitly specify the step function $$\zeta $$. We will denote the set of all non-decreasing step functions with input and output in $$[0,1] $$ by $$\mathcal{M}_{D}[0,1] $$ and the set of _all_ cumulative distribution functions with support in $$[0, 1] $$ by $$\mathcal{M}[0,1] $$.

### The Parisi PDE
We begin by first writing the Parisi PDE, which can be viewed as a Hamilton-Jacobi equation _when_ restricted to intervals $$[q_l, q_{l+1}) $$. The Parisi PDE specified below is specifically for the SK Model - For the general case of the _mixed_ $$p $$-spin models, there is a dependence on the so-called mixture polynomial $$\xi(.) $$ as well, but we will not concern ourselves with the general case.
$$\begin{equation}
  \partial_t\phi_\zeta(x,t) + \frac{\beta^2}{2}\left(\partial_{xx}\phi_\zeta(x,t) + \zeta(t)\left(\partial_x\phi_\zeta(x,t) \right)^2\right) = 0\, .
\end{equation}
$$
The initial condition to the equation is given at $$t = 1 $$ and the equation is to be solved "backwards in time",
$$\begin{equation}
  \phi_\zeta(x, 1) = \log(\cosh(\beta x))\, .
\end{equation}$$
If we restrict, however, to some interval $$[q_l, q_{l+1}) $$ (which are specificed by some $$\zeta \in \mathcal{M}_D[0,1] $$) with $$\zeta(t) = m_l $$ for $$t \in [q_l, q_{l+1}) $$, we get to solve a sequence of such PDEs,
$$\begin{equation}
  \left\{\partial_t\phi_\zeta(x, t) + \frac{\beta^2}{2}\left(\partial_{xx}\phi_\zeta(x,t) + m_l(\partial_x\phi_\zeta(x,t))^2\right) = 0\right\}_{l \in [k]}\, ,\,\forall t \in [q_l, q_{l+1}]\, .
\end{equation}$$
As it turns out, the continuity argument by Guerra essentially says that we can use a _discrete_ set of intervals to arbitrarily approximate the optimum value of the Parisi formula over the set of all possible measures. Notice, however, this argument doesn't illuminate any structural property about these intervals. For instance, the approximation could be to use two disjoint intervals that cover the entire line, or the have a countable union of infinitesmal ones, or to have a finite number of decreasing-length intervals followed by another self-limiting sequence. The point is that there is a lot of ways to choose the underlying _support_ of these non-decreasing step functions, and that is what allows one to approximate the _true_ Parisi measure arbitrarily well (in Total-Variation Distance) with these.   

### Piece-wise Hopf-Cole Variables
We will now _linearize_ the clearly non-linear Parisi PDE by applying a wonderful trick that dates back to the 1950s, invented independently by Hopf & Cole to solve the visocsity equation. The transformation is as simple as defining the following function,
$$
\begin{equation}
  g(x,t) = e^{c\cdot\phi_\zeta(x,t)}\, , c \in \mathbb{R}\, .
\end{equation}
$$
Another simple thing to do is to _invert the time_ (do a variable substitution with $$t' = (q_{l+1} - t) $$ with $$t \in [q_l, q_{l+1}] $$). The final substitution then yields,
$$
\begin{equation}
  \hat{g}(x,t) = e^{c\cdot\phi(x, q_{l+1} - t)}\, ,\, \forall l \in [k]\, ,
\end{equation}
$$
with the initial condition at $$t = q_{l} $$. We start solving the equations from the _last_ interval, and recurse backwards.

We now substitute the above in the Parisi PDE and some algebra along with the chain-rule convinces us that the Parisi PDE reduces to the following heat equation on _every_ interval $$[q_l, q_{l+1}) $$,
$$
\begin{equation}
  \partial_t\hat{g}(x,t) = \frac{1}{2}\partial^2_x\hat{g}(x,t)\, .
\end{equation}
$$
Note that for the last (and correspondingly first interval solved) this initial conditional simply turns out to be,
$$
\begin{equation}
  \hat{g}(x, 0) = e^{m_l\log(\cosh(\beta x))}\, ,
\end{equation}
$$
where $$l \in [k] $$ and the initial condition will be weighted by the mass $$m_l $$ for each heat equation depending on the interval (and corresponding heat equation) being solved.
The initial conditions for the later intervals can be calculated by repeating the same calculation. However, more critically, note that for _any_ $$t' \in [0, q_{l+1} - q_l] $$, we can compute the solution for any $$t' $$ by using the solution to the heat equation given by _gaussian interpolation_,
$$
\begin{equation}
  \hat{g}(x,t) = \mathbb{E}_{g\sim\mathcal{N}(0,1)}\left[\hat{g}(x + \sqrt{t}g, 0)\right]\, .
\end{equation}
$$
At this point we are done, as taking the inverse Hopf-Cole transform will bring us to the desired function we originally wanted to solve for (which is $$\phi_\zeta(x, t) $$):
$$
\begin{equation}
\phi(x, q_{l+1} - t) = \frac{1}{c}\log(\hat{g}(x,t)) = \frac{1}{m_l}\log\left(\mathbb{E}_{g\sim\mathcal{N}(0,1)}\left[\hat{g}(x + \sqrt{t}g, 0)\right]\right)\, .
\end{equation}
$$
Making the final substitution will give the desired form as follows,
$$
\begin{equation}
  \phi(x, t) = \frac{1}{m_l}\log\left(\mathbb{E}_{g\sim\mathcal{N}(0,1)}\left[\hat{g}(x + \sqrt{q_{l+1} - t}\cdot g, 0)\right]\right)\, , \forall t \in [q_l, q_{l+1}]\, .
\end{equation}
$$
The above expression for $$\phi(x,t) $$ defined piece-wise on each interval with $$t \in [q_l, q_{l+1}) $$ for $$l \in [k] $$ gives the so-called "Piece-wise Hopf-Cole Variables".

### The Parisi Random Variables
These are random variables that recursively encode, on intervals between every interval in the support of the Parisi Measure, the solution to the Heat Equation (after linearization). They use the fact that the solution to the heat equation is known to be an averaging over an appropriately scaled gaussian, and then apply this recursively.

**<u>Parisi Sequence</u>:** Given the functional order parameter $$\zeta $$ encoded above in the sequence $$(m_i, q_j)_{i\in[k],j\in[k+1]} $$ and independent mean 0 gaussian random variables $$(z_j)_{j \in k} $$ with variance,
$$
\begin{equation}
  \mathbb{E}[z^2_j] = \sqrt{2\beta^2(q_{j+1} - q_j)}\, ,
\end{equation}
$$
the Parisi Sequence is defined recursively as follows,
$$
\begin{align*}
  X_{k+1} = \log\left(\cosh\left(\sum_{j=0}^kz_j\right)\right)\, , \\
  X_i = \frac{1}{m_i}\log\left(\mathbb{E}_{z_{i+1}}\exp\left(m_iX_{i+1}\right)\right)\, , 0 \leq i \leq k\, .
\end{align*}
$$

Using the above Parisi sequence, the Parisi functional can be written in a very straightforward way (in a form that is different from the one that Parisi originally gave):

**<u>Parisi Functional [Cascade Representation]:</u>** The Parisi functional over $$\zeta \in \mathcal{M}_D[0,1] $$ can be represented using the Parisi Sequence as follows,
$$
\begin{equation}
  \mathcal{P}_\beta(\zeta) = X_0 - \frac{\beta^2}{2}\int_0^1t\zeta(t)dt\, .
\end{equation}
$$

 It suffices, due to the continuity argument of Guerra, to optimize the above representation over discrete step-functions to approximate (to arbitrary precision) the original Parisi formulation using the Parisi PDE over continuous measures.

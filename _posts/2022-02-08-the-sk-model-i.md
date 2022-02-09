---
layout: post
title:  "The SK Model I"
date:   2022-02-08 11:15:40
blurb: "The Aizenman-Sims-Starr Representation & Ruelle Probability Cascades"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)
<br />
The goal of this post is to introduce the _Sherrington-Kirkpatrick_ model, which is the canonical starting point in Spin-Glass Theory. This model (along with a slightly generalized family) is one of the first models for which the infamous _Parisi-Variational Principle_ was (formally) proven to be true. It is also a model that instigated the initiation of many other concepts via the _Replica-Symmetry Breaking_ ansatzen such as _Ultrametricity_, the _TAP Equations_, and the _Ghirlanda-Guerra Identities_. It is the prototypical model of a **Mean-Field Spin Glass** and has recently had algorithmic implications, which in conjunction with the growing body of work on the _Overlap-Gap Property_ have created an insteresting line of work that shed lights on the average-case complexity of a large family of optimization problems[^1].

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
Given $$n^2$$ i.i.d. standard normal ($$\mathcal{N}(0, 1)$$) variables $$\{J_{ij}\}_{i, j \in [n]}$$, we are interested in the optimal value of the following optimization problem over the hypercube,

$$ \begin{equation} \max_{\sigma \in \{\pm 1\}^n} \frac{1}{\sqrt{n}}\sum_{i, j=1}^n J_{i,j}\sigma_i\sigma_j\, . \end{equation}$$

The problem above can be seen as asking for the $$\mathsf{MAX}$$-$$\mathsf{CUT}$$ of a complete graph with i.i.d. $$\mathcal{N}(0,1)$$ weights. As mentioned, we will be interested in the _typical_ optimal value of the above problem as $$n$$ gets large. The normalized quadratic form above will be expressed as $$H_n(\sigma)$$ which, in physics language, translates to asking for the value of the hamiltonian $$H_n$$ under the configuration of spins $$\sigma \in \{\pm 1\}^n$$. Therefore, we will ask the following question instead,

$$ \begin{equation} \lim_{n \to \infty} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma) \right]\, . \end{equation}$$

While this problem seems well defined, we need to nornalize it carefully so that the limit we wish to compute does not diverge.

<br />

### Covariance and Overlaps
A simple observation (formalized below) will reveal that the fluctuations (variance) of the term above leads to a divergent limit. Therefore, we must normalize it appropariately. To know what normalization is appropriate, we compute the covariance of the underlying gaussian process explicitly,

$$ \begin{align} \mathbb{E}[H_n(\sigma^1)H_n(\sigma^2)] &= \frac{1}{n}\mathbb{E}\left[ \sum_{i_1, j_1}J_{i_1, j_1}\sigma^1_{i_1}\sigma^2_{j_1}\sum_{i_2, j_2}J_{i_2, j_2}\sigma^1_{i_2}\sigma^2_{j_2}\right] \\ &= \frac{1}{n}\mathbb{E}\left[\sum_{i_1, j_1}J^2_{i_1, j_1}\sigma^1_{i_1}\sigma^2_{i_1}\sigma^1_{j_1}\sigma^2_{j_1}\right] \\ &= \frac{1}{n}(\sum_{i =1}^n \sigma^1_i\sigma^2_i)^2 = n\cdot (\frac{1}{n}\langle \sigma^1, \sigma^2\rangle)^2\, . \end{align} $$

The above calculation reveals that the covariance of the quantity we wish to maximize is a $$O(n)$$ scaling of the _normalized overlap_ between two configurations $$\sigma^1$$ and $$\sigma^2$$. Since the latter quantity is a constant, this tells us that the extra normalization term we wish to add divides the quantity to optimize by $$n$$. This will yield the final average quantity whose limit we wish to compute precisely,

$$ \begin{equation} \lim_{n \to \infty} \frac{1}{n} \mathbb{E}\left[\max_{\sigma \in \{\pm 1\}^n} H_n(\sigma)\right]\, , \end{equation} $$

and this quantity will be termed the _Ground State Energy_ of the system. When studying such quantities, it is a standard trick in statistical physics to understand this quantity using its _smoothed_ approximation, known as the _Free Energy Density_ of the system. This quantity depends on the smoothing parameter $$\beta$$ and is defined as,

$$ \begin{equation} F_{n, \beta} := \frac{1}{n}\mathbb{E}\left[\log(\sum_{\sigma \in \{ \pm 1\}^n} e^{\beta H_n(\sigma)})\right] \, , \end{equation} $$
where the exponential summation term, 

$$ \begin{equation} Z_{n, \beta} := \sum_{\sigma \in \{\pm 1\}^n}e^{\beta H_n(\sigma)}\, , \end{equation} $$

is termed the _Partition function_ in Statistical Physics. As we shall see, this softmax trick induces a measure called the _Gibbs Measure_ which will turn out to be very closely related to the _Free Energy_ of the model and will appear as soon as a derivative of the latter is taken with respect to an appropriate parameter in the hamiltonian. Understanding the asymptotic geometric structure of the gibbs, therefore, will be of paramount importance in formally proving the Parisi Variational Principle.

<br />

### Guerra-Tonnineli Interpolation
<br />

### Gaussian Concentration
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

[^1]: There are many references to this body of work which will be given in due course. Nonetheless, the following two surverys are a nice start: [G21](https://arxiv.org/pdf/2109.14409.pdf), [B05](http://www.numdam.org/item/SB_2004-2005__47__349_0.pdf).

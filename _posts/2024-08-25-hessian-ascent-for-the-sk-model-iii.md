---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model III"
date:   2024-08-25 11:15:40
blurb: "Optimizing the SK model via free probability and convex duality"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)

This is the third (and final) blog post in a **3-part** series, themed on a recent result [[JSS24]](https://arxiv.org/abs/2408.02360) by [David](https://davidjekel.com/), [Jonathan](https://jshi.science/) and me. The first blog post can be found [here](https://juspreetsandhu.me/2024/08/08/hessian-ascent-for-the-sk-model-i), and the second one can be found [here](https://juspreetsandhu.me/2024/08/12/hessian-ascent-for-the-sk-model-ii).

A quick recap of what we accomplished in the second blog post:
1. We proved the first major component in the analysis of the PHA algorithm on the SK model: the construction and desirable properties of the projector $$P_j(D)^2 $$ into the top-eigenspace of the TAP corrected Hessian.
2. We then overviewed the second major component in the analysis of the PHA algorithm on the SK model: the proof that the empirical distribution of the coordinates of every iterate of the algorithm converge, with high probability, to the (primal) Auffinger-Chen SDE.
<br>
<br>
<br>

#### Table of Contents
1. [Connections: Potential Hessian ascent and other algorithms](#connections:-potential-hessian-ascent-and-other-algorithms)
   * [Potential Hessian ascent and high-entropy step processes](#potential-hessian-ascent-and-high-entropy-step-processes)
   * [Potential Hessian ascent and approximate-message passing](#potential-hessian-ascent-and-approximate-message-passing)
2. [Ongoing work](#ongoing-work)
   * [Parisi-type formulae for different geometries](#parisi-type-formulae-for-different-geometries)
   * [Low-degree sum-of-squares certificates for HES processes on the cube](#low-degree-sum-of-squares-certificates-for-hes-processes-on-the-cube)
   * [Potential Hessian ascent: mixed p-spin models](potential-hessian-ascent:-mixed-p-spin-models-and-bounded-degree-csps)
3. [The ultimate goal: a research program](the-ultimate-goal:-a-research-program)
3. [Footnotes](#footnotes)
<br>
<br>
<br>

## Connections: Potential Hessian ascent and other algorithms

The main objective in this section is to discuss the relationship between potential Hessian ascent (PHA) and other families of algorithms. In discussing this, we will come across various natural questions (which are still open problems) involving algorithmic unification and (possibly efficient) certificates for reasoning about the limitations of entire families of algorithms.

While the analysis for the PHA algoithm on the cube is more challenging than that of the corresponding Hessian ascent algorithm of Subag on the sphere [[Sub19]](https://arxiv.org/abs/1812.04588), the algorithmic principle is the same: follow the top-eigenspace of the Hessian of the generalized TAP energy. Consequently, our comparisons will focus on algorithmic paradigms that have found success in related average-case optimization problems but are far from obvious in how they compare with PHA.

### Potential Hessian ascent and high-entropy step processes
In a recent result [[SS24]](https://arxiv.org/abs/2401.14383), Jonathan and I introduced a sum-of-squares hierarchy to _certify_ the suprema of certain Gaussian processes. However, the sum-of-squares certificates (of low-degree) output by this hierarchy (termed the HES SoS hierarchy) do not certify the value of an instance of the problem, as is traditional for certification algorithms. Instead, they certify the value that can be achieved on an instance by a certain family of measures over the solution domain. This family of measures are termed **high-entropy step distributions**.

For the spherical spin glass problem, we demonstrate that the certified value matches (up to constant factors) the _true_ value under a technical condition. Even in the absence of said condition, the certificates obtain the correct _scaling_ of the maximum; due to known lower-bounds against the standard SoS hierarchy for the same problem [[BGL16]]()[[HPK17]](), the standard hierarchy is off by superlinear factors.

So, what relationship do HES processes have with PHA? For starters, in [[Theorem 7.1, SS24]](https://arxiv.org/abs/2401.14383) we demonstrate that a randomized version of Subag's Hessian ascent algorithm is a HES process. The PHA algorithm on the cube _also_ seems to be a HES process[^1]. Importantly, the geometry seems to not affect the description of the underlying algorithmically optimal process. In fact, we think that HES processes will capture the algorithmic threshold of **all** polynomial-time algorithms for a _larger_ class of polynomial optimization problems [[Conjecture 1.6, SS24]](https://arxiv.org/abs/2401.14383).

There are various subtleties in demonstrating that the 

### Potential Hessian ascent and approximate-message passing

## Ongoing work

### Parisi-type formulae for different geometries

### Low-degree sum-of-squares certificates for HES processes on the cube

### Potential Hessian ascent: mixed p-spin models

## The ultimate goal: a research program


#### FOOTNOTES

[^1]: This is a premeditated crime. The construction of the conditional covariances $$P_j(D)^2 $$ of the PHA algorithm, and the properties about it proved in [[Theorem 4.1, JSS24]()], have the dual aim of making it easier to demonstrate convergence to the AC SDE _and_ also satisfy the properties imposed on the conditional covariances of a HES process.

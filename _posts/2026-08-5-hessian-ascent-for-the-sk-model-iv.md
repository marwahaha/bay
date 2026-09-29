---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model IV"
date:   2026-08-05 11:15:40
blurb: "TVD sampling for the SK model via Hessian dynamics, Jarzynski's equality and Entropy contraction"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)


Continuing on the program that [Jonathan](https://www.jshi.science/) and I started with [David](https://davidjekel.com/) in [PHA 1](https://arxiv.org/abs/2408.02360) and [PHA 2](), with [Ewan](https://www.ewandavies.org/) and [Holden](https://holdenlee.github.io/), we put out a PHA 3 preprint on the arXiv [(Potential Hessian Ascent III: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1/2, DLSS26)](https://arxiv.org/abs/2605.03718). Shortly after PHA 3, Holden, Jonathan and I put out a follow-up, creatively termed [(Potential Hessian Ascent IV: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1, LSS26)](https://arxiv.org/abs/2609.30590) that builds on the cavity interpolation theory and free probability toolkit built in PHA 3 to introduce a "local" error analysis which only uses (and proves) "local" regularity properties of the underlying SDEs. 

Togethere, these works together yield a $$o_n(1) $$ TVD sampler for the SK model up to the replica-symmetric threshold $$\beta < 1 $$ which is the conjectured hardness threshold for sampling. The algorithm combines algorithmic stochastic localization (ASL) with rejection sampling over path-space via Jarzynski's equality (JE). The analysis uses local regularity properties of the TAP Hessian and the PHA framework, developing new cavity interpolation theory and an extension to the free-probability toolkit introduced in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360), which is then combined with delicate SDE error analysis and the ability to efficiently sample from sufficiently "localized" distributions, i.e., those with large enough external field induced by the stochastic localization (SL) process. 

I have decided to write, with a guest contribution from Holden, a **5-part** blog post explaining the background, the development of the algorithm, the main proof skeletons, and highlighting the main technical innovations. 

- In this first blog post, I will explain the connection between ASL and PHA, introduce Jarzynski's equality, the over-determined system of ASL, TAP and PHD, and the natural development of the ''desiderata'' that must be proved for the algorithm to succeed. I will end with three consequences (and perspectives) of the result -- one for theoretical physicists, another for analysts/probabilists, and a final one for theoretical computer scientists.   
- In the second blog post, I will spend time developing the cavity interpolation theory that gives exact moment estimates and ''better-than-naive''' overlap-concentration along the SL process. This is necessary to prove the desideratum that the algorithmic covariance (Hessian of the TAP free energy) stays close to the true covariance along the SL path. It will turn out to be the case that this theory holds for $$\beta < 1$$ even though this was not explicitly written down in PHA 3, and a couple of final details for this were not fleshed out till PHA 4. The key technical contribution here is to run a cavity interpolation for the planted SL model with SL tilt and compute "self-stability" estimates on doing Taylor expansions for bulk-deviation quantities (such as overlaps and magnetizations) to get *exact* moment estimates.
-  In the third blog post, I will spend time explaining the extension to the analysis of the non-commutative interpolation developed in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360) to control the diagonal entries of the (squared) algorithmic covariance -- this will then imply various regularity estimates and show that the ASL-TAP process closely tracks the PHD process. We will first do this for $$ \beta < 1/2 $$, and then I will then show that, in fact, one can reason about a "regularized" resolvent which agrees with the original resolvent *only* over a (possibly non-convex) subset of $$(-1,1)^n $$ and, simultaneously, *significantly* simplify the free interpolation analysis in PHA 3 by using a result of [Bandeira, Boedihardjo and van-Handel, 2023](https://arxiv.org/abs/2108.06312). This allows one to reason about the concentration and free limits of the diagonal sub-algebra of the TAP Hessian (a resolvent) all the way to $$\beta < 1 $$ where it *can* be singular in a subset of $$(-1,1)^n $$.   
-  In the fourth blog post, Holden will show how the localized distribution concentrates on a wedge of the hypercube (provided it is run for sufficiently long) and how a 2-stage decomposition on this wedge allows one to show the entropy contraction property, which is the final desideratum that needs to be proved -- this is the approach adopted in PHA 3. After that, I will briefly overview how a result of [Kumar et al, 2026](), which allows one to sample from a Gibbs distribution with quadratic potential and sufficiently strong external field on a wedge of $$\lbrace -1,1\rbrace^n $$ using an annealing-type procedure, is used in PHA 4 to directly sample from the localized distribution after invoking the fact that SL run for large time *will* cause the distribution to be heavilty concentrated inside a wedge of sufficient size.
-  In the final blog post, I will talk about *how* to get access to a "safe set" of points $$S_{A}(c) \subset \lbrace -1,1\rbrace^n $$ when one fixes a "good" event over the input $$A $$ so that, with high probability over the SL paths, the TAP Hessian evaluated at any $$m \in S_{A}(c) $$ is $$c $$-strongly convex. This is the place where, perhaps unsurpringly, we rely on a result about AMP state evolution. Essentially, we combine the fact that posterior means of the SL tilted measures are "stable" in time for most SL paths with the $$c $$-strong convexity of Hessians in small balls around AMP iterates, to build "tubes" around good SL paths that we then union to be contained in $$S_{A}(c) $$. From this, one straightforwardly obtains that, with high probability over $$ A $$ and the SL paths, uniformly in time $$t \in [0,T] $$, small enough balls around the SL-tilted means $$\lbrace m_t\rbrace_{0 \le t \le T} $$ lie in $$S_{A}(c) $$. 
<br>

#### Table of Contents
1. [Algorithm design](#algorithm-design)
   * [Stochastic Localization and Hessian Dynamics](#the-parisi-formula-and-auffinger-chen-representation)
   * [The TAP free energy](#the-generalized-tap-free-energy)
   * [Jarzynski equality](#a-primal-theory-for-the-parisi-pde-via-convex-duality)
2. [ASL, TAP and PHD]()
   * [Overdetermined system and errors](#overdetermined-system-and-error)
   * [Emergent desiderata](#emergent-desiderata)
4. [Footnotes](#footnotes)
<br>

## Algorithm Design
Stuff and things
<br>

### Stochastic Localization and Hessian Dynamics
Stuff and things
<br>

### The TAP free energy
Stuff and things
<br>

### Jarzynski equality
Stuff and things
<br>

## ASL, TAP and PHD
Stuff and things
<br>

### Overdetermined system and errors
Stuff and things
<br>

### Emergent desiderata
Stuff and things
<br>

#### FOOTNOTES

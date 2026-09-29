---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model IV"
date:   2026-08-05 11:15:40
blurb: "TVD sampling for the SK model via Hessian dynamics, Jarzynski's equality and Entropy contraction"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)


Continuing on the program that [Jonathan](https://www.jshi.science/) and I started with [David](https://davidjekel.com/) in [PHA 1]() and [PHA 2](), with [Ewan]() and [Holden]() we put out a preprint on the arXiv [(Potential Hessian Ascent III: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1/2, DLSS26)](https://arxiv.org/pdf/2605.03718). Shortly after PHA-3, Holden, Jonathan and I put out a follow-up, creatively termed [(Potential Hessian Ascent IV: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1, LSS26)]() that builds on the cavity ionterpolation theory and free probability estimates built in PHA-3 to introduce a "local" error analysis which only uses (and proves) "local" regularity properties of the underlying SDEs. 

These works together yield a $$o_n(1) $$ TVD sampler for the SK model up to $$\beta < 1 $$ using an algorithm that combines algorithmic stochastic localization with Jarzynski's equality. The analysis uses the PHA framework, developing new cavity interpolation theory and combining it with novel estimates/techniques from free-probability developed in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360) and the ability to efficiently sample from sufficiently "localized" distributions. I have decided to write, with a guest contribution from Holden, a **5-part** blog post explaining the background, the development of the algorithm, the three main proof skeletons, and highlighting the main technical innovations.

- In this first blog post, I will explain the connection between algorithmic stochastic localization and potential Hessian ascent, Jarzynski's equality, the over-determined system of ASL, TAP and PHD, and the natural development of the ''desiderata'' the analysis must prove. 
- In the second blog post, I will spend time developing the cavity interpolation theory that gives exact moment estimates and ''better-than-naive''' overlap-concentration necessary to prove the desideratum that the algorithmic covariance (Hessian of the TAP free energy) stays close to the true covariance along the SL path. It will turn out to be the case that this theory holds for $$\beta < 1$$ even though this was "implicit" (and a couple of final details for this wre not fleshed out till PHA 4).
-  In the third blog post, I will spend time explaining the extension to the analysis of the non-commutative interpolation developed in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360) to control the diagonal entries of the (squared) algorithmic covariance -- this will then imply various regularity estimates and show that the ASL-TAP process closely tracks the PHD process. I will then show that, in fact, one can reason about a "regularized" resolvent defined only over a (possibly non-convex) subset of $$(-1,1)^n $$ *significantly* simplify the free interpolation analysis in PHA 3 by using a result of [Bandeira, B.. and van-Handel, 2023]().   
-  In the fourth blog post, Holden will show how the localized distribution concentrates on a wedge of the hypercube (provided it is run for sufficiently long) and how a 2-stage decomposition on this wedge allows one to show the entropy contraction property, which is the final desideratum that needs to be proved. After that, I will write about how one can use a result of [Kumar et al, 2026]() which allows one to reason about sampling using an annealing-type procedure from a sufficiently localized distribution for *any* $$\beta > 0 $$ and this result, used in PHA 4, replaces the prior decomposition based on adopted in PHA 3.
-  The "local" regularity       
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

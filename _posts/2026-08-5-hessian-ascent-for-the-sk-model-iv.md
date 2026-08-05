---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model IV"
date:   2026-08-05 11:15:40
blurb: "TVD sampling for the SK model via Hessian dynamics, Jarzynski's equality and Entropy contraction"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)


Continuing on the program that [Jonathan](https://www.jshi.science/) and I started with [David](https://davidjekel.com/) in [PHA 1]() and [PHA 2](), with [Ewan]() and [Holden]() we put out a preprint on the arXiv [(Potential Hessian Ascent III: Sampling the Sherrington-Kirkpatrick Model at ß < 1/2, DLSS26)](https://arxiv.org/pdf/2605.03718).

This work gives a $$o_n(1) $$ TVD sampler for the SK model up to $$\beta < 1/2 $$ using an algorithm that combines algorithmic stochastic localization with Jarzynski's equality. The analysis uses the PHA framework, developing new cavity interpolation theory and combining it with novel estimates/techniques from free-probability developed in []() and entropy contraction from Markov chain theory. I have decided to write, with guest contributions from Holden, a **4-part** blog post explaining the background, three main proof structures, and main technical innovations of this result and the open questions it naturally induces.

In this first blog post, I will explain the connection between algorithmic stochastic localization and potential Hessian ascent, Jarzynski's equality, the over-determined system of ASL, TAP and PHD, and the natural development of the ''desiderata'' the analysis must prove. In the second blog post, I will spend time developing the cavity interpolation theory that gives exact moment estimates and ''better-than-naive''' overlap-concentration necessary to prove the desideratum that the algorithmic covariance (Hessian of the TAP free energy) stays close to the true covariance along the SL path. In the third blog post, I will spend time explaining the extension to the analysis of the non-commutative interpolation developed in [PHA 1, Section-4]() to control the diagonal entries of the (squared) algorithmic covariance as this will then imply various regularity estimates and show that the ASL-TAP process closely tracks the PHD process leading to the next set of desiderata being proved. In the final blog post, I will show how the localized distribution concentrates on a wedge of the hypercube (provided it is run for sufficiently long) and how a 2-stage decomposition on this wedge allows one to show the entropy contraction property, which is the final desideratum that needs to be proved.   
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

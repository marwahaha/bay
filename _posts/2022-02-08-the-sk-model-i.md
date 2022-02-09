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

$$ \begin{equation} \max_{\sigma \in \{\pm 1\}^n} \sum_{i, j=1}^n J_{ij}\sigma_i\sigma_j\, . \end{equation}$$

The problem above can be seen as asking for the $$\mathsf{MAX}$$-$$\mathsf{CUT}$$ of a complete graph with i.i.d. $$\mathcal{N}(0,1)$$ weights.
<br />

### Covariance and Overlaps
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

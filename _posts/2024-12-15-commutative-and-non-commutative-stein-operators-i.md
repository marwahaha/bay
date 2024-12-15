---
layout: post
section: "Probability, Spin Glasses and Complex Analysis"
title: "Commutative and Non-Commutative Stein Operators I"
date: 2024-12-12 11:15:40
blurb: "Stein's Method, Berry-Esseen Theorem and Normal Approximation of Lipschitz Functions"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)
<br />

In this series of two posts, I want to talk about Stein's method. In the first post, the goal will be to introduce an analytic characterization of the normal distribution via Stein's operator. This characterization leads to a set of techniques to show when distributions are "close to" normal. Naturally, this finds uses in getting quantiative variants of the CLT, and the method has relationship with Markov semigroups. In the [second post](), we will focus on a _non-commutative_ version of Stein's lemma, and a technique of Collins, Guionnet and Parraud [[CGP22]]() which allows one to "replace" sufficiently large random matrices with their idealized free-probabilistic analog with vanishing error.
<br>
<br>

#### Table of Contents
1. [Stein's Method](#stein's-method)
   * [The easy direction of Stein's lemma](#the-easy-direction-of-stein's-lemma)
   * [Stein operator for normal variables](#stein-operator-for-normal-variables)
   * [Stein identities for sums and exchangeable pairs of random variables](#stein-identities-for-sums-and-exchangeable-paris-of-random-variables)
2. [Normal approximation of Lipschitz functions](#normal-approximation-for-lipschitz-functions)
3. [Berry-Esseen theorem: CLT with quantitative convergence rates](#berry-esseen-theorem:-clt-with-quantitative-convergence-rates)
4. [Footnotes](#footnotes)
<br>

## Stein's Method

The Stein operator leads to an _exact_ characterization of the Gaussian distribution via an analytic property. Stein's lemma says that a random variable $$Z $$ over $$(\mathbb{R}, \mathcal{B}(\mathbb{R})) $$ is distributed as a Gaussian if-and-only-if the corresponding Stein operator is exactly zero on a sufficiently nice class of test functions. One direction of this proof is quite straightforward and follows by some integration-by-parts and the explicit use of the Gaussian density function. The other direction, however, uses a cleverly constructed test function that solves a particular ODE. Then, on proving various Lipschitz properties about the solution and its derivative using tail-estimates and numerous careful analytic facts, one obtains that the solution is, in fact, a valid test function. 

More generally, this leads to a technique of defining a Stein operator for many distributions. There are two ways to do this, but we will only mention a way which involves relating the generator of a particular Markov process whose stationary distribution is related to the distribution we are interested in. We will touch upon this when we show how to derive "Stein identities" for sums of independent random variables and a coupling for "exchangeable pairs" of random variables. 


In fact, the idea of approximating reasonable random variables by normal variables with quantiative precision is partially summarized in a proposition we will state later. To conclude, we will then see how Stein's method gives a quantitative version of the CLT known as the Berry-Esseen theorem.


### The easy direction of Stein's lemma


### Stein operator for normal variables


### Stein identities for sums and exchangeable pairs of random variables


## Normal approximation of Lipschitz functions


## Berry-Esseen theorem: CLT with quantitative convergence rates


#### FOOTNOTES

[^1]: 

[^2]: 

[^3]: 

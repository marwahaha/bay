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
3. [Berry-Esseen theorem](#berry-esseen-theorem)
4. [Footnotes](#footnotes)
<br>
<br>


## Stein's Method

The Stein operator leads to an _exact_ characterization of the Gaussian distribution via an analytic property. Stein's lemma says that a random variable $$Z $$ over $$(\mathbb{R}, \mathcal{B}(\mathbb{R})) $$ is distributed as a Gaussian if-and-only-if the corresponding Stein operator is exactly zero on a sufficiently nice class of test functions. One direction of this proof is quite straightforward and follows by some integration-by-parts and the explicit use of the Gaussian density function. The other direction, however, uses a cleverly constructed test function that solves a particular ODE. Then, on proving various Lipschitz properties about the solution and its derivative using tail-estimates and numerous careful analytic facts, one obtains that the solution is, in fact, a valid test function. 

More generally, this leads to a technique of defining a Stein operator for many distributions. There are two ways to do this, but we will only mention a way which involves relating the generator of a particular Markov process whose stationary distribution is related to the distribution we are interested in. We will touch upon this when we show how to derive "Stein identities" for sums of independent random variables and a coupling for "exchangeable pairs" of random variables. 

In fact, the idea of approximating reasonable random variables by normal variables with quantiative precision is partially summarized in a proposition we will state later. To conclude, we will then see how Stein's method gives a quantitative version of the CLT known as the Berry-Esseen theorem.

We begin by writing down Stein's lemma, and then spending the next two sections establishing its proof.

_<u>(Stein's Lemma)</u>:_ The following statements are true:
* If $$Z \sim \mathcal{N}(0,1) $$, then for every absolutely continuous function $$f : \mathbb{R} \to \mathbb{R} $$ with $$\mathbb{E}[f'(Z)] < \infty $$, it is the case that

$$ \begin{equation} \mathbb{E}[f'(Z)] = \mathbb{E}[Zf(Z)]\, . \end{equation} $$

* If, for every $$f: \mathbb{R} \to \mathbb{R} $$ that is absolutely continuous, bounded and continuously differentiable, it is true that

$$ \begin{equation} \mathbb{E}[f'(W)] = \mathbb{E}[Wf(W)]\, , \end{equation} $$

then $$W \sim \mathcal{N}(0,1) $$.

### The easy direction of Stein's lemma
The proof for the first direction is simple. We use the fundamental theorem of calculus along with the fact that $$xe^{-x^2/2} $$ is an odd function to evaluate the correlation between the normal variable $$Z $$ and $$f(Z) $$, 

$$ \begin{equation} \mathbb{E}[Zf(Z)] = \frac{1}{\sqrt{2\pi}}\int_x (x e^{-x^2/2})f(x)dx = \int_x \left(\int_0^x f'(y)dy\right)xe^{-x^2/2}dx\,. \end{equation} $$ 

Splitting the integral around $$0 $$ yields

$$ \begin{equation} \mathbb{E}[Z(f(Z))] = \frac{1}{\sqrt{2\pi}}\left(\int_{-\infty}^0\left(\int_x^0 f'(y)dy\right)(-xe^{-x^2/2})dx + \int_{0}^{\infty}\left(\int_0^x f'(y)dy\right)xe^{-x^2/2}dx\right)\,. \end{equation} $$

Rearranging integrals, which is justified by an invocation of Fubini's theorem and the absolutely continuity of $$f $$, gives

$$ \begin{equation} \mathbb{E}[Zf(Z)] = \frac{1}{\sqrt{2\pi}}\left(\int_{-\infty}^0f'(y)\left(\int_{-\infty}^y-xe^{-x^2/2}dx\right)dy + \int_0^{\infty} f'(y)\left(\int_y^\infty xe^{-x^/2}dx\right)dy\right)\,. \end{equation} $$  

At this point, one observes that $$d/dx (e^{-x^2/2}) = -xe^{-x^2/2} $$ and one more invocation of the fundamental theorem of calculus along with elementary limit taking simplifies the RHS of the above equation to

$$ \begin{equation} \mathbb{E}[Zf(Z)] = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^{\infty}f'(y)e^{-y^2/2}dy = \mathbb{E}[f'(Z)]\,. \end{equation} $$

### Stein operator for normal variables
We now prove the second part of Stein's lemma. To do this, we will (somewhat painfully) probe the analytic properties of a parmaterized family of functions $$\{f_{z}(w)\}_{z \in \mathbb{R}} $$ that solve a particular ODE. This ODE measures, for every $$z \in \mathbb{R} $$, the pointwise difference in the action of the Stein operator with the pointwise difference in a function that indicates being $$\le z $$ with the CDF of a normal at $$z $$. We notate this function as

$$ \begin{equation} g_z(w) := \mathbb{1}_{w \le z} - \Phi(z)\, , \end{equation} $$

where $$\Phi(z) := \mathsf{Pr}_{w \sim \mathcal{N}(0,1)}[w \le z] $$.


Note that, for every $$z \in \mathbb{R} $$,

$$ \begin{equation} \mathbb{E}_{w \sim \mathcal{N}(0,1)}[g_z(w)] = \mathbb{E}_{w \sim \mathcal{N}(0,1)}[\mathbb{1}_{w \le z}] - \mathbb{E}_w[\Phi(z)] = \Phi(z) - \Phi(z) = 0\,. \end{equation} $$ 

Therefore, this implies that the average of $$g_z(w) $$ under $$w \sim \mathcal{N}(0,1) $$ is $$0 $$ for every $$z $$. The parameterized families of functions, therefore, will solve an ODE of the Stein operator acting on $$f_z(w) $$ for _every_ $$z $$ in a way that the RHS is $$0 $$ on average. This begs the question: what is the Stein operator?

_<u>(Stein operator)</u>_: The Stein operator $$\mathcal{A} $$ maps every bounded, continuously differentiable function $$f: \mathbb{R} \to \mathbb{R} $$ to a new function denoted as

$$ \begin{equation} \mathcal{A}[f](w) := f'(w) - wf(w)\,. \end{equation} $$

Therefore, in order to "memorize" the standard normal distribution via the Stein operator, we want it to be the case that,

$$ \begin{equation} \mathbb{E}_{w}[\mathcal{A}[f](w)] = 0 = \mathbb{E}_w[g_z(w)]\, . \end{equation} $$

This naturally induces the Stein ODE, to which the parametrized family of functions $$\{f_z(w)\}_{z \in \mathbb{R}} $$ present a solution. The ODE reads

$$ \begin{equation} \mathcal{A}[f_z](w) = \mathbb{1}_{w \le z} - \Phi(z)\,. \end{equation} $$

We are given that $$\mathbb{E}_w[\mathcal{A}[f](w)] = 0 $$ for every "nice enough" function $$f $$. Consequently, proving the second part of Stein's lemma now reduces to showing that the solutions $$\{f_z(w)\}_{z \in \mathbb{R}} $$ to the Stein ODE are bounded, continuous, and piecewise differentiable. If so, then we immediately have that

$$ \begin{equation} \mathbb{E}_{w}[\mathcal{A}[f](w)] = 0 = \mathsf{Pr}_w[w \le z] - \Phi(z)\, , \end{equation} $$

which implies that $$w $$ is distributed as standard normal.

We will first give an explicit expression for $$f_z(w) $$ when it uniquely solves the Stein ODE, and then we will show that this solution is bounded, continuous and piecewise continuously differentiable.

**<u>Unique solution to the Stein ODE</u>**:  


**<u>Boundedness & differentiability of the solution</u>**:

### Stein identities for sums and exchangeable pairs of random variables


## Normal approximation of Lipschitz functions


## Berry-Esseen theorem


#### FOOTNOTES

[^1]: 

[^2]: 

[^3]: 

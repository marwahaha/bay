---
layout: post
section: "Paper Expositions"
title:  "Hessian Ascent for the SK Model IV"
date:   2026-08-05 11:15:40
blurb: "TVD sampling for the SK model via Hessian dynamics, Jarzynski's equality and localized sampling"
og_image: /assets/img/content/post-example/Banner.jpg
---

[//]: # (<img src="{{ "/assets/img/content/post-example/Banner.jpg" | absolute_url }}" alt="bay" class="post-pic"/>)


Continuing on the program that [Jonathan](https://www.jshi.science/) and I started with [David](https://davidjekel.com/) in [PHA 1](https://arxiv.org/abs/2408.02360) and [PHA 2](), with [Ewan](https://www.ewandavies.org/) and [Holden](https://holdenlee.github.io/), we put out a PHA 3 preprint on the arXiv [(Potential Hessian Ascent III: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1/2, DLSS26)](https://arxiv.org/abs/2605.03718). Shortly after PHA 3, Holden, Jonathan and I put out a follow-up, creatively termed [(Potential Hessian Ascent IV: Sampling the Sherrington-Kirkpatrick Model at $$\beta $$ < 1, LSS26)](https://arxiv.org/abs/2609.30590) that builds on the cavity interpolation theory and free probability toolkit developed in PHA 3 to introduce a "local" SDE error analysis which only uses (and proves) "local" regularity properties for the driver and drift terms underlying the matricial functions of certain quantities that drive the SDEs. 

Together, these works yield a $$o_n(1) $$ TVD sampler for the SK model up to the replica-symmetric threshold $$\beta < 1 $$, above which sampling from the Gibbs measured is conjectured to be hard. The algorithm combines algorithmic stochastic localization (ASL) with rejection sampling over path-space via Jarzynski's equality (JE). The analysis uses local regularity properties of the TAP Hessian and the PHA framework, developing new cavity interpolation theory and an extension to the free-probability toolkit introduced in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360), which is then combined with delicate SDE error analysis and the ability to efficiently sample from sufficiently "localized" distributions, i.e., those with large enough external field induced by a stochastic localization (SL) process that's been run long enough. 

I have decided to write, with a guest contribution from Holden, a **5-part** blog post explaining the conceptual background and development of the algorithm, the main proof skeletons, and the main technical innovations. 

- In this first blog post, I will explain the connection between SL and Hessian dynamics (HD), introduce Jarzynski's equality, discuss the over-determined system of ASL, TAP and PHD SDEs, and the natural development of the ''desiderata'' that must be proved for the algorithm to succeed. I will end with three consequences (and perspectives) of the result -- one for theoretical physicists, another for analysts/probabilists, and a final one for theoretical computer scientists.   
- In the second blog post, I will spend time developing the cavity interpolation theory that gives exact moment estimates and ''better-than-naive''' overlap-concentration along the SL process. This is necessary to prove the desideratum that the algorithmic covariance (Hessian of the TAP free energy) stays close to the true covariance along the SL path. It will turn out to be the case that this theory holds for $$\beta < 1$$ even though this was not explicitly written down in PHA 3, and a couple of final details for this were not fleshed out till PHA 4. The key technical contribution here is to run a cavity interpolation for the planted SL model with SL tilt and compute "self-stability" estimates on doing Taylor expansions for bulk-deviation quantities (such as overlaps and magnetizations) to get *exact* moment estimates.
-  In the third blog post, I will spend time explaining the extension to the analysis of the non-commutative interpolation developed in [PHA 1, Section-4](https://arxiv.org/pdf/2408.02360) to control the diagonal entries of the (squared) algorithmic covariance -- this will then imply various regularity estimates and show that the ASL-TAP process closely tracks the PHD process. We will first do this for $$ \beta < 1/2 $$, and then I will then show that, in fact, one can reason about a "regularized" resolvent which agrees with the original resolvent *only* over a (possibly non-convex) subset of $$(-1,1)^n $$ and, simultaneously, *significantly* simplify the free interpolation analysis in PHA 3 by using a result of [Bandeira, Boedihardjo and van-Handel, 2023](https://arxiv.org/abs/2108.06312). This allows one to reason about the concentration and free limits of the diagonal sub-algebra of the TAP Hessian (a resolvent) all the way to $$\beta < 1 $$ where it *can* be singular in a subset of $$(-1,1)^n $$.   
-  In the fourth blog post, Holden will show how the localized distribution concentrates on a wedge of the hypercube (provided it is run for sufficiently long) and how a 2-stage decomposition on this wedge allows one to show the entropy contraction property, which is the final desideratum that needs to be proved -- this is the approach adopted in PHA 3. After that, I will briefly overview how a result of [Kumar et al, 2026](), which allows one to sample from a Gibbs distribution with quadratic potential and sufficiently strong external field on a wedge of $$\lbrace -1,1\rbrace^n $$ using an annealing-type procedure, is used in PHA 4 to directly sample from the localized distribution after invoking the fact that SL run for large time *will* cause the distribution to be heavilty concentrated inside a wedge of sufficient size.
-  In the final blog post, I will talk about *how* to get access to a "safe set" of points $$S_{A}(c) \subset \lbrace -1,1\rbrace^n $$ when one fixes a "good" event over the input $$A $$ so that, with high probability over the SL paths, the TAP Hessian evaluated at any $$m \in S_{A}(c) $$ is $$c $$-strongly convex. This is the place where, perhaps unsurpringly, we rely on a result about AMP state evolution. Essentially, we combine the fact that posterior means of the SL tilted measures are "stable" in time for most SL paths with the $$c $$-strong convexity of Hessians in small balls around AMP iterates, to build "tubes" around good SL paths that we then union to be contained in $$S_{A}(c) $$. From this, one straightforwardly obtains that, with high probability over $$ A $$ and the SL paths, uniformly in time $$t \in [0,T] $$, small enough balls around the SL-tilted means $$\lbrace m_t\rbrace_{0 \le t \le T} $$ lie in $$S_{A}(c) $$. 
<br>

#### Table of Contents
1. [Algorithm design](#algorithm-design)
   * [Stochastic localization and Hessian dynamics](#the-parisi-formula-and-auffinger-chen-representation)
   * [Algorithmic surrogates via the TAP free energy](#algorithmic-surrogates-via-the-tap-free-energy)
   * [ASL-TAP and boosting to a KL divergence bound](#asl-tap-and-boosting-to-a-kl-divergence-bound)
2. [Jarzynski equality](#jarzynski-equality)
   * [Rejection sampling over path space](#rejection-sampling-over-path-space)
   * [Final desiderata](#final-desiderata)
3. [Conclusions]()
   * [Certifying the overlap distribution](certifying-the-overlap-distribution)
   * [Weak functional inequalities and resolvent driven SDEs](weak-functional-inequalities-and-resolvent-driven-SDEs)
   * [Sampling without functional inequalities]()
4. [Footnotes](#footnotes)
<br>

## Algorithm Design
Let us first reason about how to bound the error between two continuous diffusion processes run for large constant time. Think of one process as representing an "ideal" process that, if run for infinite time, will localize on a point in the support of the Gibbs measure, and think of the other as one that attempts to proxy it with a covariance that can be algorithmically computed (but incurs approximation errors). 

Denote the ideal process as $$dm_t = Q(m_t)dB_t $$ and the algorithmic process as $$d\hat{m}_t = \hat{Q}(\hat{m}_t)dB_t $$. Then, by a simple argument using the fact that the Wasserstein-$$2 $$ ($$W_2 $$) distance between two Gaussians is upper bounded by the Frobenius norm of the difference between their covariances, we have

$$
W_2(\mathsf{dist}(m_t),\mathsf{dist}(\hat{m}_t)) \le \mathbb{E}\|\hat{Q}(\hat{m}_t) - Q(m_t)\|^2_F\, .
$$

This already hints at why one of the most crucial desiderata in the entire program is going to be an error estimate tracking how well a "surrogate" covariance $$\hat{Q}(\cdot) $$ tracks the "ideal" covariance $$Q(\cdot) $$ across SL paths $$\lbrace m_t \rbrace_t $$. 

It turns out, we will need to compute the error between two sets of Ito processes. The ideal processes will be a set of two coupled processes, and the algorithmic case will be similar (though not necessarily coupled). These processes are SDEs that also have drift terms. In addition to the diffusion process before, the ideal process will also include $$ dy_t = m_tdt + dB_t $$ and the algorithmic estimate will likewise include $$ d\hat{y}_t = \hat{f}_t(\hat{m}_t)dt + dB_t$$. Now, we set $$ y_0 = m_0 = \hat{y}_0 = \hat{m}_0 = 0^n $$ (which is the relevant initialization for us) and, assuming that the solutions exist and are well-posed, obtain the following representation for the ideal and algorithmic processes

$$
d\begin{pmatrix} y_t \\ m_t \end{pmatrix} = \begin{pmatrix} m_t \\ 0^n \end{pmatrix} dt + \begin{pmatrix} I_n \\ Q(m_t) \end{pmatrix}dB_t\,,
$$

and

$$
d\begin{pmatrix} \hat{y}_t \\ \hat{m}_t \end{pmatrix} = \begin{pmatrix} \hat{f}_t(\hat{m}_t) \\ 0^n \end{pmatrix} dt + \begin{pmatrix} I_n \\ \hat{Q}(\hat{m}_t) \end{pmatrix}dB_t\,.
$$

We can now upper-bound the $$W_2 $$ distance between the ideal and algorithmic processes using the cumulative difference process $$\text{err}_t = \begin{pmatrix} y_t - \hat{y}_t \\ m_t - \hat{f}_t(\hat{m}_t) \end{pmatrix} $$. Indeed, a simple application of Ito's lemma tells us that the rate of the average instantaneous error between the two processes is

$$
\frac{d}{dt}\mathbb{E}\|\text{err}_t\|^2_2 = 2\mathbb{E}\langle\text{err}_t, (\hat{f}_t(\hat{m}_t),0^n)-(m_t,0^n)\rangle + \mathbb{E}\|Q(m_t) - \hat{Q}(\hat{m}_t)\|^2_F\,,
$$

whereupon an application of a $$c $$-weighted AM-GM inequality on the first term followed by some triangle inequalities and the fact that $$(a+b)^2 \le 2a^2+2b^2 $$ tells us that 

$$
\frac{d}{dt}\mathbb{E}\|\text{err}_t\|^2_2 \le c\mathbb{E}\|\text{err}_t\|^2_2 + \frac{2}{c}\left(\underbrace{\mathbb{E}\|\hat{f}_t(\hat{m}_t)-\hat{f}_t(m_t)\|^2_2}_{\hat{f}_t \text{ Lipschitz error}} + \underbrace{\mathbb{E}\|\hat{f}_t(m_t)-m_t\|^2_2}_{\text{PHD-TAP drift error}}\right) + 2\underbrace{\mathbb{E}\|\hat{Q}(\hat{m}_t)-\hat{Q}(m_t)\|^2_F}_{\hat{Q}(\cdot)\text{ Lipschitz error}} + 2\underbrace{\mathbb{E}\|\hat{Q}(m_t)-Q(m_t)\|^2_F}_{\text{covariance estimtate error}}\,.
$$

At this point, it is clear that if we are only to run our algorithmic processes for finite time $$ T $$ and the final four terms in the bound above are of order $$O_{t,\beta}(1) $$ at every $$0 \le t \le T $$, a Gronwall's inequality bound will immediately give that $$\mathbb{E}\|\text{err}_T\|^2_2 = O_{T,\beta}(1) $$. 

Now, if we stop our algorithmic processes at some large constant time $$ T $$, we *must* either assert that some deterministic and efficient function $$h(\hat{m}_T,\hat{y}_T) $$ outputs a sample with $$o_n(1) $$ TVD error, **or** we can use the pair $$(\hat{m}_T,\hat{y}_T) $$ as a "warm start" to another efficient algorithm which outputs a sample $$\sigma \in \lbrace -1,1\rbrace^n $$ that is $$o_n(1) $$ close in TVD error to the target measure. We will choose the latter approach.

Between the bounds on the four final four quantities for the cumulative error process, the required sampler with the warm start mentioned above, and the fact that the algorithmic surrogate $$\hat{Q}(\cdot) $$ *must* be a valid covariance, we already have a list of the desiderata that the algorithm requires, at least to have $$W_2 $$ error that is $$O_{T,\beta}(1) $$[^1]:
- $$\hat{Q}(\cdot) $$ is a regular and valid covariance, namely $$c(\beta)I_n \preceq \hat{Q}(m) \preceq C(\beta)I_n $$ at all $$ m \in \lbrace -1,1\rbrace^n $$[^2].
- The covariance errors are small, that is $$\mathbb{E}\|\hat{Q}(m_t) - Q(m_t)\|^2_F \le O_{t,\beta}(1) $$.
- The functions $$\hat{f}_t(\cdot) $$ are $$C_{\beta,t}$$-Lipschitz with respect to $$\|\cdot\|_2 $$ inside the solid cube.
- The error PHD-TAP drift error for magnetization is bounded, that is $$\mathbb{E}\|\hat{f}_t(m)-m\|^2_2 \le O_{t,\beta}(1) $$ for every $$m \in \lbrace -1,1\rbrace^n $$[^3].
- After running the (discretized) algorithmic processes for $$(\hat{m}_t,\hat{y}_t) $$ for time $$T $$, there is a sampler that samples from a certain simpler ("stochastically localized") distribution with $$o_n(1) $$ TVD error in polynomial time.

While it is not clear right now, the covariance and PHD-TAP drift error estimates will only permit error of the desired order when $$ \beta < 1 $$ -- this will be discussed in the second and third blog post. Additionally, the PSD-ness property of the covariance will hold pointwise inside the cube *only* for $$ \beta < 1/2 $$. To make the PSD-ness property, and the PHD-TAP drift error estimate (which relies on it in an indirect way), work for $$\beta < 1 $$, we will need to refine these desiderata to hold *only* at a certain "safe" set $$S_A(c) \subset \lbrace -1, 1\rbrace^n $$ -- the PHD-TAP error part of this will be the subject of the third blog post, and obtaining PSD-ness over the "safe" set will be established in the final blog post. Lastly, the Lipschitz error for $$\hat{Q}(\cdot) $$ will be a consequence of the Loewner regularity of $$ \hat{Q}(\cdot) $$ combined with the resolvent structure it has based on the *explicit* choice we use -- see [(1.2)](algorithmic-surrogates-via-the-tap-free-energy).

We are now two steps away from getting a $$o_n(1) $$-TVD sampler, if we can show the desiderata outlined above. First, we need to reason about *one* more source of $$L_2 $$-error to go from $$O_{T,\beta}(1) $$-$$W_2 $$ error to $$O_{T,\beta}(1) $$-KL divergence error via an application of Girsanov's theorem. At that point, Pinsker's inequality yields a $$O_{T,\beta}(1) $$-TVD error, but that is still not $$o_n(1) $$. The second step, which achieves this, is to use rejection sampling (but over path space) to suppress the error further to $$o_{n}(1) $$-TVD error -- this is where Jarzynski's equality enters the picture. We will now introduce the ideal processes (SL/HD), then define the TAP free energy and use it to derive the algorithmic surrogates. At that point, we will be able to complete the ASL-TAP and PHD $$L_2 $$-error bound needed to apply Girsanov's theorem. We will then move on to defining Jarzynski's equality, see how it allows us to do rejection sampling, and briefly overview how it suppresses the TVD error further. Doing the last step will incur *one* more desiderata, at which point we will conclude with our final list of desiderata.
<br>

### Stochastic localization and Hessian dynamics
Given a measurable space $$(\Omega, \mathcal{B}(\Omega)) $$ a localization process is a stochastic process $$\lbrace\nu_t(\cdot)\rbrace_{t\ge 0} $$ over the space of probability measures on $$(\Omega, \mathcal{B}(\Omega)) $$. The process has two distinct properties: 1) for any $$A \in \mathcal{B}(\Omega) $$, $$\lim_{t\to\infty} \nu_t(A) \in \lbrace0,1\rbrace$$, and 2) $$\nu_t(\cdot) $$ is a martingale.

The first property tells us that the stochastic process "localizes" at one particular event $$A $$ in the measure space, and the fact that it is a martingale implies that the process is "stationary" upon averaging. The latter property allows us to think of the localization process as a convex decomposition of the measure $$\nu_0(\cdot) $$, weighted by the "tilted" measures along the way. For us, there are two important consequences of this:
1. If we can algorithmically simulate the process by "following" the tilts that generate the sequence of tilted measures for a very large (but constant) time $$T $$, then we will come close to a sample from the target measure, and 
2. We can choose *any* measure-valued process that localizes and is a martingale based on its convenience for simulating efficiently, and analzying its paths to prove our desiderata.

There is just one remaining subtlety -- while we can choose a localization process that we can analyze and simulate, where do we start it? As we will see, we choose the "linear-tilt" localization scheme which essentially starts at a sample $$x_0 $$ drawn from the target measure, and consists of noisy observations through time $$t $$ that ultimately become more informative and localize as a Dirac measure $$\delta_{x_0} $$[^4]. This scheme starts at $$y_0 = 0^n $$ and reveals itself at time $$t $$ as

$$
y_t = tx_0  + B_t\, ,
$$

where $$x_0 \sim \nu_0 $$ and $$B_t \sim \mathcal{N}(0,t) $$. The measure $$\nu_t(\cdot) $$ induced by this scheme is

$$
\nu_t(\sigma) \propto e^{\langle y_t,\sigma \rangle}\nu_0(\sigma)\,,
$$

and it is easy to see that $$\lim_{t \to\infty} \nu_t(\cdot) \to \delta_{x_0} $$ almost-surely and that $$\lbrace\nu_t(\cdot)\rbrace_{t\ge 0} $$ is a martingale process since they form a Doob martingale. At this point, using the many known equivalent characterizations of stochastic localization (see [[Sections 1 & 2, STZ26]](https://arxiv.org/pdf/2510.04460)) one can also rewrite the linear-tilt process as

$$
dy_t = m_t dt + dB_t\,,
$$

where $$m_t = \mathbb{E}_{x\sim\nu_t}\left[x\right] $$. A main conceptual insight in our work is that there is **yet** another equivalent rewrite for a stochastic localization process, and this rewrite tracks the evolution of the magnetizations/averages $$\lbrace m_t\rbrace_{t\ge 0} $$ of the "tilted" measures $$\nu_t(\cdot) $$ under the linear-tilt localization scheme $$\lbrace y_t\rbrace_{t\ge 0} $$. This process is

$$
dm_t = \mathsf{Cov}(\nu_t)dB_t\,\qquad \text{Hessian dynamics (HD)}\,, 
$$

where $$\mathsf{Cov}(\nu_t) $$ is the covariance matrix for the measure $$\nu_t $$. It is not difficult to see that this is the SDE that the magnetization process should obey from its definition as the mean of the tilted measures plus an application of Ito's lemma in conjunction with the fact that the magnetizations must form a martingale. It will turn out to be the case (for convenience more than anything else) that, though the actual algorithm will not run an algorithmic version of *this* particular version of SL, it will still be very useful to analyze it to prove some of the desiderata[^5]. The algorithmic version of this process is

$$
d\hat{m}_t s= \hat{Q}(\hat{m}_t)dB_t\,,
$$

and now a critical part of proving our desiderata and writing down the final algorithm relies on gaining access to an efficiently computable sequence of matrix-valued functions $$\hat{Q} : [-1,1]^n \to M_n(\mathbb{R})_{\text{sa}} $$ -- this is exactly where the contiguity to a planted model *and* the TAP free energy will be of assistance. 
<br>

### Algorithmic surrogates via the TAP free energy
Before we introduce the TAP free energy for the SK model, as is relevant in the high-temperature regime, let us quickly remind ourslves of the (random) Gibbs measure associated with it. The target/Gibbs measure for the SK model is given as

$$
	\mu_{\beta A}(\sigma) := \frac{e^{\frac{1}{2}\left\langle\sigma,\beta A\sigma\right\rangle}}{Z}\,,
$$

where $$ Z $$ is the normalizing constant called the partition function, and the density is defined for every $$\sigma \in \Sigma_n := \lbrace -1,1 \rbrace^n $$.

The SK model has a particularly nice form for its free energy in the high-temperature regime ($$\beta < 1 $$) whose validity was rigorously established in a series of papers -- see, for instance [[CP19]](https://arxiv.org/abs/1709.03468). The free energy is equivalent to the largest value of the TAP functional evaluated as a supremum over all possible magnetizations. The TAP functional for the SK model at some magnetization $$m \in [-1,1]^n $$ is given as

$$
	\mathcal{F}_{\mathsf{TAP}}(m) = -\frac{\beta}{2}\langle \sigma, A \sigma \rangle - \sum_{i=1}^n h(m_i) - n\left(\frac{\beta^2}{4}\left(1-\frac{\|m\|^2_2}{n}\right)\right)\,,
$$

where $$ h(\cdot) : [-1,1] \to [0,1]$$ is the entropy functional for a Bernoulli random variable (see [PHA 3, (2.2)](https://arxiv.org/pdf/2605.03718)). This structure of the TAP free energy is closely related to taking a Legendre transform of the standard definition of the free energy, and consists of adding the final term as a "correction term" to the Gibbs variational principle over product measures. For our purposes, we don't just need the TAP free energy of the standard SK model, but we need some analogue of it for any sequence of tilts $$ y_t $$ generated by the linear-tilt localization scheme we have -- this is easily accomplished by adding an external field term. Namely, setting $$ \nu_0 = \mu_{\beta A} $$ and using the fact that $$ \nu_t(m) \propto e^{\langle y_t, m \rangle}\mu_{\beta A} $$ gives us the following modification to the TAP free energy, given a fixed tilt $$ y \in \mathbb{R}^n $$

$$
	\mathcal{F}_{\mathsf{TAP}}(y,m) := -\frac{\beta}{2}\langle \sigma, A \sigma \rangle - \langle y, m \rangle - \sum_{i=1}^n h(m_i) - n\left(\frac{\beta^2}{4}\left(1-\frac{\|m\|^2_2}{n}\right)\right)\,. 
$$ 

We will be interested in viewing this as a Fenchel legendre transform, and given some algorithmic estimate $$\hat{y}_t $$ of $$ y_t $$ (with sufficiently good approximation) we will be able to compute $$ \hat{m}_t $$ using the strong convexity of this TAP free energy (since that will guarnatee the uniqueness of $$ \hat{m}_t$$ ), which will in turn allow us to compute the *next* tilt $$\hat{y}_{t+\delta t} $$. 

This is fantastic for the algorithmic procedure, and it also suggests a natural covariance $$ \hat{Q}(m) $$ that can be defined everywhere[^6] in $$[-1,1]^n $$ -- the inverse of the Hessian of the TAP free energy. This is so since $$ y_t $$ and $$ m_t $$ are dual to each other, and the Crouzeix identity tells us that the second-derivative with respect to one variable ($$y_t $$) is the functional inverse of the other ($$m_t $$) and this gives

$$
	\hat{Q}(m) := \left(\nabla^2 \mathcal{F}_{\mathsf{TAP}}(y,m)\right)^{-1} = \left(\beta^2\mathsf{tr}_n[D^{-1}(m)]I_n -\beta A + D(m)-\frac{2\beta^2}{n}mm^T\right)^{-1}\,.
$$   
 
A pleasant consequence of the choice of $$\hat{Q}(m) $$ is that it is actually independent of the tilt $$y $$. Another immediate consequence is that, under the Loewner order sandwich on $$\hat{Q}(\cdot) $$ (and consequently $$D(\cdot)\hat{Q}(\cdot) $$) required by the first desideratum for the covariance matrix, one can easily obtain that

$$
\begin{aligned}
\|\hat{Q}(a) -\hat{Q}(b)\|^2_F &= \|\hat{Q}(a)D(a)\left(D^{-1}(a)-D^{-1}(b)\right)D(b)\hat{Q}(b) -\hat{Q}(a)(R(a)-R(b))\hat{Q}(b)\|^2_F  \\
&\le \|D(a)\hat{Q}(a)\|_\infty^2\|D^{-1}(a) - D^{-1}(b)\|_F^2 + \|\hat{Q}(a)\|^2_\infty\|R(a) - R(b)\|^2_F \\
&\le C(\beta)\|a-b\|^2_2 + C'(\beta)\|R(a)-R(b)\|^2_F\, ,
\end{aligned}
$$

where $$R(a) = \left(\beta^2\mathsf{tr}_n[D^{-1}(a)]I_n-\beta A - \frac{2\beta^2}{n}aa^T\right)^{-1} $$. Then, conditioning on the event that $$\|A\|_\infty \le 2+\delta_\beta $$ and doing some elementary estimates using the fact that $$a,b \in (-1,1)^n $$, yields that $$\|R(a)-R(b)\|^2_F \le C(\beta)\|a-b\|^2_2 $$. Substituting this into the bound above shows that $$\hat{Q}(\cdot) $$ is Lipschitz with respect to $$\|\cdot\|_F $$ and proves one of the four regularity properties in the desiderata (simply as a consequence of the first desideraturm and surrogate choice of covariance).
<br>

### ASL-TAP and boosting to a KL divergence bound
The $$O(1) $$ error estimate for the $$W_2 $$ distance between the ideal and algorithmic SDEs for the magnetization and tilt we developed at the beginning will, unfortunately, not be sufficient to conlcude a $$O(1) $$ KL-divergence error via an application of Girsanov's theorem. This is because in the SL process, the magnetization $$m_t $$ is coupled with the tilt $$y_t $$, whereas this is simply not true for the algorithmic process for $$ \hat{m}_t $$ and $$ \hat{y}_t $$ when $$ \hat{Q}(\hat{m}_t) = \left(\nabla^2\mathcal{F}_{\mathsf{TAP}}(\hat{m}_t)\right)^{-1} $$. Consequently, we must account for *one* additional source of error, and that is the expected squared error between the ASL-TAP and PHD process. The final desiderata to obtain $$O(1) $$ KL divergence error will come from this part of the argument.
<br>

## Jarzynski equality
Stuff and things
<br>

### Rejection sampling over path space
Stuff and things
<br>

### Final desiderata
Stuff and things
<br>

## Conclusions
Stuff and things
<br>

### Certifying the overlap distribution
Stuff and things
<br>

### Weak functional inequalities and resolvent driven SDEs
Stuff and things
<br>

### Sampling without functional inequalities
Stuff and things
<br>

#### FOOTNOTES

[^1]: Note that even having a sampler with $$W_2 $$-error of order $$O_{T,\beta}(1) $$ is an improvement over the prior state of the art result, if it works for the entire replica-symmetric regime of $$0 <\beta< 1 $$. As we shall see, one only needs to argue that the PHD process stays close to the ASL-TAP process, which is derived in [(1.3)](#asl-tap-and-boosting-to-a-kl-divergence-bound), to obtain $$O_{T,\beta}(1) $$ error. Jarzynski's equality only enters as a final step, to boost the 

[^2]: Note that the fact that $$\hat{Q}(\cdot) $$ must be $$C$$-Lispchitz is not immediately implied by the upper bound on the Loewner order, but requires using the resolvent identity and definition of the exact choice of $$\hat{Q}(\cdot) $$ used in the algorithmic process.

[^3]: It will shortly become clear why I am calling this a "PHD-TAP drift error".

[^4]: This seems cyclical, since we are starting at a sample drawn from the measure we wish to eventually sample from. However, it will turn out that the dependency between the "tilted" measures along the localization process at every time $$t $$ and the input randomness can be decoupled -- this is because, at high-temperature ($$\beta < 1 $$), our measure turns out to be contiguous with respect to drawing the initial sample $$x_0 \sim \mathsf{Unif}\left(\lbrace-1,1\rbrace^n\right) $$ and then running the process, albeit with a "plant" term that independently gets added to the input randomness. This obviously changes the structure of our surrogate covariance $$\hat{Q}(\cdot) $$, but the added technical burden can be dealt with and is substantially easier that dealing with a localization process where the tilted measures cannot be decoupled (made independent) from the randomness of the instance. See [[Section 2, EAMS'22]](https://arxiv.org/abs/2203.05093) for more details about the planted model. 

[^5]: For the particular linear-tilt scheme that we are using as our localization process, we will be able to assert that the algorithmic process we use allows us to efficiently estimate $$\hat{m}_t $$ from $$\hat{y}_t $$. 

[^6]: Once again, this is only meaningfully true when the TAP free energy is strongly convex *everywhere* inside $$[-1,1]^n $$. This ceases to be true when $$1/\sqrt{2} < \beta < 1 $$ and we will need to "regularize" our choice of $$\hat{Q} $$ for that situation -- this is one of the key conceptual upgrades in PHA 4. We will discuss this in the final blog post. 
---
layout: post
title:  "Representations of Finite Groups and Character Theory"
date:   2023-11-14 11:15:40
blurb: "Irreducible representations, Schur's orthogonality relations & Group Algebra"
og_image: /assets/img/content/post-example/Banner.jpg
---

<br />
We will ...
<br />

#### Table of Contents
1. [Representations and Subrepresentations](#representations-and-subrepresentations)
   * [Types of representations and equivalences](#types-of-representations-and-equivalences)
   * [Irreducibility and Schur's lemma](#irreducibility-and-schur's-lemma)
2. [Character Theory](#character-theory)
    * [Characters of representations](#characters-or-representations)
    * [Schur's orthogonality relations for characters](#Schur's-orthogonality-relations-for-characters)
    * [Decompositions of representations](#decompositions-of-representations)
    * [Induced representations](#induced-representations)
3. [Group Algebras](#ruelle-probability-cascades)
    * [Decomposition of the group algebra](#decomposition-of-the-group-algebra)
    * [Center of the group algebra](#center-of-the-group-algebra)
    * [Integrality properties of characters](#integrality-properties-of-characters)
4.  [Mackey's Irreducibility Criterion](#mackey's-irreducibility-criterion)
    * [Characters of induced representations](#characters-of-induced-representations)
    * [Frobenius reciprocity](#frobenius-reciprocity)
    * [Degrees of irreducible representations]()
5. [Footnotes](#footnotes)

<br />

## Representations and Subrepresentations
The primary concern of this section is to introduce the notion of representations (for finite groups) formally and then demonstrate that these representations are unique (up to certain isomorphisms), illustrating the three primary types of representations. We will then introduce irreducible representations as those that have no proper ''stable'' subrepresentations, and show that subspaces of representations that are ''stable'' under the group action have similar ones in their orthogonal complements. We will conclude by proving Schur's lemma, which allows one to argue that there exists a unique (up to isomorphism) decomposition of a representation into finitely many irreducible ones.

### Types of representations and equivalences
To briefly enumerate certain common representations, we first define a representation.

<u><strong>(Representation of G)</strong></u>: Given a finite group $$G $$ and a vector space $$V $$ over $$\mathbb{C} $$, a representation $$\rho : G \to \mathsf{GL}(V) $$ is a group homomorphism.

The notation $$\mathsf{GL}(V) $$ simply represents the set of _invertible_ linear operators $$T : V \to V $$. For finite dimensional vector spaces $$V $$, this is equivalent to the set of invertible square matrices. The main idea is to encode structural information about the group $$G $$ into decompositions of linear operators, which can be studied transparently using linear algebra.

Representations need not be unique, and therefore, it is important to glance at the various types of representations one can have and which are considered canonical. To do this, we first define some notion of ''equivalence'' between representations of a group.

<u><strong>(Equivalent Representations)</strong></u>: Given representations $$\rho_1: G \to \mathsf{GL}(V) $$, $$\rho_2: G \to \mathsf{GL}(W) $$ and a linear map $$\phi: V \to W $$,
$$ \begin{equation}\rho_1 \sim \rho_2 \iff \forall g \in G,\, \rho_2(g)\cdot \phi = \phi \cdot \rho_1(g)\, . \end{equation} $$

The above definition simply asserts that two representations are equivalent when there exists an invertible linear map $$\phi $$ between $$V $$ and $$W $$ that is evidence for the fact that $$\rho_1(g) $$ and $$\rho_2(g) $$ are similar matrices, for every $$g \in G$$. Equivalently,
$$ \begin{equation} \rho_2(g) = \phi \cdot \rho_1(g) \cdot \phi^{-1},\, \forall g \in G. \end{equation} $$
This makes it clear that two representations are equivalent if they preserve the group structure embedded up to conjugation by some invertible linear transformations on their various linear maps. As similarity between finite dimensional linear maps is an equivalence relation, so is equivalence between representations.

Three common types of representations are:
1. *Characters:*
2. *Regular representation:*
3. *Permutation representation:*

### Irreducibility and Schur's Lemma

# PHYSICS & COMPUTATION READING GROUP

This will be a reading group across Tufts, Harvard and MIT. The focus will be on analytic and probabilistic tools used to study the mathematical properties and algorithmic tractability of approximating (and sampling from statistics of) low-energy states of certain families of classical hamiltonians (spin glasses) and quantum hamiltonians (bosonic/fermionic systems).

**Meeting Time**: Thursdays, 11 AM - 1 PM EST

**Location**: [Zoom room](https://harvard.zoom.us/j/2185567693?pwd=TzBGeHZIRWJRWXp3OVdOdzF6Qk1uZz09)

**Organizers**: [Saeed Mehraban](https://sites.google.com/view/saeedmehraban/about), [Juspreet Singh Sandhu](https://juspreetsandhu.me)

# Talk Schedule

- Lecture-1: Motivation & Introduction - Approximability and hardness in physics and computation 
  - Speaker(s): Saeed, Juspreet
  - Date: 02/02/2023

- Lecture-2: Barvinok's method
  - Speaker: Saeed
  - Date: 02/09/2023

- Lecture-3: Deterministic counting and approximation of the permanent and partition functions I 
  - Speaker: Saeed
  - Date: 02/16/2023

- Lecture-4: Deterministic counting and approximation of the permanent and partition functions II
  - Speaker: Saeed
  - Date: 02/23/2023

- Lecture-5: Overlap concentration in spin glasses and sparse random Max-CSPs I 
  - Speaker: Juspreet
  - Date: 03/02/2023

- Lecture-6: Overlap concentration in spin glasses and sparse random Max-CSPs II
  - Speaker: Juspreet
  - Date: 03/09/2023

- Lecture-7: Matrix sum-of-squares optimization for spherical spin glasses
  - Speaker: Juspreet
  - Date: 03/16/2023

- Lecture-8: Stability of low-degree polynomials & Langevin dynamics
  - Speaker: Juspreet
  - Date: 03/23/2023

- Lecture-9: _Open Discussion Session_: Bridging the gap between different approximation techniques
  - Moderators: Saeed, Juspreet
  - Date: 03/30/2023

- Lecture-10: TBA
  - Speaker: TBA
  - Date: TBA


# Relevant Papers
For a survey of Barvinok's method, a good reference is the monograph by Barvinok [[Bar]](https://link.springer.com/book/10.1007/978-3-319-51829-9) with Chapter-2 giving an overview of the technical toolkit, which uses complex analysis to study the roots of real-valued polynomials. Another source are the talks on "Computing Partition Functions" from the Simon's workshop on the [geometry of polynomials](https://www.youtube.com/watch?v=TUjCLXPqW2Y&list=PLgKuh-lKre13XzHXH_rnq0ptd3ahU5TfB). A good historical reference for the overlap-gap property and its use in obstructing algorithms on random optimization problems is the survey by Gamarnik [[G21]](https://arxiv.org/pdf/2109.14409.pdf). A more technical survey is the one by Auffinger, Montanari and Subag [[AMS22]](https://arxiv.org/pdf/2206.10217.pdf). The three most accessible and comprehensive texts on the proof of the Parisi formula along with the applications of these techniques to solve other problems in mathematical and statistical physics are the ones by Talagrand ([[Vol. 1]](https://link.springer.com/book/10.1007/978-3-642-15202-3),[[Vol. 2]](https://link.springer.com/book/10.1007/978-3-642-22253-5)) and Panchenko [[Pan13]](https://link.springer.com/book/10.1007/978-1-4614-6289-7).

## Deterministic Counting and Approximation of Quantum Hamiltonians
- Statistical Theory of Equations of State and Phase Transitions. I. Theory of Condensation [[YL52]](https://journals.aps.org/pr/abstract/10.1103/PhysRev.87.404).
- Completely Analytical Interactions: Constructive Description [[DS87]](https://link.springer.com/article/10.1007/BF01011153).
- Polynomial-Time Approximation Algorithms for the Ising Model [[JS93]](https://www.math.cmu.edu/~af1p/Teaching/MCC17/Papers/JSIsing.pdf).
- A Polynomial-Time Approximation Algorithm for the Permanent of a Matrix with Nonnegative Entries [[JSV]](https://faculty.cc.gatech.edu/~vigoda/Permanent.pdf).
- The Ising Partition Function: Zeros and Deterministic Approximation [[LSS19]](https://link.springer.com/article/10.1007/s10955-018-2199-2).
- Deterministic Polynomial-Time Approximation Algorithms for Partition Functions and Graph Polynomials [[PR17]](https://arxiv.org/pdf/1607.01167.pdf).
- Classical Algorithms, Correlation Decay, and Complex Zeros of Partition Functions of Quantum Many-Body Systems [[HMS19]](https://arxiv.org/pdf/1910.09071.pdf). 
- Approximating the Determinant of Well-Conditioned Matrices by Shallow Circuits [[AEM19]](https://arxiv.org/pdf/1912.03824.pdf).
- Approximating the Permanent of a Random Matrix with Vanishing Mean [[EM18]](https://arxiv.org/pdf/1711.09457.pdf).
- The Computational Complexity of Linear Optics [[AA10]](https://arxiv.org/pdf/1011.3245.pdf).
- The Computational Hardness of Counting in Two-Spin Models on d-Regular Graphs [[SS12]](https://arxiv.org/pdf/1203.2602.pdf).

## Analysis and Probability in Mean-Field Spin-Glasses
### Interpolations, Ultrametricity and the TAP Approach
- The Thermodynamic Limit in Mean Field Spin Glass Models [[GT02]](https://arxiv.org/pdf/cond-mat/0204280.pdf).
- Broken Replica Symmetry Bounds in the Mean Field Spin Glass Model [[Gue02]](https://arxiv.org/pdf/cond-mat/0205123.pdf).
- The Parisi Formula [[Tal06]](https://annals.math.princeton.edu/wp-content/uploads/annals-v163-n1-p04.pdf).
  - On the Proof of the Parisi Formula by Guerra and Talagrand [[Bol05]](http://www.numdam.org/item/SB_2004-2005__47__349_0.pdf).
- Ghirlanda-Guerra Identities and Ultrametricity: An Elementary Proof in the Discrete Case [[Pan11]](https://arxiv.org/pdf/1106.3984.pdf).
- The Parisi Ultrametricity Conjecture [[Pan13]](https://arxiv.org/pdf/1112.1003.pdf).
  - The Parisi Formula for Mixed p-Spin Models [[Pan14]](https://arxiv.org/pdf/1112.4409.pdf).
- On the Energy Landscape of the Mixed Even p-Spin Model [[CHL18]](https://arxiv.org/pdf/1609.04368.pdf).
- Thouless-Anderson-Palmer Equations for Generic p-Spin Glasses [[AJ18]](https://arxiv.org/pdf/1612.06359.pdf).
  - The Thouless-Anderson-Palmer Equation in Spin Glass Theory [[Bol08]](https://anr-malin.sciencesconf.org/data/pages/Aussois_2.pdf).
- The Generalized TAP Free Energy I, II [[CPS18]](https://arxiv.org/pdf/1812.05066.pdf), [[CPS21]](https://arxiv.org/pdf/1903.01030.pdf).

### Analytic Properties of the Parisi Formula
- On Differentiability of the Parisi Formula [[Pan08]](https://arxiv.org/pdf/0709.1514.pdf).
- On Properties of Parisi Measures [[AC13]](https://arxiv.org/pdf/1303.3573.pdf).
- The Parisi Formula has a Unique Minimizer [[AC14]](https://arxiv.org/pdf/1402.5132.pdf).
- Variational Representations for the Parisi Functional and the Two-Dimensional Guerra-Talagrand Bound [[C16]](https://arxiv.org/pdf/1501.06635.pdf).
- The SK model is Infinite Step Replica Symmetry Breaking at Zero Temperature [[ACZ20]](https://arxiv.org/pdf/1703.06872.pdf).

## Overlap Concentration, Low-Degree Stability and Overlap-Gap Properties: Algorithmic Hardness
- Limits of Local Algorithms over Sparse Random Graphs [[GJ14]](https://arxiv.org/pdf/1304.1831.pdf).
- Local Algorithms For Independent Sets Are Half-Optimal [[VR17]](https://arxiv.org/pdf/1402.0485.pdf).
- Suboptimality of Local Algorithms for a Class of Max-Cut Problems [[CGPR19]](https://arxiv.org/pdf/1707.05386.pdf).
- The Overlap Gap Property and Approximate Message Passing Algorithms for p-Spin Models [[GJ21]](https://projecteuclid.org/journals/annals-of-probability/volume-49/issue-1/The-overlap-gap-property-and-approximate-message-passing-algorithms-for/10.1214/20-AOP1448.short).
- Hardness of Random Optimization Problems for Boolean Circuits, Low-Degree Polynomials, and Langevin Dynamics [[GJW20]](https://arxiv.org/pdf/2004.12063.pdf).
- Tight Lipschitz Hardness for Optimizing Mean Field Spin Glasses [[HS21]](https://arxiv.org/pdf/2110.07847.pdf).
- Limitations of Local Quantum Algorithms on Random Max-k-XOR and Beyond [[CLSS22]](https://arxiv.org/pdf/2108.06049.pdf).
- Random Max-CSPs Inherit Algorithmic Hardness from Spin Glasses [[JMSS22]](https://arxiv.org/pdf/2210.03006.pdf).

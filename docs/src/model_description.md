!!! warning "Warning"
    This page is under construction.
# Introduction 

Prior research has found that the judged likelihood of a hypothesis often depends on the order in which evidence is presented. In other words, the final judgment that a hypothesis is true is different for evidence sequence $e_i, e_j$ and evidence sequence $e_j, e_i$. The goal of this tutorial is to describe a quantum order effect model (QOEM) as it applies to a medical diagnosis task. The basic model can be adapted to other tasks in which evidence is evaluated sequentially.  

# Medical Diagnosis Task

Medical doctors routinely perform medical diagnosis as part of their job requirements. Medical diagnosis is the process of assessing the probability a patient has a disease based on symptoms and medical tests. Below, we use the QOEM to understand a how a person evaluates the hypothesis that a fictitous patient has a disease (d) after evidence from two informatino sources are presented sequentially. The judgments are made in two conditions, which vary the presentation order of the evidence. In one condition, a person makes the following judgments:

1. the probability of disease based on initial symptoms 
2. update the probability of disease after the medical history $S_i$ provides positive evidence for the disease
3. update the probability of disease after the laboratory test $S_j$ provides negative evidence for the disease

The procedure for the other condition is identical except the order of the information sources $S_i,S_j$ is reversed: 

1. the probability of disease based on initial symptoms 
2. update the probability of disease after the laboratory test $S_j$ provides negative evidence for the disease
3. update the probability of disease after the medical history $S_i$ provides positive evidence for the disease

# Order Effect 

An order effect occurs when the final probability judgment of disease depends on the order in which evidence is presented. An order effect can be stated formally as: 

$\Pr(D=d \mid S_i=1, S_j=-1) \ne \Pr(D=d \mid S_j=-1, S_i=1).,$

where $1$ indicates positive evidence and $-1$ indicates negative evidence for the disease. Order effects are challenging to explain for classical probability models because evidence is commutative. One workaround for classical probability is to augment the sample space with an additional event representing evidence order. However, augmenting the sample space comes at the cost of creating additional parameters which need to be specified and justified, or estimated from limited data.

# Quantum Order Effect Model 

According to the QOEM, positive and negative evidence form *incompatible* events, meaning they cannot be considered simultaneously. As a consequence, a person cannot represent the full joint distribution of events over disease status (present vs. absent), evidence type (positive vs. negative), and information source (medical history vs. laboratory test). Instead, incompatible events are represented in a lower dimensional space by different bases, which must be evaluated sequentially. Bases correspond to different perspectives of a situation. In the medical diagnosis task, the QOEM assumes that positive and negative evidence are viewed one at a time from different perspectives. Order effects arise from this process because the linear algebra operations described below are non-commutative. 

## Basis

The first in defining a quantum model is to determine which events are compatible and which are incompatible. The QOEM assumes information source is incompatible, but disease status and evidence type are compatible. Consequentially, the basis of the QOEM consists of four states corresponding to the all possible combinations of disease status and evidence type:

1. disease present (p) and positive evidence for disease (p)
2. disease present (p) and negative evidence for disease (n)
3. disease absent (a) and positive evidence for disease (p)
4. disease absent (a) and negative evidence for disease (n)

where the values in the parentheses correspond to indices. In the notation used below, the first index corresponds to the presence or absence of the disease, and the second index corresponds to the type of evidence for the disease (positive or negative). Information source is represented as different bases within the reduced 4 dimensional space. Each basis for information source is related to the other bases through a rotation factor, which corresponds to the idea of viewing the diagnosis from different perspectives. More formally, the basis is given by the following orthonormal vectors in the standard position:  

```math
\mathbf{B} = \{\ket{\textrm{B}_{p,p}},\ket{\textrm{B}_{p,n}},\ket{\textrm{B}_{a,p}},\ket{\textrm{B}_{a,n}}\},
```
where each basis vector consists of a 1 with all other elements equal to zero, e.g., 

$\ket{\textrm{B}_{p,p}} = \begin{bmatrix}
	1 \\ 
	0 \\ 
	0 \\ 
	0\\ 
\end{bmatrix}$
 
Combining the basis vectors into a single matrix, we get the identity matrix:

```math
\mathbf{I}_4 = \begin{bmatrix}		
	1 & 0 & 0 & 0\\
	0 & 1 & 0 & 0\\
	0 & 0 & 1 & 0\\
	0 & 0 & 0 & 1\\
\end{bmatrix}
```
## States

The state of the cognitive system is a superposition (i.e. linear combination) over basis states:

$\ket{\boldsymbol{\psi}} = \alpha_{p,p} \ket{B_{p,p}} + \alpha_{p,n} \ket{\textrm{B}_{p,n}}+ \alpha_{a,p}\ket{\textrm{B}_{a,p}}+ \alpha_{a,n} \ket{\textrm{B}_{a,n}},$

where $\lVert\ket{\boldsymbol{\psi}} \rVert = 1$. The coefficients can be written as:

$\boldsymbol{\alpha} = \begin{bmatrix}
	\alpha_{p,p} \\ 
	\alpha_{p,n} \\ 
	\alpha_{a,p} \\ 
	\alpha_{a,n} \\ 
\end{bmatrix}.$

The initial state is is based on the mean probability judgment after the first symptoms are described: 

$\ket{\boldsymbol{\psi}} = \begin{bmatrix}
	\sqrt(\frac{.676}{2}) \\ 
	\sqrt(\frac{.676}{2}) \\ 
	\sqrt(\frac{.324}{2}) \\ 
	\sqrt(\frac{.324}{2}) \\ 
\end{bmatrix}.$

## Projectors 

The QOEM makes repeated use of three projection matrices. The following projector is used to evaluate the probability of disease:

$\mathbf{P}_d = \ket{\textrm{B}_{p,p}} \bra{\textrm{B}_{p,p}} + \ket{\textrm{B}_{p,n}} \bra{\textrm{B}_{p,n}} = \begin{bmatrix}		
	1 & 0 & 0 & 0\\
	0 & 1 & 0 & 0\\
	0 & 0 & 0 & 0\\
	0 & 0 & 0 & 0\\
\end{bmatrix}.$

Notice it spans the 2D sub-space in which the disease is present. The next projector is used to evaluate the probability of positive evidence:

$\mathbf{P}_p = \ket{\textrm{B}_{p,p}} \bra{\textrm{B}_{p,p}} + \ket{\textrm{B}_{a,p}} \bra{\textrm{B}_{a,p}} = \begin{bmatrix}		
	1 & 0 & 0 & 0\\
	0 & 0 & 0 & 0\\
	0 & 0 & 1 & 0\\
	0 & 0 & 0 & 0\\
\end{bmatrix}.$

Similarly, the projector spans the 2D sub-space in which positive evidence is discovered. Finally, the following projector evaluates the probability of negative evidence:

$\mathbf{P}_n = \ket{\textrm{B}_{p,n}} \bra{\textrm{B}_{p,n}} + \ket{\textrm{B}_{a,n}} \bra{\textrm{B}_{a,n}} = \begin{bmatrix}		
	0 & 0 & 0 & 0\\
	0 & 1 & 0 & 0\\
	0 & 0 & 0 & 0\\
	0 & 0 & 0 & 1\\
\end{bmatrix}.$
As before, the projector spans the 2D sub-space in which negative evidence is discovered. 

## Hamiltonian Matrices

Hamiltonian matrices govern the decision dynamics of the model. The Hamiltonian matrix $\mathbf{H}$ consists of two components: $\mathbf{H}_A$ is sensitive to the payoff matrix, and $\mathbf{H}_B$ is sensitive to cognitive dissonance between beliefs and actions. The component $\mathbf{H}_A$ is defined as follows: 

$\mathbf{H}_A = \begin{bmatrix}		
	\mathbf{H}_{A_d} & \mathbf{0}\\
	\mathbf{0} & \mathbf{H}_{A_c}\\
\end{bmatrix},$
where

$\mathbf{H}_{A_k} = \frac{1}{\sqrt{1 + \mu_k^2}}\begin{bmatrix}		
	\mu_k & 1\\
	1 & -\mu_k\\
\end{bmatrix}.$


## Evidence Evaluation

This selection describes the process of selecting an action and determining the defection probability. The time evolution is governed by the unitary transformation matrix which is given by:

$\mathbf{U}_h= e^{-i \cdot t \cdot  \mathbf{H}_p},$

$\mathbf{U}_l = e^{-i \cdot t \cdot  \mathbf{H}_n},$


## QOEM Predictions

First, we will compute the probability the disease is present before evidence from the tests is presented. The probability is the same for both orders. The computation involves projecting the intial state onto the sub-space spanning disease present and squaring the magnitude of the projection, which is computed as:

$\Pr(D=d) = \lVert \mathbf{P}_d \ket{\psi} \rVert^2.$

Next, suppose we update the assessment of the disease given positive evidence from test $i$ is found. The process involves three computations: (1) projecting the initial state onto the positive evidence sub-space, (2) normalizing the state vector, and (3) projecting the new state onto the sub-space spanning disease present. To update the state, we first project onto the sub-space representing positive evidence.

$\ket{\psi_p^\prime} = \mathbf{P}_p \mathbf{U}_h\ket{\psi}.$ 

Note that he unitary matrix $\mathbf{U}_p$ changes the system to the basis for positive evidence. After projecting on to the sub-space for positive evidence, the state collapses to the positive evidence sub-space because the evidence from the first test is \emph{known}$ to be positive. As a result, the new state must be normalized such that the squared magnitude is 1:

$\ket{\psi_p} = \frac{\ket{\psi_p^\prime}}{\lVert \ket{\psi_p^\prime} \rVert}.$ 

The last step involves projecting onto the sub-space spanning disease present and squaring the magnitude to obtain the revised probability judgment for the disease:

$\Pr(D=d \mid E_i = 1) = \lVert \mathbf{P}_d \ket{\psi_p} \rVert^2.$

Next, negative evidence is presented from the second test. A similar sequence of steps is followed. The first step involves projecting onto the sub-space spanned by negative evidence: 

$\ket{\psi_{pn}^\prime} = \mathbf{P}_n \mathbf{U}_l \mathbf{U}_p^\dagger \ket{\psi_p}.$ 

Notice that the process of changing bases involves an extra step. Because the model is currently in the positive evidence basis, it is necessary to perform the inverse operation with the conjugate transpose $\mathbf{U}_p^\dagger$ and then switch to the basis for negative evidence with $\mathbf{U}_l$. As before, the state collapses to negative evidence because negative evidence was observed. Thus, the state must be normalized:

$\ket{\psi_{pn}} = \frac{\ket{\psi_{pn}^\prime}}{\lVert \ket{\psi_{pn}^\prime} \rVert}.$ 

Now that the state is in the sub-space for negative evidence, the last step involves projecting onto the 
disease present sub-space to obtain a probability judgment for the presense of the disease:

$\Pr(D=d \mid E_i = 1, E_j=-1) = \lVert \mathbf{P}_d \ket{\psi_{pn}} \rVert^2.$



|  Condition   | Formula    | Data | Model |
| :-- | :-- | :-- | :-- |
|  1   |  $\Pr(R_1=d \mid R_2=d)$   | .84| .81|
|  2  |   $\Pr(R_1=d \mid R_2=c)$   | .66| .65|
|  3   |  $\Pr(R_1=d)$  | .55 | .57|

The code used to generate the predictions can be viewed by expanding the code block below:

### Dynamics 

The plot below shows the dynamics of the model for each condition.

### Interference Effects

The plot below shows the interference effect as a function of $\mu$ for multiple values of $\gamma$. In the simulations below, we fix $t=\frac{\pi}{2}$.

# References

Trueblood, J. S., & Busemeyer, J. R. (2011). A quantum probability account of order effects in inference. Cognitive science, 35(8), 1518-1552.
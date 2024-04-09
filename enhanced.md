---
title: Enhanced sampling
exports:
    - format: pdf
      template: plain_latex
---

```{note}
In this chapter I will use *reaction* as a general term for a microscopic process that transforms the system of interest (be it a macromolecule, an extended many-body system, a collection of atoms, *etc.*) between two well-defined states, $A$ and $B$.
```

## Reactions and rare events

Consider a system that can switch, possibly reversibly, between two macrostates, $A$ and $B$. Here the term macrostate is used loosely to indicate ensembles of microstates where the system resides for times that are much larger than the microscopic characteristic time; in thermodynamic parlance, $A$ and $B$, which are sometimes called *basins*, should be either metastable or equilibrium states, and therefore separated by a free-energy barrier $\Delta F_b$ larger than the thermal energy.

```{hint} Some examples
Examples relevant to computational biophysics are processes involving protein folding and unfolding, nucleic acid hybridisation, or switching between different conformations of the same (macro)molecule.
```

In this context the free-energy barrier[^activation_energy] between $A$ from $B$, $\Delta F_b^{A \to B} = F_{\rm max} - F_A$, is defined as the difference between the free energy of $A$, $F_A$ and that of the transition state, $F_{\rm max}$, which is the highest free-energy point along the reaction pathway connecting $A$ to $B$[^delta_F]. Note that $\Delta F_b^{A \to B}$ controls not only the probability of jumping from $A$ to $B$, but also the the rate of the reaction, which is proportional to $e^{-\beta \Delta f_b}$ (see *e.g.* [](https://doi.org/10.1063/1.1749604
) and [](https://doi.org/10.1039/TF9353100875)). See [](#fe_barrier) for a graphical definition of these quantities. 

But what is a "reaction pathway"? While a more precise answer will be given in [](#reaction-coordinates), here I just mention that it is possible to construct a *reaction coordinate* that can be used to gauge the progress (extent) of a reaction. Such a reaction coordinate, for which I will use the symbol $\xi$, is in general a function of the microscopic degrees of freedom of the system, $\dofs$, and its definition is not unique.

```{figure} figures/fe_barrier.png
:name: fe_barrier 
:align: center

An example of the free energy landscape of a system displaying two basins separated by a transition state.
```

It is often the case that what interests us is the *reaction* itself rather than the $A$ and $B$ states, which are often known (and possibly have been characterised) beforehand. In this case, simulations starting from one of the two states, say $A$, would remain in $A$ for sometime, then quickly jump to state $B$, where it would again reside for some time before switching basin once again, and so on. Moreover, if the free-energy barrier between the two basins is large ($\delta F_b \gt 10 k_BT$), the number of transitions from $A$ to $B$ and back will be very small. Therefore, using unbiased simulations[^unbiased] to sample the transition itself, for instance to evaluate the free-energy landscape as in [](#fe_barrier), requires a large computational effort which is mostly wasted in sampling uninteresting parts of the phase space. 

In this part we will understand how the sampling of the transition from $A$ to $B$ can be enhanced by using advanced computational techniques (collectively known as [rare event sampling techniques](https://en.wikipedia.org/wiki/Rare_event_sampling)).

[^activation_energy]: Sometimes also called *activation (free) energy*.
[^delta_F]: Note that, per this definition, $\Delta F_b^{A \to B} \neq \Delta F_b^{B \to A} = F_{\rm max} - F_B$.
[^unbiased]: In unbiased simulations the system evolves according to the fundamental laws of physics (Newton's laws of motion in this context), without any artificial bias or constraints imposed.

(reaction-coordinates)=
## Reaction coordinates

Before we delve into the methods themselves, we need to introduce a framework that makes it possible to describe the reaction we are interested in. We do so by defining one or more *reaction coordinates* (also known as *order parameters* in physics lingo), which are *collective variable* that are function of (all or a subset of) the microscopic degrees of freedom $\dofs$. 

In the context of molecular simulations, a reaction coordinate (RC) is a parameter or set of parameters that describe the progress of a given reaction. It serves as a way to quantify and track the changes occurring within a system as it evolves through the reaction pathway. The primary purpose of a reaction coordinate is to provide a simplified description of the system's thermodynamics, making it possible to monitor and analyze the progress of a reaction in terms of a single or a few variables: by using a reaction coordinate we are reducing the complexity of a many-body system with many degrees of freedom to obtain a simplified description that can be used to investigate the reaction itself, effectively applying a *dimensionality reduction* procedure. This simplification is essential for understanding the microscopic underpinnings of the reaction of interest. Defining a reaction coordinate makes it possible to draw a diagram such as the one shown in [](#fe_barrier), which are often called free-energy profiles or landscapes, where the variation of the free energy along a particular reaction coordinate or collective variable is plotted.

The choice of the RC depends on the specific process being studied and it is not, in general, unique. It could be a simple geometric parameter such as bond length, bond angle, or dihedral angle, or it could be a more complex collective variable that captures the overall structural changes in the system, such as the distance between two key functional groups, the position or coordination number of a particular atom, or the solvent-accessible surface area of a biomolecule.

```{note}
In general, a reaction coordinate can be multidimensional, and sometimes it is useful, or even necessary, to define such complex variables. However, for the sake of simplicity here I will use unidimensional reaction coordinates, which will be designated by the symbol $\xi = \xi(\dofs)$.
```

Once we have chosen a RC to characterise the reaction of interest, which is in general not an easy task, and an active area of research in itself, we can use it to describe the reaction. We can formally see how if we calculate the *partially-integrated partition function* (see *e.g.* [](10.33011/livecoms.4.1.1583)) by integrating the Boltzmann factor over all the degrees of freedom at constant $\xi$, *viz.*:

$$
Q(\xi) = \int_V e^{-\beta H(\dofs)} \delta(\xi - \xi(\dofs)) d\dofs,
$$ (Q_unbiased)

where $\delta(\cdot)$ is the Dirac delta distribution function. Note that the same procedure has been already carried out in the context of [coarse graining](./coarse_grained.md), where the partial integration on the phase space made it possible to obtain a simplified description of the system, showing that the process of dimensionality reduction we apply is strictly the same in the two cases.

In turn, $Q(\xi)$ can be used to obtain a free-energy profile (sometimes called free-energy surface is $\xi$ is multidimensional) such as the one presented in [](#fe_barrier), which is defined as

$$
F(\xi) = -k_B T \ln (Q(\xi)).
$$

Now consider an observable that can be written as a function of $\xi$. Its ensemble average can be formally written as

$$
\langle O(\xi(\dofs) \rangle = \frac{\int_V e^{-\beta H(\dofs)} O(\xi(\dofs)) d\dofs}{\int_V e^{-\beta H(\dofs)} d\dofs},
$$

which can be simplified by using Eq. [](#Q_unbiased) to

$$
\langle O(\xi) \rangle = \frac{\int_{\xi_{\rm min}}^{\xi_{\rm max}} O(\xi) Q(\xi) d\xi}{\int_{\xi_{\rm min}}^{\xi_{\rm max}} Q(\xi) d\xi} = \int_{\xi_{\rm min}}^{\xi_{\rm max}} O(\xi) P(\xi) d\xi,
$$

where $\xi_{\rm min}$ and $\xi_{\rm max}$ correspond to the minimum and maximum values of $\xi$, and we have defined the *marginal probability density*

$$
P(\xi) = \frac{Q(\xi)}{\int_{\xi_{\rm min}}^{\xi_{\rm max}} Q(\xi) d\xi} = \frac{Q(\xi)}{Q},
$$ (marginal_P)

where we have used the partition function $Q = \int_{\xi_{\rm min}}^{\xi_{\rm max}} Q(\xi) d\xi$.


```{hint} A simple example
Consider a system formed by two DNA strands made of $N_1$ and $N_2$ units (atoms or coarse-grained beads), respectively, that can hybridise: a possible (although not exactly ideal) reaction coordinate could be the distance between the centres of mass of the two strands, *viz.*:

$$
\xi = \left| \vec{R}^{\rm cm}_1 - \vec{R}^{\rm cm}_2 \right| = \left| \frac{1}{N_1} \sum_{i \in N_1} \vec{r}_i - \frac{1}{N_2} \sum_{i \in N_2} \vec{r}_i \right|.
$$
```

## Umbrella Sampling

The first technique I will present is the venerable *umbrella sampling* (US), which was introduced in the Seventies by [Torrie and Valleau](https://doi.org/10.1016/0021-9991(77)90121-8).

The basic idea behind umbrella sampling is to bias the system along a chosen reaction coordinate by adding a so-called *biasing potential*[^umbrella] that confines the system to different regions (also known as *windows*) along that coordinate. By running multiple simulations with biasing potentials centred on different points along the reaction coordinate, the entire range of interest can be sampled. After sufficient sampling is done in each window, the bias introduced by the additional potentials can be removed to obtain the unbiased free energy profile (or any other observable of interest) along the reaction coordinate.

A typical umbrella sampling simulation thus comprises several steps, which I will discuss separately.

### 1. Choosing the reaction coordinate

This is arguably the most important step, since choosing a sub-optimal RC can sometimes massively increase the required simulation time. Fortunately, most of the times the choice is either obvious (*e.g.* the concentration of the product in a chemical reaction), or dictated by the observable(s) of interest (see below for an [example](#a-real-world-example)).

### 2. Selecting a biasing potential

The role of the biasing potential $V^{\rm bias}(\xi)$ is to confine a system within a (usually rather narrow) region of the reaction coordinate. As such it must be a function of the reaction coordinate(s) only, without any explicit dependence on any of the microscopic $\dofs$. The most common choice is a harmonic potential, whose shape gives the method its name and usually takes the form

$$
V^{\rm bias}(\xi) = \frac{1}{2} K (\xi - \bar\xi)^2,
$$

where $\bar\xi$ is the position of the minimum of the potential and $K$ is the strength of the resulting spring force. Other choices are possible (see *e.g.* [here](https://doi.org/10.1039/C4SM02218A) or [here](https://doi.org/10.1063/1.1739216) for examples of biases that are not differentiable and therefore can only be used in Monte Carlo simulations).

### 3. Partitioning the reaction coordinate into windows

Next, we need to split the range of interest, $[\xi_{\rm min}, \xi_{\rm max}]$, into windows. The most common strategy is to divide the reaction coordinate into equispaced windows centred on $\xi_1, \xi_2, \xi_3, \ldots$, with $i \in [1, N]$, where $N$ is the total number of windows (and hence of independent simulations). The distance between two neighbouring windows, $\Delta \xi_i = \xi_{i + 1} - \xi_i$, which is often taken as a constant, should be chosen carefully: on one hand it should be as large as possible to make $N$ as small as possible; on the other hand, $\Delta \xi$ should be chosen so that there is some overlap between adjacent windows to prevent discontinuities in the free energy profile. This is to ensure that the neighboring windows provide sufficient sampling for accurate reweighting. An often good-enough first estimate can be made by assuming that $P(R_i) \approx P(R_{i + 1})$, and then by choosing a $\Delta \xi$-value for which $V^{\rm bias}(\Delta \xi_2)$ is of the order ot $k_B T$, *i.e.* that the value of the biasing potential calculated in the midpoint separating two neighbouring windows is of the order of the thermal energy.

In practice, the number, size and spacing of the windows depends on the curvature of the free energy profile along the reaction coordinate, which is not known beforehand. Smaller windows may be needed in regions with steep gradients or large energy barriers, while larger windows may suffice in more gradually changing regions. Fortunately, given the independent nature of the simulations that run in each window, the partitioning can be improved upon *a posteriori*: if one realises that the explored range is not sufficient, it can be extended by adding simulations with biasing potentials centred beyond $\xi_{\rm min}$ and/or $\xi_{\rm max}$. Sampling can also be improved by adding simulations in regions of the RC where the $P(R)$ is steeper.

At the end of this procedure, each window will be assigned a biasing potential $V^{\rm bias}_i(\xi)$.

```{hint} Adaptive sampling
There are more advanced methods, where the size and placement of windows are adjusted dynamically based on the evolving free energy landscape observed during the simulation (see *e.g.* [](https://doi.org/b9k6dv) or [](https://doi.org/10.1021/jp972280j)). This can help to focus computational resources on regions of interest and improve sampling efficiency.
```

### 4. Sampling

Molecular dynamics or Monte Carlo simulations are performed within each window, allowing the system to equilibrate and sample configurations consistent with the biasing potential. In general ensuring that a given window has converged is not necessarily straightforward, but can be done by techniques such as [block averaging](https://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html).

### 5. Reweighting and combining the data

In the final step we gather the data from each window and combine it together to calculate the unbiased quantities of interest. I will first show how to unbias the data from each window, and then how to join all the results together.

In analogy with Eq. [](#marginal_P) we can defined a *biased* marginal probability density for the $i$-th windows, $P^b_i(\xi)$, as

$$
P^b_i(\xi) = \frac{Q^b_i(\xi)}{\int_{\xi_{\rm min}}^{\xi_{\rm max}} Q^b_i(\xi) d\xi} = \frac{\int_V e^{-\beta (H(\dofs) + V^{\rm bias}_i(\xi))} \delta(\xi - \xi(\dofs)) d\dofs}{\int_V e^{-\beta (H(\dofs) + V^{\rm bias}_i(\xi(\dofs))} d\dofs},
$$ (marginal_P_b)

where the biased partially-integrated partition function $Q^b_i(\xi)$ has also been defined. We note that the biasing factor $e^{-\beta V^{\rm bias}_i(\xi)}$ depends only on $\xi$ and therefore, since integration is performed on all degrees of freedom but $\xi$, can be moved outside of the integral. If we do so and then multiply and divide by $Q$ we obtain

$$
\begin{aligned}
P^b_i(\xi) & = e^{-\beta V^{\rm bias}_i(\xi)} \frac{\int_V e^{-\beta H(\dofs)} \delta(\xi - \xi(\dofs)) d\dofs}{\int_V e^{-\beta (H(\dofs) + V^{\rm bias}_i(\xi(\dofs))} d\dofs} \frac{Q_i}{Q_i}\\
& = e^{-\beta V^{\rm bias}_i(\xi)} \frac{\int_V e^{-\beta H(\dofs)} \delta(\xi - \xi(\dofs)) d\dofs}{Q} \frac{Q}{\int_V e^{-\beta (H(\dofs) + V^{\rm bias}_i(\xi(\dofs))} d\dofs}\\
& = e^{-\beta V^{\rm bias}_i(\xi)} P_i(\xi) \left\langle \frac{1}{e^{-\beta V^{\rm bias}_i(\xi)}} \right\rangle,
\end{aligned}
$$

where $\langle \cdot \rangle$ represents an *unbiased* ensemble average and $P_i(\xi)$ is the marginal probability density of the $i$-th window. Note that, being an ensemble average, $\left\langle \frac{1}{e^{-\beta V^{\rm bias}_i(\xi)}} \right\rangle$ does not depend on $\xi$, and therefore it is a (in general unknown) constant[^constant]. As a consequence, we can obtain the unbiased marginal probability density up to a multiplicative constant:

$$
\mathcal{P}_i(\xi) = P^b_i(\xi) e^{\beta V^{\rm bias}_i(\xi)} \propto P_i(\xi).
$$ (unbiasing)

where I use the symbol $\mathcal{P}_i(\xi)$ in place of $P_i(\xi)$, since the former is unnormalised and therefore not a proper probability density. This procedure is known as *unbiasing*, and applying it yields $N$ functions $\mathcal{P}_i(\xi)$ that are shifted relative to each other because of the unknown constant. The *total* $\mathcal{P}(\xi)$ can be recoverd by stitching together all the $\mathcal{P}_i(\xi)$, utilising the regions of the $\xi$-space where each pair of windows overlap significantly to find the unknown multiplying constants. There are several methods available to perform this task. Here I will present two such methods: a simple least-squares fit and the (much more powerful) WHAM.

```{attention} On the discrete nature of $P_i(\xi)$
The derivation above has been carried out by considering continuous functions for the sake of clarify. However, the simulation output is always a *histogram*, *i.e.* $P^b_i(\xi_k)$ (and, equivalently, $\mathcal{P}_i(\xi_k)$), where $\xi_k$ is a equispaced discrete variable. In the following derivations I will use this latter notation.
```

#### Least-squares method

Consider two windows $i$ and $j$ (with $|j - i| = 1$), whose unnormalised marginal probability densities overlap in a $\xi$-region $\lbrace \xi_o \rbrace$. We want to find the constant $C_{ij}$ that, multiplying $\mathcal{P}_j(\xi_k)$, minimises the mean-squared error between the two overlapping portions of the histograms, which is defined as

$$
{\rm MSE}_{ij} = \sum_{\zeta \in {\xi_o}} \left(\mathcal{P}_i(\zeta) - C_{ij} \mathcal{P}_j(\zeta) \right)^2.
$$

Imposing $\frac{d {\rm MSE}_{ij}}{d c_{ij}} = 0$ we find

$$
C_{ij} = \frac{\sum_{\zeta \in {\xi_o}} \mathcal{P}_i(\zeta) \mathcal{P}_j(\zeta)}{\sum_{\zeta \in {\xi_o}} \mathcal{P}_j^2(\zeta)}.
$$ (least-squares)

```{tip} Exercise
From a numerical point of view, it is often better to stitch the free-energy profiles $F_i(\xi) = -k_B T \ln \mathcal{P}_i(\xi)$ rather than the bare $\mathcal{P}_i(\xi)$, since in the former case the unknown constant is additive rather than multiplicative. Try to derive an expression for such an additive constant $A_{ij}$.
```

In practice, with this method the $0$-th window data are unchanged, while all the subsequent ones are rescaled one after the other by repeatedly applying Eq. [](#least-squares). Once the final $\mathcal{P}(\xi_k)$ is obtained, we can either normalise it (if we need a proper probability density), or used to compute the associated free-energy profile

$$
F(\xi) = -k_B T \ln \mathcal{P}(\xi) + {\rm const},
$$

where I made explicit the fact that classical free energies are always specified up to an additive constant, which can be chosen freely. If we are interested in a reaction between states $A$ and $B$, it is common to set $F(\xi_A)$ or $F(\xi_B)$ to 0, while for potentials of mean force (see for instance the [example](#us-example) below) it is customary to set $\lim_{\xi \to \infty} F(\xi) = 0$.

#### The Weighted Histogram Analysis Method (WHAM)

[WHAM](https://doi.org/10.1002/jcc.540130812) is a widely used reweighting technique for combining data from multiple biased simulations to obtain an unbiased estimate of the free energy profile. The basic idea behind WHAM is to reweight the probability distributions obtained from each window simulation such that they are consistent with each other and with the unbiased distribution. The reweighting process involves applying a set of equations that account for the biasing potentials applied in each window and the overlap between adjacent windows. WHAM simultaneously solves a set of self-consistent equations to iteratively refine the estimates of the unbiased probability distribution and the corresponding free energy profile.

[^constant]: This constant will take different values in different windows, since it is an ensemble average of a window-dependent quantity, $V^{\rm bias}_i(\xi)$.

(us-example)=
### A real-world example

As discussed in the chapter on [coarse-grained force fields](./coarse_grained.md), the effective interaction between two objects composed of multiple interacting units (atoms, molecules or beads) can be estimated in the dilute limit (*i.e.* at low density) as 

$$
U_{\rm eff}(R) = -k_B T \ln g(R),
$$ (umbrella_example_U_eff)

where $R$ is the distance between the two objects (defined *e.g.* as the distance between the two centres of mass) and $g(R)$ is the associated radial distribution function. For complicated objects composed of many parts, estimating $g(R)$ with sufficient accuracy through means of unbiased simulations requires an ungodly amount of computation time and is therefore unfeasible. Here I show an example where this issue has been overcome by using umbrella sampling.

Using the language introduced in this section, $R$ is the reaction coordinate and $g(R) = P(R) / 4 \pi R^2$ the observable of interest, where $P(R)$ is the marginal probability density.

```{figure}
:name: umbrella_example
:align: center

![A simulation snapshot showing the two chains (coloured differently).](figures/U_eff_conf.png)
![The biased radial distribution functions, with each colour corresponding to a different simulation.](figures/gr_biased.png)
![The unbiased radial distribution functions, with each colour corresponding to a different simulation.](figures/gr_unbiased.png)
![The reconstructed effective interaction](figures/U_eff.png)

The results of umbrella sampling simulations: the raw data is unbiased and then combined together to yield the final free-energy profile. Here the reaction coordinate $R$ is the distance between the centres of mass of two polymers.
```

[](#umbrella_example) shows the results of umbrella sampling simulations of a system composed of two polymer chains, where the chosen reaction coordinate is the distance between the two centres of mass, $R$, and the final output is the effective chain-chain interaction as a function of $R$. [](#umbrella_example-a) shows a snapshot of the two chains, [](#umbrella_example-b) shows the raw (biased) $g^b_i(R)$ data for all the windows $i$, and [](#umbrella_example-c) contains the $g^u_i(R), unbiased according to Eq. [](#unbiasing). Finally, [](#umbrella_example-d) contains the effective interaction, obtained with the WHAM method and shifted so that it vanishes at large distances.



[^umbrella]: The bias often takes the form of a harmonic potential whose shape, resembling an umbrella, gives the method its name.

## Thermodynamic & Hamiltonian integration

## Metadynamics

## Forward-flux sampling


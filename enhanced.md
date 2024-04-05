---
title: Enhanced sampling
---

```{note}
In this chapter I will use *reaction* as a general term for a microscopic process that transforms the system of interest (be it a macromolecule, an extended many-body system, a collection of atoms, *etc.*) between two well-defined states, $A$ and $B$.
```

## Reactions and rare events

Consider a system that can switch, possibly reversibly, between two macrostates, $A$ and $B$. Here the term macrostate is used loosely to indicate ensembles of microstates where the system resides for times that are much larger than the microscopic characteristic time; in thermodynamic parlance, $A$ and $B$, which are sometimes called *basins*, should be either metastable or equilibrium states, and therefore separated by a free-energy barrier $\Delta F_b$ larger than the thermal energy.

Examples relevant to computational biophysics are processes involving protein folding and unfolding, nucleic acid hybridisation, or switching between different conformations of the same (macro)molecule.

In this context the free-energy barrier between $A$ from $B$, $\Delta F_b^{A \to B} = F_{\rm max} - F_A$, is defined as the difference between the free energy of $A$, $F_A$ and that of the transition state, $F_{\rm max}$, which is the highest free-energy point along the reaction pathway[^reaction_pathway] connecting $A$ to $B$. Note that, per this definition, $\Delta F_b^{A \to B} \neq \Delta F_b^{B \to A} = F_{\rm max} - F_B$. See [](#fe_barrier) for a graphical definition of these quantities.From a kinetic perspective, the free-energy barrier influences the rate of the reaction, which is proportional to $e^{-\beta \Delta f_b}$.


```{figure} figures/fe_barrier.png
:name: fe_barrier 
:align: center

An example of the free energy landscape of a system displaying two basins separated by a transition state.
```

It is often the case that what interests us is the *reaction* itself rather than the $A$ and $B$ states, which are often known (and possibly have been characterised) beforehand. In this case, simulations starting from one of the two states, say $A$, would waste a lot of time in $A$, then quickly jump to state $B$, where it would again reside for a long time before switching basin once again, and so on. Therefore, using unbiased simulations[^unbiased] to sample the transition itself, for instance to evaluate the free-energy landscape as in [](#fe_barrier), requires a large computational effort which is mostly wasted in sampling uninteresting parts of the phase space.

In this part we will understand how the sampling of the transition from $A$ to $B$, a so-called [rare event](https://en.wikipedia.org/wiki/Rare_event_sampling), can be enhanced by using advanced computational techniques.

[^reaction_pathway]: See [below](#reaction-coordinates) for the precise meaning of "reaction pathway".
[^unbiased]: In unbiased simulations the system evolves according to the fundamental laws of physics (Newton's laws of motion in this context), without any artificial bias or constraints imposed.

(reaction-coordinates)=
## Reaction coordinates

Before we delve into the methods themselves, we need to introduce a framework that makes it possible to describe the reaction we are interested in. We do so by defining one or more *reaction coordinates* (also known as *order parameters* in physics lingo). In the context of molecular simulations, a reaction coordinate is a parameter or set of parameters that describe the progress of a given reaction. It serves as a way to quantify and track the changes occurring within a system as it evolves through the reaction pathway. The primary purpose of a reaction coordinate is to provide a simplified description of the system's thermodynamics, thus making it possible to monitor and analyze the progress of a reaction in terms of a single or a few variables rather than the whole multidimensional phase space. This simplification is essential for understanding the microscopic underpinnings of the reaction of interest. Defining a reaction coordinate makes it possible to draw a diagram such as the one shown in [](#fe_barrier), which are often called free-energy profiles or landscapes, where the variation of the free energy along a particular reaction coordinate or collective variable is plotted.

The choice of reaction coordinate depends on the specific process being studied. It could be a simple geometric parameter such as bond length, bond angle, or dihedral angle, or it could be a more complex collective variable that captures the overall structural changes in the system, such as the distance between two key functional groups, the position or coordination number of a particular atom, or the solvent-accessible surface area of a biomolecule. Here *collective variable* means a quantity that is a function of multiple degrees of freedom. As a simple example, consider a system formed by two DNA strands made of $N_1$ and $N_2$ units (atoms or coarse-grained beads), respectively, that can hybridise: a possible (although not exactly ideal) reaction coordinate could be the distance between the two strands' centres of mass, *viz.*:

$$
\xi = \left| \vec{R}^{\rm cm}_1 - \vec{R}^{\rm cm}_2 \right| = \left| \frac{1}{N_1} \sum_{i \in N_1} \vec{r}_i - \frac{1}{N_2} \sum_{i \in N_2} \vec{r}_i \right|
$$

```{note}
In general, a reaction coordinate can be multidimensional, and sometimes it is useful, or even necessary, to define such complex variables. However, for the sake of simplicity here I will always try to use unidimensional reaction coordinates, which will be designated by the symbol $\xi$.
```

## Umbrella Sampling

## Thermodynamic & Hamiltonian integration

## Metadynamics

## Forward-flux sampling


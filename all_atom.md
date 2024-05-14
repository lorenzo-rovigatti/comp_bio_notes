---
title: Classical molecular dynamics and all-atom simulations
exports:
    - format: pdf
      template: plain_latex
---

Quantum mechanics provides a rigorous framework for describing the behavior of molecules that explicitly includes the effect of the electrons. However, the calculations are very heavy and thus the evolution of a system can be followed for short times only. In addition, the computational complexity of simulating quantum systems often scales super-linearly with the size of the problem[^scaling]. This poses significant challenges when attempting to model large-scale many-body systems accurately and for long times.

One of the fundamental distinctions between quantum and classical approaches lies in the treatment of electron contributions. In quantum calculations, the interactions and motions of electrons are computed explicitly, accounting for their wave-like nature and the intricacies of their interactions with nuclei and other electrons. This level of detail enables precise predictions of electronic structure and properties but comes at a high computational cost, particularly as the number of electrons and/or atoms increases.

In contrast, classical force fields adopt a simplified approach where the contributions of electrons are averaged out. Instead of explicitly modeling individual electron behaviors, classical force fields approximate the interactions between atoms using simplified mathematical models based on classical mechanics. While this approach sacrifices the quantum mechanical rigor of true wavefunction-based simulations, it significantly reduces computational complexity, making it feasible to simulate large-scale many-body systems over longer timescales and larger spatial dimensions.

[^scaling]: The complexity ranges from $\mathcal{O}(e^N)$ for brute-force implementations, to $mathcal{N^3}$ for many DFT codes, but can be linear in some cases (see *e.g.* [](doi:10.1088/0034-4885/75/3/036503)).

# From quantum to classical mechanics

```{warning}
TODO
```

# Classical force fields

Classical (or empirical) force fields are foundational tools in computational chemistry and physics, providing a versatile and efficient framework for simulating the behavior of systems at the atomic and molecular scales. At their core, classical force fields describe the interactions between atoms using simplified mathematical models based on quantum mechanics calculations, empirical observations, and physical principles from classical mechanics. Unlike quantum methods, which rigorously account for the wave-like nature of particles and their interactions through principles like superposition and entanglement, classical force fields offer a pragmatic approach that balances computational feasibility with accuracy.

The central concept behind classical force fields is the representation of interatomic interactions through potential energy functions, which quantify the energy associated with specific configurations of atoms and molecules. These potential energy functions typically consist of terms that describe various types of interactions, including covalent bonds, van der Waals (dispersion) forces and electrostatic interactions. By summing these contributions over all pairs or groups of atoms, classical force fields make it possible to evaluate the total potential energy of a system, as well its derivatives with respect to the system's degrees of freedom, which can be used to investigate its behaviour and dynamics.

One of the key advantages of classical force fields lies in their computational efficiency. Indeed, unlike quantum simulations, which require solving the Schrödinger equation or approximations thereof for each configuration of the system, classical force fields involve straightforward mathematical operations that scale linearly[^classical_scaling] with the number of particles. Despite their simplicity, classical force fields can exhibit remarkable predictive power and have been successfully applied across diverse scientific disciplines.

In a classical force field (FF from now on), the main assumption is that of *additivity*, which makes it possible to write down the total energy of a system as a sum of several contributions. These are classified into terms that account for atoms that are linked by covalent bonds (*bonded* or *local* terms), and terms that account for noncovalent interactions (*non-bonded* or *non-local* terms), so that the total energy of the system is written as

$$
E_{\rm tot}(\dofs) = E_{\rm bonded}(\dofs) + E_{\rm non-bonded}(\dofs),
$$ (eq:tot-energy)

where I made explicit the fact that the energy (and its contributions) are functions of the set of positions of the atoms, $\dofs$.

```{note}
The separation of contributions into bonded and non-bonded terms is not only a conceptual one, but it has also practical implications. For instance, bonded interactions act on pairs or groups of atoms that do not change during the course of the simulation, and therefore the data structures that keep track of these interacting atoms do not need to be updated. Moreover, non-bonded interactions are usually much less steep than bonded ones, which makes it possible to implement multiple-timestep protocols to speed-up simulations (see *e.g.* Chapter 14 of @schlick2010molecular).
```

The functional forms of the functions (implicitly) defined in Eq. [](#eq:tot-energy) depend on the force field employed, and are parametrised so as to reproduce the experimental structure and dynamics of target molecular systems.

[^classical_scaling]: As discussed below, linear scaling in systems interacting through *long-range interactions* such as the Coulomb potential require sophisticated techniques such as [Fast Multipole Method](https://doi.org/10.1016/0021-9991(87)90140-9). A more easy-to-achieve scaling is $\mathcal{O}(N \log N)$.

## Bonded interactions

In the context of molecular systems, *bonded interactions* refer to the forces that arise between atoms due to the formation of chemical bonds. These interactions occur when atoms are directly connected to each other through covalent bonds, which involve the sharing of electrons.

```{caution}
Note that the term *bonded* is often used for other mechanism through which atoms or molecules pair beyond covalend bonding (*e.g.* hydrogen bonding). It should be clear from the context what its implied meaning is, but be cautious!
```

Bonded interactions typically include several contributions due to bond stretching, angle bending, dihedral and improper torsions, so that the total bonded energy can be written as

$$
E_{\rm bonded} = \sum_{b \in \rm bonds} E_{\rm stretch}(b) + \sum_{\theta \in \rm angles} E_{\rm bend}(\theta) + \sum_{\varphi \in \rm dihedrals} E_{\rm tors}(\varphi),
$$ (eq:bonded-energy)

where the three sums run over all the covalent bonds, all bond angles[^bond_angle], and proper and improper dihedral angles[^dihedral], respectively.

Note that in writing Eq. [](#eq:bonded-energy) we have made an additional assumption on the interactions, which is commonly (but not always) verified in many FFs. Namely, there are no *cross terms* that couple different internal coordinates (such as bond lengths and angles). We will take a cursory look at these terms [below](#sec:cross_terms).

[^bond_angle]: A bond angle is the angle between two bonds that include a common atom.
[^dihedral]: Angles defined by groups of four atoms that are connected by covalent bonds.

### Bond stretching

This interaction occurs when atoms connected by a covalent bond move closer together or farther apart, resulting in the stretching or compression of the bond. The energy associated with bond stretching is described by a potential energy function that typically resembles a harmonic oscillator potential, with the energy increasing quadratically with bond elongation or compression.

### Angle bending

When three atoms are connected by two consecutive covalent bonds, the angle between these bonds can change, causing the atoms to deviate from a linear configuration. Angle bending interactions describe the energy associated with such deformations and are often modeled using harmonic potentials, where the energy increases quadratically with the deviation of the angle from its equilibrium value.

### Dihedral rotation

Dihedral interactions arise when four atoms are connected by three consecutive covalent bonds, forming a torsion angle or dihedral angle. Rotating one group of atoms around the axis defined by the central bond changes the conformation of the molecule, leading to different energy minima corresponding to different dihedral angles. Dihedral interactions are typically described using periodic potentials that capture the periodicity of the energy as a function of the dihedral angle.

### Improper rotation

Improper torsion interactions occur when four atoms are bonded in a way that one atom is not directly in line with the other three. This configuration is often necessary to maintain the correct geometry or chirality of a molecule. Improper torsions are used to model the energy associated with deviations from the ideal geometry. The potential energy associated with improper torsions is typically described by a harmonic potential or a cosine function, depending on the force field parameters and the desired behavior of the molecule.

(sec:cross_terms)=
### Cross terms

## Non-bonded interactions

*Non-bonded interactions* refer to the energy contributions that arise between atoms or molecules that are not directly connected by chemical bonds, and are typically categorized into two main types:

### Van der Waals Interactions

Van der Waals interactions are weak forces that arise due to fluctuations in the electron density of atoms or molecules. These interactions include both attractive forces, arising from dipole-dipole interactions and induced dipole-induced dipole interactions (van der Waals dispersion forces, see @israelachvili2011intermolecular for a derivation of these terms), and repulsive forces, resulting from the overlap of electron clouds at close distances. Van der Waals interactions are described by empirical potential energy functions, such as the Lennard-Jones potential, which accounts for both attractive and repulsive components.

### Electrostatic Interactions

Electrostatic interactions arise from the attraction or repulsion between charged particles, such as ions or polar molecules. These interactions are governed by Coulomb's law, which states that the force between two charged particles is proportional to the product of their charges and inversely proportional to the square of the distance between them. In molecular systems, electrostatic interactions include interactions between charged atoms or molecules (ion-ion interactions), as well as interactions between charged and neutral atoms or molecules (ion-dipole interactions or dipole-dipole interactions).

### Hydrogen Bonding

Represents the specific type of non-covalent interaction between a hydrogen atom covalently bonded to an electronegative atom (e.g., nitrogen, oxygen, or fluorine) and a nearby electronegative atom. Hydrogen bonding interactions are often modeled using empirical potentials that account for the directionality and strength of hydrogen bonds.

## Transferability

```{warning}
TODO
```

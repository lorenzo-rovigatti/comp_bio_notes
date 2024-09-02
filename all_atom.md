---
title: Classical molecular dynamics and all-atom simulations
exports:
   - format: pdf
---

```{tip}
The main references for this part are @frenkel2023understanding, @schlick2010molecular and @leach2001molecular.
```

Quantum mechanics provides a rigorous framework for describing the behavior of molecules that explicitly includes the effect of the electrons. However, the calculations are very heavy and thus the evolution of a system can be followed for short times only. In addition, the computational complexity of simulating quantum systems often scales super-linearly with the size of the problem[^scaling]. This poses significant challenges when attempting to model large-scale many-body systems accurately and for long times.

One of the fundamental distinctions between quantum and classical approaches lies in the treatment of electron contributions. In quantum calculations, the interactions and motions of electrons are computed explicitly, for instance with DFT methods, as discussed in the previous Chapter. This level of detail enables precise predictions of electronic structure and properties but comes at a high computational cost, particularly as the number of electrons and/or atoms increases.

In contrast, classical force fields adopt a simplified approach where the contributions of electrons are averaged out. Instead of explicitly modeling individual electron behaviors, classical force fields approximate the interactions between atoms using simplified mathematical models based on classical mechanics. Note that there exist also "mixed" methods, where some degrees of freedom are accounted for by using a quantum-mechanical treatment, while others are treated classically, *e.g.* the nuclear degrees of freedom in the CPMD approach.

[^scaling]: The complexity ranges from $\mathcal{O}(e^N)$ for brute-force implementations, to $\mathcal{N^3}$ for many DFT codes, but can be linear in some cases (see *e.g.* [](doi:10.1088/0034-4885/75/3/036503)).

# Molecular dynamics

In a molecular dynamics simulation, the atoms, or particles, follow Newton's equations of motions. As a result, the dynamics happens on a hypersurface defined by $E = \text{const}$, where $E$ is the energy of the system. In practice, the equations of motions are solved iteratively by discretising time: the quantities of interest (position $r$, velocity $v$ and force $F$) at time $t$ are used to obtain those at time $t + \Delta t$, where $\Delta t$ is the integration *time step*. The flow of a simple MD program is:

1. The simulation parameters (*e.g.* temperature, density, time step) are read and initialised.
2. The initial configuration (*i.e.* the initial positions and velocities) is read or generated.
3. The simulation runs, iterating the following steps:
   1. The forces (and possibly torques) on all particles are computed.
   2. The positions and velocities are updated according to Newton's equations.
4. The simulation stops when some condition is met (*e.g.* number of steps run).

The following snippet shows a pseudo-code implementation of foregoing algorithm:

:::{code} pseudocode
:label: code:MD_simple
:caption: Pseudo-code for a simple molecular dynamics code.

CALL initialise_parameters()
CALL initialise_system()

WHILE t is smaller than t_max
   CALL compute_force()
   CALL integrate()
   t += delta_t
   CALL sample_observables()
:::

## The force calculation

This is the most time-consuming part of an MD simulation. It consists of evaluating the force (and, for rigid bodies, the torque) acting on each individual particle. Since, in principle, every particle can interact with every other particle, for a system of $N$ objects the number of interacting contributions will be equal to the number of pairs, $N (N - 1) / 2$. Therefore, the time required to evaluate all forces scales as $\mathcal{O}(N^2)$. Although, as we will see later, there are methods that can bring this complexity down to $\mathcal{O}(N)$, here I present the simplest case.

## Integrating the equations of motion

:::{important}
For the sake of simplicity, I will derive equations for 1D systems, but these remain valid for 2D and 3D systems, as the dynamics is decoupled along the three dimensions.
:::

Consider a particle moving along the $ x $-axis under the influence of a force $F(t)$. Newton's second law gives $m \frac{d^2 x(t)}{dt^2} = F(t)
$, so that the acceleration $ a(t) $ is $a(t) = \frac{F(t)}{m}$

Expanding the position $ x(t) $ around the time $t$ using a Taylor series we obtain

\begin{align}
x(t + \Delta t) = x(t) + v(t) \Delta t + \frac{1}{2} a(t) \Delta t^2 + \frac{1}{6} \frac{da(t)}{dt} \Delta t^3 + \mathcal{O}(\Delta t^4)\\
x(t - \Delta t) = x(t) - v(t) \Delta t + \frac{1}{2} a(t) \Delta t^2 - \frac{1}{6} \frac{da(t)}{dt} \Delta t^3 + \mathcal{O}(\Delta t^4)
\end{align}

Adding the two Taylor expansions yields

$$
x(t + \Delta t) + x(t - \Delta t) = 2x(t) + a(t) \Delta t^2 + \mathcal{O}(\Delta t^4),
$$

which, if we neglect the higher-order terms $ \mathcal{O}((\Delta t)^4)$ becomes

$$
x(t + \Delta t) = 2x(t) - x(t - \Delta t) + a(t) \Delta t^2.
$$

This is the Verlet algorithm, which allows us to calculate the position $x(t + \Delta t)$ at the next time step using the current position $x(t)$, the previous position $ x(t - \Delta t)$, and the current acceleration $a(t)$. Although this basic Verlet method does not explicitly involve velocity, we can compute it as

$$
v(t) = \frac{x(t + \Delta t) - x(t - \Delta t)}{2\Delta t} + \mathcal{O}(\Delta t^2).
$$

Note that this is accurate only up to order $\Delta t^2$. We can modify the basic Verlet approach to include an explicit update for the velocity, leading to the Velocity Verlet method. Instead of relying on the positions from the previous and current time steps, the Velocity Verlet algorithm updates the position and velocity in a two-step process.

First, we use the current velocity and acceleration to update the position at time $ t + \Delta t $. This is done similarly to the basic Verlet method but with the velocity term explicitly included:

$$
x(t + \Delta t) = x(t) + v(t) \Delta t + \frac{1}{2} a(t) \Delta t^2
$$ (eq:velocity_verlet_x)

This equation uses the current position $ x(t) $, the current velocity $ v(t) $, and the current acceleration $ a(t) $ to compute the new position $ x(t + \Delta t) $. Next, after updating the position, we need to compute the new acceleration at time $ t + \Delta t $ because the force (and hence the acceleration) may have changed due to the updated position. The new acceleration is given by:

$$
a(t + \Delta t) = \frac{F(t + \Delta t)}{m}
$$

With this new acceleration in hand, we can update the velocity. Instead of using just the current acceleration, the Velocity Verlet method uses the average of the current and new accelerations to update the velocity:

$$
v(t + \Delta t) = v(t) + \frac{1}{2} \left[ a(t) + a(t + \Delta t) \right] \Delta t
$$ (eq:velocity_verlet_v)

This velocity update equation accounts for the change in acceleration over the time step, providing a more accurate velocity update than simply using the current acceleration. The Velocity Verlet method is the *de-facto* standard for MD codes. The common way to implement it is to split the velocity integration step in two, so that a full MD interation becomes

1. Update the velocity, first step, $v(t + \Delta t / 2) = v(t) + \frac{1}{2} a(t) \Delta t$.
2. Update the position, $r(t + \Delta t) = r(t) + v(t + \Delta t / 2)\Delta t = r(t) + v(t) \Delta t + \frac{1}{2} a(t) \Delta t^2$ (*e.g.* Eq. [](#eq:velocity_verlet_x)).
3. Calculate the force (and therefore the acceleration) using the new position, $r(t + \Delta t) \to a(t + \Delta t) = F(t + \Delta t) / m$.
4. Update the velocity, second step, $v(t + \Delta t) = v(t + \Delta t / 2) + \frac{1}{2} a(t + \Delta t) = v(t) + \frac{1}{2} \left[ a(t) + a(t + \Delta t) \right] \Delta t$ (*e.g.* Eq. [](#eq:velocity_verlet_v)).

An MD implementation based on [](#code:MD_simple) that implements a Velocity Verlet algorithm is

:::{code} pseudocode
:label: code:MD_velocity_verlet
:caption: Pseudo-code for an MD code implementing the Velocity Verlet algorithm.

CALL initialise_parameters()
CALL initialise_system()

WHILE t is smaller than t_max
   CALL integrate_first_step()
   CALL compute_force()
   CALL integrate_second_step()
   
   t += delta_t
   CALL sample_observables()
:::

## Reduced units

## Energy conservation

## Thermostats

(sec:andersen_thermostat)=
### Andersen

### Nose-Hoover

### Bussi-Donadio-Parrinello

## Barostats

## Tricks of the trade

(sec:neighbour_lists)=
### Neighbour lists

### Long-range interactions

## Some observables

(sec:compute_pressure)=
### Pressure

### Radial distribution function

### Mean-squared displacement

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

Dihedral or torsional interactions arise when four atoms are connected by three consecutive covalent bonds, forming a torsion angle or dihedral angle. Rotating one group of atoms around the axis defined by the central bond changes the conformation of the molecule, leading to different energy minima corresponding to different dihedral angles.

As noted [before](./proteins.md#sec:molecular_vibrations), in proteins the most important dihedral (torsional) angles are those associated to rotations around the $N - C^\alpha$ and $C^\alpha - C$ bonds, $\phi$ and $\psi$, which feature free-energy barriers of the order of the thermal energy. However, rotations around the peptide bond (dihedral $\omega$) and around sp$^3$-sp$^3$ bonds such as those found in aliphatic said chains (dihedral $\chi$) can also play a role in dictating the flexibility of proteins.

Dihedral interactions are typically described using periodic potentials that capture the periodicity of the energy as a function of the dihedral angle. A generic torsional interaction associated to a dihedral angle $\varphi$ takes the form

$$
E_t(\varphi) = \sum_{n} \frac{V_n}{2}[1 + \cos(n\varphi - \varphi_0)],
$$

where $n$ is an integer that determines the periodicity of the barrier of height $V_n$, and $\varphi_0$ is a reference angle, which is often set to $0$ or $\pi$.

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

## GROMACS

GROMACS, an acronym for Groningen Machine for Chemical Simulations, is an open-source software package designed primarily for molecular dynamics simulations of biomolecular systems, providing insights into the structural dynamics of proteins, lipids, nucleic acids, and other complex molecular assemblies. GROMACS can run simulations efficiently on a wide range of hardware platforms, from single processors to large parallel computing clusters, and implements advanced algorithms and techniques, such as domain decomposition, particle-mesh Ewald summation for long-range electrostatics, and multiple time-stepping schemes. 

Interestingly (and differently from other simulation engines) GROMACS also offers a comprehensive suite of analysis tools, making it possible to perform tasks such as trajectory analysis, free energy calculations, and visualization of simulation results.

* Use `mdp` files to create a `tpr` (which is a binary, non-human-readable file): at this step we need a topology file (`topol.top`)[^topol_file]
* Index files (`.ndx`) are used to assign atoms to categories that can be used to easily find them. Use `gmx make_ndx -f file.tpr -o output` to make a new index file.
* The `-deffnm` flag tells Gromacs to use the input file as a template for the filenames of output quantities
* The default trajectory file (`.trr`) contains all the info (coordinates, velocities, *etc.*). By contrast, `xtc` files are lighter since they store only the atom coordinates.

```{warning} Berendsen thermostat
On newer versions, Gromacs will spit out a warning if you use the Berendsen thermostat (`tcoupl = berendsen`). However, by default `gmx grompp` considers a single warning as a fatal error. Use the `-maxwarn` flag to raise the number of acceptable warnings (*e.g.* `-maxwarn `).
```

[^topol_file]: in Gromacs lingo, a topology file contains the details about the force field

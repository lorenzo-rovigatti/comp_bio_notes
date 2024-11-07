---
title: Molecular quantum mechanics
exports:
   - format: pdf
---

```{warning}
This is the part that is the farthest from my own expertise, which is about classical simulations. Use this chapter with caution!
```

```{tip}
The main references for this part are @giustino2014materials and @bottcher2021computational.
```

In this chapter I briefly present the background required to understand the basic algorithms used to run *ab-initio* molecular dynamics (AIMD) simulations. The methods presented here makes it possible to simulate the dynamics of atoms and molecules from first principles, which is what "ab initio" means.

# A minimalistic introduction to quantum mechanics

Quantum mechanics is a theoretical framework that describes the behavior of matter and energy at the smallest scales. Unlike classical mechanics, where objects have definite positions and velocities, quantum mechanics introduces the concept of wave-particle duality, where particles exhibit both wave-like and particle-like properties.

The full information about the quantum state of a given quantum system is contained in a complex-valued function of the spatial and time coordinates $\vec{r}$ and $t$, called the wavefunction, $\Phi(t, \vec{r})$. Max Born's interpretation of the wavefunction provides a probabilistic view of quantum mechanics. According to Born's rule, the probability of finding a particle at position $\vec{r}$ at time $t$ is proportional to $|\Phi(t, \vec{r})|^2$. Thus, the wavefunction must be normalized, meaning that the total probability of finding the particle anywhere in space must equal 1:

$$
\int_V |\Phi(t, \vec{r})|^2 d\vec{r} = \int_V |\Psi(\vec{r})|^2 d\vec{r} = 1,
$$ (eq:QM_normalisation)

where we used the factorisation $\Phi(t,\vec{r}) = \phi(t) \Psi(\vec{r})$ which will be justified in a moment, and the integral is carried out over the system's volume, $V$. This probabilistic interpretation is one of the fundamental departures from classical mechanics, where particles have deterministic trajectories. Note that, in order for Eq. [](#eq:QM_normalisation) to make sense, wave functions need to be "square-integrable" or, in other words, to belong to the Hilbert space $L^2(\mathbb{R}^3)$. This means that it is also possible to define an inner product between any two wave functions $\phi$ and $\psi$:

$$
\langle \phi | \psi \rangle = \int \phi^*(\vec r) \psi(\vec r) d\vec r,
$$

where $\phi^*(\vec{r})$ is the complex conjugate of $\phi(\vec r)$, and we also introduced the bra-ket notation.

In quantum mechanics, physical quantities (observables) such as energy, momentum, and position are represented by operators that act on the wavefunction. For example, the position operator is $\hat{\vec{r}} = \vec{r}$, and the momentum operator, is $\hat{\vec{p}} = -i\hbar \nabla$, where $\hbar \equiv h / 2 \pi$ is the reduced Planck constant. The expectation value of an observable over a wavefunction $\Psi$, which is the average value measured in an experiment, is given by:

$$
\langle \hat{A} \rangle = \langle \Psi | \hat A | \Psi \rangle = \int \Psi^*(\vec{r}) \hat{A} \Psi(\vec{r}) d\vec{r},
$$

where $\hat{A}$ is the operator corresponding to the observable.

The *Hamiltonian operator*, or simply *Hamiltonian*, of a system, $\hat{H}$, represents the total energy of the system (kinetic and potential) and determines its evolution. The Hamiltonian typically takes the form:

$$
\hat{H} = \frac{\hat p^2}{2m} + V(\vec{r}) = -\frac{\hbar^2}{2m} \nabla^2 + V(\vec{r})
$$ (eq:hamiltonian)

where $m$ is the mass of the particle, $\hat p =  -i\hbar \nabla$ is the particle's momentum, $\nabla^2$ is the Laplacian operator, representing the kinetic energy contribution, and $V(\vec{r})$ is the potential energy as a function of position. For a system of multiple particles, the Hamiltonian becomes more complex, accounting for the kinetic and potential energies of all particles as well as their interactions.

Given a Hamiltonian, the fundamental equation that determines how the wavefunction of a system evolves over time is the time-dependent Schrödinger equation:

$$
i\hbar \frac{d \Phi(t, \vec{r})}{dt} = \hat H \Phi(t, \vec{r}).
$$ (eq:schroedinger_time)

If $\hat H$ does not depend explicitly on time, the dependence on time and spatial variables of the wavefunction can be separated, *i.e.* $\Phi(t, \vec{r}) = \Psi(\vec{r}) \phi(t)$. The time-dependence of Eq. [](#eq:schroedinger_time) can be explicitly integrated, yielding

$$
\phi(t) = e^{-iEt / \hbar},
$$

where $E$ is the expectation value of the Hamiltonian and will be further discussed below. The time-independent part of the Schrödinger equation is instead

$$
\hat{H} \Psi(\vec{r}) = E \Psi(\vec{r}).
$$ (eq:schroedinger)

This is a partial differential eigenvalue equation, where an operator acts on a function, called eigenfunction, and returns the same function multiplied by the eigenvalue associated to the eigenfunction. Here the eigenfunction is the wavefunction, and the eigenvalue is its energy, $E$. Since the kinetic part of the Hamiltonian always contains $\nabla^2$ terms, the Schrödinger equation is a second-order differential equation.

## The many-body problem

For systems with multiple interacting particles (*e.g.*, an atom, a molecule or larger objects), the Schrödinger equation becomes a many-body problem. The Hamiltonian for such a system includes terms for the kinetic energy of each particle, the potential energy due to external fields, and the interaction energy between the particles. For instance, for a system of $N$ charged particles, the many-body Hamiltonian can be written as:

$$
\hat{H} = \sum_{i=1}^{N} \left( -\frac{\hbar^2}{2m_i} \nabla_i^2 + V(\vec{r}_i) \right) + \frac{1}{2} \sum_{i} \sum_{i \neq j} \frac{q_i q_j}{|\vec{r}_i - \vec{r}_j|},
$$ (eq:many_body_H)

where $m_i$ and $q_i$ are the mass and charge of the $i$-th particle, respectively. Here the first term represents the kinetic energy of the particles, the second term represents the potential energy due to an external field (*e.g.*, from atomic nuclei if we are considering electrons), and the third term represents the Coulomb interaction. Here we set the Coulomb constant $1/4\pi\epsilon_0$ to 1 to simplify the notation.

```{note} Atomic units
A common convention is to use atomic units (a.u.), where the reduced Planck constant $\hbar$, the electron charge $e$, the electron mass $m_e$, and the Coulomb constant are all set to 1. This simplifies many of the equations and is a common convention in computational quantum chemistry.
```

Solving the Schrödinger equation for many-body systems exactly is extremely challenging, if not impossible, for all but the simplest cases, such as the hydrogen atom. Therefore, approximations and numerical methods are essential. This is where *ab initio* methods like Density Functional Theory (DFT) come into play, offering a way to handle the complexities of interacting electrons.

(sec:born-oppenheimer)=
# The Born-Oppenheimer approximation

When dealing with molecules or solids, one faces a many-body problem involving both the nuclei and the electrons. The total Hamiltonian of a molecular system includes terms for the kinetic energy of the nuclei, the kinetic energy of the electrons, and the potential energy due to interactions between all particles (nuclei-nuclei, nuclei-electrons, and electron-electron interactions). This problem is, in general, analitically unsolvable due to the coupling between the electronic and nuclear degrees of freedom.

To make the problem tractable, we introduce the Born-Oppenheimer (BO) approximation. This approximation exploits the large mass difference between nuclei and electrons ($\approx 3$ orders of magnitude) to decouple their motions, making it possible to simplify the Schrödinger equation and compute electronic structures for fixed nuclear configurations.

For a molecule composed of $N$ electrons of mass $m_e$ and $M$ nuclei, each of mass $M_I$, the total Hamiltonian can be written as:

$$
\hat{H} = \hat{T}_\text{nuc} + \hat{T}_\text{el} + \hat{V}_\text{nuc-nuc} + \hat{V}_\text{nuc-el} + \hat{V}_\text{el-el},
$$ (eq:full_H)

where $\hat{T}_\text{nuc} = -\sum_{I=1}^{M} \frac{\hbar^2}{2M_I} \nabla_I^2$ and $\hat{T}_\text{el} = -\sum_{i=1}^{N} \frac{\hbar^2}{2m_e} \nabla_i^2$ are the kinetic energy operators for the nuclei and electrons, respectively, while $\hat{V}_\text{nuc-nuc}$, $\hat{V}_\text{nuc-el}$ and $\hat{V}_\text{el-el}$ represent the nucleus-nucleus, electron-nucleus, and electron-electron Coulomb interactions, respectively.

The total wavefunction $\Psi(\{\vec{R}_I\}, \{\vec{r}_i\})$ depends on both the positions of the nuclei $\{\vec{R}_I\}$ and of the electrons $\{\vec{r}_i\}$. Solving for the full wavefunction is complicated because the nuclear and electronic motions are coupled through the potential terms. However, the mass disparity between the nucleus and the electrons suggests that, on the timescale of nuclear motion, the electrons can adjust almost instantaneously to any change in nuclear positions. This allows us to decouple the nuclear and electronic motions by assuming that the electronic wavefunction can be solved for a fixed configuration of nuclei.

We assume that the total wavefunction $\Psi(\{\vec{R}_I\}, \{\vec{r}_i\})$ can be factored into a product of a nuclear wavefunction $\chi(\{\vec{R}_I\})$ and an electronic wavefunction $\psi(\{\vec{r}_i\}; \{\vec{R}_I\})$:

$$
\Psi(\{\vec{R}_I\}, \{\vec{r}_i\}) \approx \chi(\{\vec{R}_I\}) \psi(\{\vec{r}_i\}; \{\vec{R}_I\})
$$

Here, the electronic wavefunction $\psi(\{\vec{r}_i\}; \{\vec{R}_I\})$ is parameterized by the fixed nuclear positions $\{\vec{R}_I\}$, while $\chi(\{\vec{R}_I\})$ describes the motion of the nuclei.

## Solving the Schrödinger equation

With the nuclei fixed, we can solve the electronic Schrödinger equation for a given set of nuclear positions:

$$
\hat{H}_\text{el} \psi(\{\vec{r}_i\}; \{\vec{R}_I\}) = E_\text{el}(\{\vec{R}_I\}) \psi(\{\vec{r}_i\}; \{\vec{R}_I\}).
$$

The electronic Hamiltonian $\hat{H}_\text{el}$ includes the kinetic energy of the electrons and the potential energy due to electron-electron and electron-nuclei interactions:

$$
\hat{H}_\text{el} = \hat{T}_\text{el} + \hat{V}_\text{el-el} + \hat{V}_\text{nuc-el}.
$$ (eq:electronic_H)

The result of solving this equation is the electronic energy $E_\text{el}(\{\vec{R}_I\})$, which depends parametrically on the nuclear positions. This energy forms a potential energy surface (PES) on which the nuclei move.

Once the electronic problem is solved, the nuclear dynamics can be described by an effective nuclear Hamiltonian, where the potential energy surface $E_\text{el}(\{\vec{R}_I\})$ from the electronic calculation acts as the potential for the nuclei:

$$
\hat{H}_\text{nuc} \chi(\{\vec{R}_I\}) = \left( \hat{T}_\text{nuc} + \hat{V}_\text{nuc-nuc} + E_\text{el}(\{\vec{R}_I\}) \right) \chi(\{\vec{R}_I\})
$$

This equation governs the motion of the nuclei on the potential energy surface created by the electrons. In practice, depending on the temperature and energy scales, this motion can be treated either classically or quantum mechanically.

```{note} Adiabaticity
In quantum mechanics, a process is called *adiabatic* if the system remains in the same quantum state (or subspace) throughout the process, provided that the external conditions (like the positions of the nuclei) change sufficiently slowly. Specifically, for the BO approximation, this means that as the nuclei move, the electrons remain in their instantaneous ground state at all times. 

More formally, if the nuclei evolve adiabatically, the system stays in the ground-state wavefunction of the electronic Hamiltonian corresponding to the instantaneous positions of the nuclei: the nuclear motion is slow enough to allow the electrons to adjust without ever being excited to higher electronic states.

The assumption of adiabaticity—and therefore the validity of the BO approximation—can break down in certain situations, leading to non-adiabatic effects. These effects occur when the nuclear motion becomes fast enough that the electrons cannot adjust instantaneously. In such cases, the system may undergo transitions between different electronic states, leading to a mixing of electronic and nuclear dynamics.

Situations where non-adiabatic effects are significant include:

- **Conical intersections:** Points in the nuclear configuration space where two or more potential energy surfaces intersect. Near conical intersections, the energy gap between electronic states becomes very small, and even slow nuclear motion can lead to electronic transitions.
- **High-temperature or high-energy regimes:** When the nuclei move with enough kinetic energy, the assumption of slow nuclear motion relative to electronic adjustments may no longer hold.
- **Photochemical reactions:** When a molecule absorbs light, it can be excited to a higher electronic state. The subsequent dynamics often involve non-adiabatic transitions between electronic states as the molecule relaxes back to the ground state.

In these cases, the BO approximation fails, and more sophisticated methods, such as [time-dependent DFT](https://en.wikipedia.org/wiki/Time-dependent_density_functional_theory), are required to accurately describe the system.
```

## Implications and limitations

The Born-Oppenheimer approximation is widely used in *ab initio* simulations because it allows the separation of electronic structure calculations from nuclear motion, significantly reducing the computational complexity. It underpins most of the computational chemistry methods, including Hartree-Fock theory, post-Hartree-Fock methods, and density functional theory (DFT). However, the approximation has some limitations. For instance, in systems where nuclear and electronic motions are strongly coupled (*e.g.*, near conical intersections or in the presence of non-adiabatic effects, see box above), the Born-Oppenheimer approximation breaks down, and one must account for the simultaneous evolution of electrons and nuclei. Moreover, the approximation often treats nuclei as classical particles moving on a potential energy surface, which may neglect important quantum effects like tunneling or zero-point energy, particularly in light atoms like hydrogen.

# The Hohenberg-Kohn theorems

The Hohenberg-Kohn theorems form the foundation of Density Functional Theory (DFT), a widely used approach in *ab initio* simulations. DFT simplifies the quantum many-body problem by shifting the focus from the many-body wavefunction to the electron density $n(\vec{r})$. These theorems, established by [Pierre Hohenberg and Walter Kohn in 1964](doi:10.1103/PhysRev.136.B864), provide the theoretical justification for this approach, proving that all ground-state properties of a system are uniquely determined by its electron density. By leveraging this result, it is possible to reformulate quantum mechanics in terms of the electron density, which depends on only three spatial coordinates, rather than the many-body wavefunction, which depends on the coordinates of all electrons. This leads to a significant reduction in complexity and makes DFT one of the most efficient methods for electronic structure calculations in large systems.

:::{important}
The Hohenberg-Kohn theorems apply to non-degenerate ground states only! Here "ground state" means that the energy corresponding to this state is the lowest possible.
:::

Expliciting the terms of an $N$-electron Hamiltonian (see Eqs. [](#eq:many_body_H) and [](#eq:electronic_H)) one obtains

$$
\hat H_\text{el} = \sum_i \frac{p^2_i}{2 m_e} + \frac{1}{2} \sum_i \sum_{i \neq j} \frac{e^2}{|\vec r_i - \vec r_j|} + \sum_i V(\vec r_i),
$$

where the terms can be grouped to yield $\hat H_\text{el} = \hat F + \hat V$, where $F$ are the kinetic and electron-electron interaction energies, respectively, and $V$ is the "external" potential energy. The ground state $\Psi_0$ is uniquely determined by $N$ and $V(\{\vec r\})$, and we normalise it so that $\langle \Psi_0 | \Psi_0 \rangle = N$. Its associated electron density is

$$
n_0(\vec r) = \langle \Psi_0 | \rho(\vec r) | \Psi_0 \rangle = N \int |\Psi_0(\vec r, \vec r_2, \ldots, \vec r_N)|^2 d\vec r_2 \ldots d\vec r_N,
$$

where $\rho(\vec r) = \sum_i \delta(\vec r - \vec r_i)$ is the particle density operator. By using this definition, we can write the total potential energy as

\begin{align}
\langle \Psi_0 | \hat V | \Psi_0 \rangle &= \int \sum_i V(\vec r_i) |\Psi_0(\vec r_1, \vec r_2, \ldots, \vec r_N)|^2 d\vec r_1 d\vec r_2 \ldots d\vec r_N \\
& = N \int V(\vec r) |\Psi_0(\vec r, \vec r_2, \ldots, \vec r_N)|^2 d\vec r d\vec r_2 \ldots d\vec r_N \\
& = \int n_0(\vec r) V(\vec r) d\vec r.
\end{align}


:::{prf:theorem} The First Hohenberg-Kohn Theorem: Existence Theorem
For any system of interacting electrons in an external potential $V(\vec{r})$, the external potential is uniquely determined by the ground-state electron density $n_0(\vec{r})$, except for a trivial additive constant.
:::

:::{prf:proof}
The theorem is proven by contradiction. Assume that two different external potentials, $V_1(\vec{r})$ and $V_2(\vec{r})$, and therefore two different Hamiltonians, $\hat H_1 = \hat F + \hat V_1$ and $\hat H_2 = \hat F + \hat V_2$, lead to the same ground-state electron density $n_0(\vec{r})$. The corresponding ground-state wavefunctions, $\Psi_1$ and $\Psi_2$, as they correspond to different potentials. Their associated energies, which should also be different, are

\begin{align}
E_1 = \langle \Psi_1 | \hat H_1 | \Psi_1 \rangle\\
E_2 = \langle \Psi_2 | \hat H_2 | \Psi_2 \rangle,
\end{align}

Since both $\Psi_1$ and $\Psi_2$ are ground states, $E_1$ and $E_2$ are the lowest possible with the corresponding Hamiltonian. Therefore

\begin{align}
E_1 < \langle \Psi_2 | \hat H_1 | \Psi_2 \rangle & = \langle \Psi_2 | \hat H_2 | \Psi_2 \rangle + \langle \Psi_2 | (\hat H_1 - \hat H_2) | \Psi_2 \rangle =\\
&= E_2 + \int n_0(\vec r) [V_1(\vec r) - V_2(\vec r)] d\vec r\\
E_2 < \langle \Psi_1 | \hat H_2 | \Psi_1 \rangle & = \langle \Psi_1 | \hat H_1 | \Psi_1 \rangle + \langle \Psi_1 | (\hat H_2 - \hat H_1) | \Psi_1 \rangle =\\
&= E_1 + \int n_0(\vec r) [V_2(\vec r) - V_1(\vec r)] d\vec r.
\end{align}

Adding the two inequalities yeilds

$$
E_1 + E_2 < E_2 + E_1,
$$

which is clearly impossible. Therefore, the assumption must be false, and the external potential must be uniquely determined by the electron density (up to a constant).
:::

This theorem implies that the many-body wavefunction $\Psi(\vec{r}_1, \vec{r}_2, \dots, \vec{r}_N)$ is no longer the central object in quantum mechanics. Instead, the electron density $n(\vec{r})$ can be used to determine all ground-state properties. Specifically, the ground-state energy $E_0$ is a functional of the electron density, $E[n]$, and any other observable that depends on the external potential, such as forces on nuclei or response functions, is also determined by the electron density. Therefore, if one can find the correct ground-state electron density, one can in principle reconstruct the entire external potential and solve for all other properties of the system.

:::{prf:theorem} Second Hohenberg-Kohn Theorem: The Variational Principle
There exists a universal energy functional $E[n]$ of the electron density such that the correct ground-state density $n_0(\vec{r})$ minimizes this functional[^rayleigh-ritz], yielding the ground-state energy $E_0$.

This energy functional can be written as:

$$
E[n] = F[n] + \int V(\vec{r}) n(\vec{r}) \, d\vec{r},
$$ (eq:hohenberg_kohn_2)

where $F[n] = T[n] + U[n]$ includes the kinetic energy of the electrons and the electron-electron interaction energy, and $\int V(\vec{r}) n(\vec{r}) \, d\vec{r}$ is the external potential energy, representing the interaction with the external field (*e.g.* the electron-nucleus interaction potential in the Born-Oppenheimer approach). The Hohenberg-Kohn functional $F[n]$ is sometimes called "universal functional" since it does not depend on the specific system under study.
:::

:::{prf:proof}
The second Hohenberg-Kohn theorem is a straightforward application of the variational principle, whose proof I also report here for the sake of completeness. Let $\psi_n$ and $E_n$ be the eigenstates and eigenvalues of a Hamiltonian $\hat H$, with $E_0 < E_1 < E_2 < \ldots$, so that $\psi_0$ is the ground state, $\psi_1$ is the first excited state, *etc*. The eigenstates are orthonormal, *i.e.* $\langle \psi_i | \psi_j \rangle = \delta_{ij}$. Now consider a (normalised) generic wavefunction $\psi$, which can be expressed in the basis of the eigenstates as

$$
\psi = \sum_i c_i \psi_i,
$$

where $\sum_i |c_i|^2 = 1$, which implies that

$$
|c_0|^2 = 1 - \sum_{i > 0} |c_i|^2.
$$ (eq:variational_coeff)

The expectation value of the Hamiltonian on $\psi$ is

\begin{align}
\langle \psi | \hat H | \psi \rangle & = \langle \sum_i c_i \psi_i | \hat H | \sum_j c_j \psi_j \rangle = \sum_{i,j} c_i^* c_j \langle \psi_i | \hat H | \psi_j \rangle = \\
& = \sum_{i,j} c_i^* c_j E_j \langle \psi_i | \psi_j \rangle = \sum_j |c_j|^2 E_j = |c_0|^2 E_0 + \sum_{j > 0} |c_j|^2 E_j,
\end{align}

which, by applying Eq. [](#eq:variational_coeff), becomes

$$
\langle \psi | \hat H | \psi \rangle = E_0 + \sum_{j > 0} |c_j|^2 (E_j - E_0) > E_0,
$$

since $|c_j|^2 \geq 0$ and $E_j - E_0 > 0$ $\forall j > 0$.

Back to the Hohenberg-Kohn theorem, we know that the electron density $n(\vec r)$ fully determines the potential $V(\vec r)$ and the ground state $\Psi_0$. Consider a trial electron density $n'(\vec{r})$ that is not the true ground-state density $n_0(\vec{r})$. This will be associated to an external potential $V'$ and to a wave function $\Psi'$. By applying the variational principle, we find that the expectation value of the true Hamiltonian on this wave function is

$$
\langle \Psi' | \hat H_\text{el} | \Psi' \rangle = E[n'] \leq E[n_0] = \langle \Psi | \hat H_\text{el} | \Psi \rangle,
$$

where the equality holds if and only if $n'(\vec r) = n_0(\vec r)$. 
:::

The second Hohenberg-Kohn theorem provides a practical way to find the ground-state electron density by minimizing the energy functional $E[n]$. This is the basis of modern DFT numerical calculations:

1. Start with a trial electron density $n(\vec{r})$.
2. Calculate the total energy functional $E[n]$.
3. Adjust the electron density to minimize the energy.
4. The density that minimizes the energy is the true ground-state density, and the corresponding energy is the ground-state energy.

This approach bypasses the need to solve the many-body Schrödinger equation directly, instead focusing on finding the electron density that minimizes the energy functional. The challenge in DFT is to find an accurate expression for the universal functional $F[n] = T[n] + U[n]$, as this functional is not known explicitly for interacting electrons. As we will now see, various approximations, such as the Local Density Approximation (LDA) and Generalized Gradient Approximation (GGA), have been developed to approximate this functional.

[^rayleigh-ritz]: This is analogous to the [Rayleigh-Ritz variational principle for wavefunctions](https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Ritz_method#In_quantum_physics).

# The Kohn-Sham approximation

The Hohenberg-Kohn theorems laid the foundation for Density Functional Theory (DFT) by demonstrating that all ground-state properties of a many-electron system are determined by the electron density $n(\vec{r})$. However, the exact form of the universal energy functional $E[n]$, which includes the kinetic energy and electron-electron interaction energy, is unknown for interacting electrons. This is the key challenge in practical DFT calculations.

The Kohn-Sham (KS) approximation, introduced by [Walter Kohn and Lu Jeu Sham in 1965](doi:10.1103/PhysRev.140.A1133), provides a solution to this problem by mapping the interacting electron system onto an auxiliary system of non-interacting electrons that has the same ground-state electron density. This makes it possible to treat the complex interactions in a computationally efficient way while still capturing the essential physics of the many-body problem.

## The Kohn-Sham equations

The central idea of the Kohn-Sham approach is to introduce a fictitious system of non-interacting electrons that reproduces the exact ground-state electron density of the real, interacting system. For this non-interacting system, the kinetic energy can be computed exactly in terms of single-particle orbitals, which in this context are called the Kohn-Sham orbitals, $\psi_i$. The electron density is also constructed from the Kohn-Sham orbitals as:

$$
n(\vec{r}) = \sum_{i=1}^{N} |\psi_i(\vec{r})|^2.
$$

Instead of attempting to directly approximate the total energy functional $E[n]$ for the interacting system, we decompose it into parts that can be computed exactly and parts that require approximations. The total energy functional in the Kohn-Sham formalism is written as:

$$
E[n] = T_s[n] + \int V(\vec{r}) n(\vec{r}) \, d\vec{r} + U[n] + E_{\text{xc}}[n],
$$ (eq:kohn-sham_functional)

where $T_s[n]$ is the kinetic energy of the non-interacting electrons, $U[n]$ is the Hartree (Coulomb) energy, representing the classical electrostatic interaction between electrons, and $E_{\text{xc}}[n]$ is the exchange-correlation energy, which includes all the many-body effects not captured by the other terms (*e.g* the electron-electron effective repulsion due to the Pauli exclusion principle).

To find the ground-state density, the Kohn-Sham approach involves solving a set of single-particle equations, known as the Kohn-Sham equations. These equations describe the motion of non-interacting electrons in an effective potential $v_{\text{eff}}(\vec{r})$, and can be obtained by using the fact that the functional in Eq. [](#eq:kohn-sham_functional) is minimised by the ground-state electron density. Therefore, $\delta E / \delta n = 0$. If we carry out the functional derivative with the constraint that the Kohn-Sham orbitals are orthonormal, we obtain (see Appendix B of @giustino2014materials for the full derivation)

$$
\left( -\frac{\hbar^2}{2m} \nabla^2 + v_{\text{eff}}(\vec{r}) \right) \psi_i(\vec{r}) = \epsilon_i \psi_i(\vec{r}),
$$ (eq:kohn-sham)

$\epsilon_i$ is the single-particle energy corresponding to $\psi_i$, and the effective potential, which includes the effects of the external potential, the Hartree potential, and the exchange-correlation potential, is given by

$$
v_{\text{eff}}(\vec{r}) = V(\vec{r}) + \int \frac{n(\vec{r'})}{|\vec{r} - \vec{r'}|} d\vec{r'} + v_{\text{xc}}(\vec{r})
$$

The Kohn-Sham equations must be solved self-consistently:

1. Start with an initial guess for the electron density $n(\vec{r})$.
2. Construct the effective potential $v_{\text{eff}}(\vec{r})$.
3. Solve the Kohn-Sham equations to obtain the orbitals $\psi_i(\vec{r})$.
4. Calculate a new electron density from the orbitals.
5. Repeat until the electron density converges.

## The exchange-correlation functional

The exchange-correlation functional $E_{\text{xc}}[n]$ plays a crucial role in DFT, as it captures the many-body effects of electron exchange and correlation that are not accounted for in the other terms of the Kohn-Sham functional. This term is responsible for the accuracy of DFT calculations, and its exact form is unknown.

Several approximate forms for the exchange-correlation functional have been developed, each offering a different balance between accuracy and computational cost. Here I list the most common approximations.

### Local Density Approximation

In the Local Density Approximation (LDA), the exchange-correlation energy density at each point in space is assumed to depend only on the local electron density $n(\vec{r})$:

$$
E_{\text{xc}}^{\text{LDA}}[n] = \int \epsilon_{\text{xc}}(n(\vec{r})) n(\vec{r}) d\vec{r}.
$$

The LDA is based on the exchange-correlation energy of a homogeneous electron gas. While it can be surprisingly accurate for systems with slowly varying densities (*e.g.*, bulk solids), it fails, sometimes dramatically, in cases where the electron density varies rapidly, such as in molecules or surfaces.

### Generalized Gradient Approximation

The Generalized Gradient Approximation (GGA) improves upon the LDA by incorporating not only the local electron density but also its gradient $\nabla n(\vec{r})$:

$$
E_{\text{xc}}^{\text{GGA}}[n] = \int \epsilon_{\text{xc}}(n(\vec{r}), \nabla n(\vec{r})) n(\vec r) d\vec{r}.
$$

GGA functionals, such as the [Perdew-Burke-Ernzerhof](doi:10.1103/PhysRevLett.77.3865) functional, are more accurate for systems with significant density variations, such as molecules, surfaces, and low-dimensional materials. They are widely used in practical DFT calculations.

### Hybrid Functionals

Hybrid functionals go beyond the local and semi-local approximations by mixing a portion of exact exchange energy from Hartree-Fock theory with the exchange-correlation energy from DFT. The most well-known hybrid functional is [B3LYP](doi:10.1021/j100096a001), which is popular in quantum chemistry for molecular systems. Hybrid functionals generally provide higher accuracy but at a significantly increased computational cost.

## Implications and limitations

The Kohn-Sham approximation is a powerful and practical method for implementing Density Functional Theory. By mapping the interacting electron problem onto a non-interacting system with the same electron density, it makes DFT calculations computationally tractable. In fact, DFT with the Kohn-Sham approximation can be applied to systems with hundreds or even thousands of atoms, from molecules to solids, and can provide a good balance between computational cost and accuracy, particularly with carefully chosen exchange-correlation functionals. Unfortunately, the choice of the functional greatly impacts the accuracy of the DFT calculations, and no single functional is universally accurate. 

Moreover, for some systems (*e.g.*, strongly correlated materials), DFT may fail to provide reliable results. For instance, the most approximate functionals (*e.g.*, LDA and GGA) suffer from self-interaction errors, where an electron incorrectly interacts with itself. This can lead to inaccuracies, particularly in systems with localized states, such as transition metal oxides and open-shell molecules. Moreover, DFT is primarily a ground-state theory, and its application to excited states is limited. [Time-dependent DFT](https://en.wikipedia.org/wiki/Time-dependent_density_functional_theory) extends DFT to excited states, but it also inherits the limitations of the exchange-correlation functionals used.

# The Hellmann-Feynman Theorem

In the sections above I have sketched a way of solving the electronic problem and obtain the ground-state electron density. However, in *ab initio* simulations we have to also evolve the positions of the nuclei. How do we connect $n(r)$ to the forces acting on the nuclear degrees of freedom? The Hellmann-Feynman theorem provides a convenient and efficient way to compute these forces within the framework of quantum mechanics, without requiring the explicit calculation of the derivatives of the wavefunction.

The theorem states that when a system is in an eigenstate of its Hamiltonian, the forces acting on a nucleus can be calculated directly from the electron density and the external potential. This result greatly simplifies force calculations in Density Functional Theory (DFT) and other quantum mechanical methods, making it possible to efficiently perform simulations of atomic motion.

:::{prf:theorem} The Hellmann-Feynman Theorem

Consider a quantum system described by a Hamiltonian $\hat{H}(\lambda)$, which depends on a parameter $\lambda$ (*e.g.*, the position of a nucleus in a molecule or solid). If $\psi(\lambda)$ is an eigenstate of the Hamiltonian with eigenvalue $E(\lambda)$, then the derivative of the eigenvalue with respect to the parameter $\lambda$ is given by the expectation value of the derivative of the Hamiltonian:

$$
\frac{dE(\lambda)}{d\lambda} = \left\langle \psi(\lambda) \left| \frac{d\hat{H}(\lambda)}{d\lambda} \right| \psi(\lambda) \right\rangle.
$$
:::

:::{prf:proof}
:label: prf:hellman-feynman
Consider the time-independent Schrödinger equation for a system with a Hamiltonian $\hat{H}(\lambda)$ that depends on a parameter $\lambda$:

$$
\hat{H}(\lambda) \psi(\lambda) = E(\lambda) \psi(\lambda)
$$

Taking the derivative of both sides of this equation with respect to $\lambda$ we obtain[^no_lambda]

$$
\frac{d\hat{H}}{d\lambda} \psi + \hat{H} \frac{d\psi}{d\lambda} = \frac{dE}{d\lambda} \psi + E \frac{d\psi}{d\lambda}.
$$

If we multiply both sides by $\psi^*$ and integrate over all space we get

$$
\left\langle \psi \left| \frac{d\hat{H}}{d\lambda} \right| \psi \right\rangle + \left\langle \psi \left| \hat{H} \right| \frac{d\psi}{d\lambda} \right\rangle = \frac{dE}{d\lambda} + E \left\langle \psi \left| \frac{d\psi}{d\lambda} \right. \right\rangle.
$$

Since $\psi$ is an eigenstate of $\hat{H}$, $\hat{H}\psi = E\psi$, so that the second term can be simplified:

$$
\left\langle \psi \left| \frac{d\hat{H}}{d\lambda} \right| \psi \right\rangle + E \left\langle \psi \middle| \frac{d\psi}{d\lambda} \right\rangle = \frac{dE}{d\lambda} + E \left\langle \psi \left| \frac{d\psi}{d\lambda} \right. \right\rangle.
$$

The terms involving the derivative of the wavefunction cancel out, yielding

$$
\frac{dE(\lambda)}{d\lambda} = \left\langle \psi(\lambda) \left| \frac{d\hat{H}(\lambda)}{d\lambda} \right| \psi(\lambda) \right\rangle.
$$

[^no_lambda]: In the derivation I drop the explicit dependence on $\lambda$ for clarity.
:::

In the context of DFT, the Hellmann-Feynman theorem is particularly useful for calculating the forces on nuclei during molecular dynamics or geometry optimization. The force $\vec{F}_I$ on nucleus $I$ is given by the negative gradient of the total energy with respect to the position of that nucleus:

$$
\vec{F}_I = -\nabla_{\vec{R}_I} E_{\text{total}}(\{\vec{R}\})
$$

Consider the full Hamiltonian of a molecular system, Eq. [](#eq:full_H): the only two terms that depend on the positions of nuclei $\{\vec{R}_I\}$ are the nuclear-nuclear and nuclear-electron interactions. Using the Hellmann-Feynman theorem, these two contributions can be explicitly calculated:

$$
\vec{F}_I =  -\int n(\vec{r}) \nabla_{\vec{R}_I} v(\vec{r}; \{\vec{R}\}) \, d\vec{r} - \nabla_{\vec{R}_I} V_\text{nuc-nuc}(\{\vec{R}\}),
$$

where we used the expression of the electron-nucleus interaction in terms of the electron density $n(r)$ (see Eq. [](#eq:hohenberg_kohn_2)).
This result shows that the force on a nucleus depends only on the electron density and the derivative of the external potential (*e.g.*, the Coulomb potential from the other nuclei). Importantly, it does not require calculating the derivatives of the wavefunctions with respect to the nuclear positions, which is computationally demanding.

## Implications and limitations

The Hellmann-Feynman theorem is widely used in *ab initio* simulations for calculating forces on atoms, particularly in DFT and other quantum chemistry methods. The theorem simplifies force calculations by avoiding the need for wavefunction derivatives. This makes molecular dynamics and geometry optimization in DFT feasible for large systems. Interestingly, the forces calculated using the Hellmann-Feynman theorem are exact as long as the wavefunction is an exact eigenstate of the Hamiltonian. In practice, errors can arise if the wavefunction is not fully converged, but these can often be minimized with proper numerical techniques. Moreover, the Hellmann-Feynman theorem also applies to systems where the parameter $\lambda$ represents something other than nuclear positions, such as the strength of an external field or the variation of a coupling constant, making it a versatile tool.

However, note that in some cases, additional corrections beyond the Hellmann-Feynman theorem are needed. For example, when using approximate wavefunctions or basis sets, so-called [Pulay forces](doi:10.1080/00268976900100941) arise from using a non-complete basis functions that depend on the nuclear positions, contributing to the total force calculation in methods where the basis set is not fixed (*e.g.*, Gaussian-type orbitals).

This can be demonstrated by using the following derivation (heavily inspired by [this answer by Michael F. Herbst](https://mattermodeling.stackexchange.com/a/900/674)). Consider a non-complete basis $\{ f \}$, which we use to express our wavefunction:

$$
\Psi = \sum_{i} c_i f_i,
$$

where $c_i$ are constants chosen so that the resulting $\Psi$ is a good approximation of the ground state, and it is associated to an energy $E$. Here "good" means that each coefficient $c_i$ obeys the variational principle:

$$
0 = \frac{dE}{dc_i} = 2 \left\langle \Psi \middle| H \middle| \frac{d\Psi}{dc_i} \right\rangle = 2 \left\langle \Psi \middle| H \middle| f_i \right\rangle
$$ (eq:pulay-variational)

Since we are interested in forces, instead of taking the derivative of the energy with respect to a generic parameter $\lambda$, we derive the energy with respect to the nuclear positions:

\begin{equation}
\frac{dE}{d\vec{R}}=\left\langle\Psi \middle|\frac{dH}{d\vec{R}} \middle| \Psi\right\rangle
+ 2 \left\langle \Psi \middle| H \middle| \frac{d\Psi}{d\vec{R}} \right\rangle.
\end{equation}

As we have seen in [](#prf:hellman-feynman), the second ("Pulay") term vanishes if $\Psi$ is an eigenstate. Let's see what happens in this case. Consider as an example the derivative with respect to $R_1$. Its Pulay term is

$$
\left\langle \Psi \middle| H \middle| \left( \sum_i c_i \frac{df_i}{dR_1} \right)\right\rangle
$$

When does this term vanish? There are two possibilities:

1. The derivatives $df_i/dR_1$ are all zero, *i.e.* if the basis functions are independent of atomic positions, *e.g.* plane waves.
2. If the derivatives $\frac{df_1}{dR_1}$ and $\frac{df_2}{dR_1}$ are themselves basis functions or can be exactly represented by the basis. In this case,
$$\frac{df_i}{dR_1} = \sum_j k_{ij} f_i$$
for some values of $k_{ij}$. Now the Pulay term can be rewritten by leveraging Eq. [](#eq:pulay-variational) as
\begin{align}
\left\langle \Psi \middle| H \middle| \left( \sum_i c_i \frac{df_i}{dR_1} \right)\right\rangle & = \left\langle \Psi \middle| H \middle| \left( \sum_{i,j} c_i k_{ij} f_j \right)\right\rangle = \\
& = \sum_{i,j} c_i k_{ij} \langle \Psi | H | f_j \rangle = 0.
\end{align}

# The Car-Parrinello method

The Car-Parrinello method, introduced by [Roberto Car and Michele Parrinello in 1985](doi:10.1103/PhysRevLett.55.2471), revolutionized the field of *ab initio* molecular dynamics (AIMD) by combining classical molecular dynamics with electronic structure calculations based on DFT. The key innovation of this method is its ability to simultaneously evolve both nuclear and electronic degrees of freedom within a unified framework.

As discussed [earlier](#sec:born-oppenheimer), in traditional Born-Oppenheimer molecular dynamics, the forces on the nuclei are obtained by solving the electronic structure problem for each configuration of the nuclei, possibly with DFT-like approaches. The nuclei are then moved according to classical equations of motion, and this process is repeated iteratively. Although this approach is accurate, it is computationally expensive because a self-consistent field (SCF) calculation must be performed at every MD time step to obtain the forces.

In contrast, the Car-Parrinello molecular dynamics (CPMD) method evolves both the nuclear positions and the electronic wavefunctions simultaneously. This avoids the need for repeated SCF calculations, making the method more efficient while retaining accuracy. The key idea is to treat the Kohn-Sham orbitals as fictitious dynamical variables governed by a classical-like equation of motion, allowing them to follow the nuclear motion without needing to re-optimize the electronic structure at each step.

## The equations of motion

The Car-Parrinello method is based on a Lagrangian formulation that includes both the nuclear and electronic degrees of freedom. The total Lagrangian $L$ for the system is given by:

$$
L = \frac{1}{2} \sum_I M_I \dot{\vec{R}}_I^2 + \frac{\mu}{2} \sum_i \langle \dot{\psi}_i | \dot{\psi}_i \rangle - E_\text{tot}[\{\psi\}, \{\vec{R}\}] + \sum_{ij} \lambda_{ij} (\langle \psi_i | \psi_j \rangle - \delta_{ij}),
$$ (eq:cp_lagrangian)

where $\vec{R}_I$ are the nuclear positions, $M_I$ are the nuclear masses, $\dot{\vec{R}}_I$ are the nuclear velocities, $\psi_i$ are the Kohn-Sham orbitals, $\mu$ is a fictitious mass-like parameter associated with the electronic degrees of freedom, $\dot{\psi}_i$ are the time derivatives ("velocities") of the Kohn-Sham orbitals, and $E_\text{tot}[\{\psi\}, \{\vec{R}\}] = E[\{\psi\}, \{\vec{R}\}] + E_\text{nuc,nuc}(\{\vec{R}\})$ is the sum of the DFT (Kohn-Sham) energy functional and of the nuclear-nuclear interaction potential.

The first term in Eq. [](#eq:cp_lagrangian) represents the kinetic energy of the nuclei, the second term represents the kinetic energy of the fictitious electronic degrees of freedom, the third term $E_\text{tot}[\{\psi_i\}, \{\vec{R}_I\}]$ is the total potential energy, which includes contributions from the nuclear-nuclear, nuclear-electronic, and electron-electron interactions, while the fourth term enforces the orthonormality of the Kohn-Sham orbitals.

From this Lagrangian, one can derive the equations of motion for both the nuclei and the electronic wavefunctions using the Euler-Lagrange formalism:

\begin{align}
M_I \ddot{\vec{R}}_I & = -\frac{\partial E_\text{tot}[\{\psi_i\}, \{\vec{R}_I\}]}{\partial \vec{R}_I} + \sum_{i,j} \lambda_{ij} \frac{\partial}{\partial \vec R_I} \langle \psi_i | \psi_j \rangle\\
\mu \ddot{\psi}_i(\vec{r}) & = -\frac{\delta E_\text{tot}[\{\psi_i\}, \{\vec{R}_I\}]}{\delta \psi_i^*(\vec{r})} + \sum_{j} \lambda_{ij} \psi_j
\end{align}

Here, the fictitious mass parameter $\mu$ controls the dynamics of the electronic wavefunctions. These equations ensure that the electronic wavefunctions follow the nuclear motion, remaining close to the instantaneous ground state of the electronic structure, without the need for explicit SCF iterations at each time step. Note that in the functional derivative that appears in the electronic equations of motions, only the terms of the non-interacting (Kohn-Sham) Hamiltonian survives, since the nuclear-nuclear interaction does not depend on the $\psi_i$ orbitals.

The constant of motion (*i.e.* the quantity that is conserved) corresponding to the time-independence of the Car-Parrinello Langrangian is

$$
E = \frac{1}{2} \sum_I M_I \dot{\vec{R}}_I^2 + \frac{\mu}{2} \sum_i \langle \dot{\psi}_i | \dot{\psi}_i \rangle + E_\text{tot}[\{\psi\}, \{\vec{R}\}] \equiv E_\text{phys} + T_e.
$$

As well put in an excellent review,

> As the kinetic energy $T_e$ of the electronic degrees of freedom will vary with time, so will $E_\text{phys}$, the total energy of the relevant Born-Oppenheimer system. For a physically meaningful simulation, it is therefore necessary that $T_e$ only performs bound oscillations around a constant small value. The interpretation of this state of the system is that the total CP dynamical system consists of two adiabatically decoupled subsystems, the cold electronic degrees of freedom and the nuclear degrees of freedom at the relevant physical temperature. The adiabatic separation will never be complete and it will be necessary to choose the mass parameter $\mu$ in a way that the two systems stay decoupled over the full range of the simulation. However, the choice of $\mu$ also directly affects the efficiency of the simulation, as the maximal time step of the numerical integration is proportional to $\sqrt{\mu}$. The choice of an optimal value of the free parameter $\mu$ in a CP simulation needs care and has been discussed at length [in the literature].
>
> -- [](doi:10.1002/wcms.90)

Indeed, the success of the Car-Parrinello method relies on a careful choice of the fictitious mass parameter $\mu$. This parameter determines the timescale on which the electronic wavefunctions evolve. To ensure that the electrons remain adiabatically "slaved" to the nuclei, $\mu$ must be small enough that the electronic motion remains much faster than the nuclear motion. However, $\mu$ cannot be too small, or the time step required for stable integration of the electronic equations of motion would become prohibitively small. Therefore, $\mu$ should be chosen such that the electronic degrees of freedom evolve rapidly enough to stay close to the instantaneous ground state, but not so fast that the computational efficiency is compromised (which is easier said than done).

## Implications and limitations 

The Car-Parrinello method offers several advantages over traditional approaches (such as those based on the Born-Oppenheimer approximation). By evolving the electronic wavefunctions dynamically, CPMD avoids the need for repeated SCF calculations, leading to significant savings in computational cost, particularly for large systems. Moreover, the method retains the accuracy of DFT-based molecular dynamics since the electronic wavefunctions are kept close to the ground state throughout the simulation. Finally, the treatment of both nuclear and electronic degrees of freedom within a single Lagrangian formalism provides a natural and elegant way to simulate systems where electronic structure and nuclear motion are strongly coupled.

However, the introduction of the fictitious kinetic energy for the electronic wavefunctions means that energy is not strictly conserved over long simulation times. Although this can be controlled by adjusting the fictitious mass $\mu$, it remains a source of error that requires careful monitoring.
 The need to integrate both the nuclear and electronic equations of motion imposes restrictions on the time step size, which can be smaller than in traditional Born-Oppenheimer molecular dynamics, depending on the choice of $\mu$. The Car-Parrinello method assumes that the electronic wavefunctions remain close to the ground state during nuclear motion. However, in systems where non-adiabatic effects (such as electronic excitations) are important, more advanced methods may be required. Extensions of the Car-Parrinello method have been developed to address some of these limitations. For example, extended Lagrangian methods introduce additional degrees of freedom to improve energy conservation, while [time-dependent DFT](https://en.wikipedia.org/wiki/Time-dependent_density_functional_theory) can be used to model systems where excited-state dynamics play a significant role.

# Running simulations

While I will not provide any example, here I mention some of the most used open-source software packages for performing *ab initio* molecular dynamics simulations:

* [CP2K](https://www.cp2k.org/)
* [SIESTA](https://departments.icmab.es/leem/siesta/)
* [Quantum ESPRESSO](https://www.quantum-espresso.org/)

Most of these packages can be used to simulate systems at varying (and even mixed) levels of details.

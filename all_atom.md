---
title: Classical molecular dynamics and all-atom simulations
exports:
   - format: pdf
---

```{tip}
The main references for this part are @frenkel2023understanding, @vsulcintroduction (which can be downloaded [here](https://www.public.asu.edu/~psulc/myimages/chapter.pdf)), @schlick2010molecular and @leach2001molecular.
```

Quantum mechanics provides a rigorous framework for describing the behavior of molecules that explicitly includes the effect of the electrons. However, the calculations are very heavy and thus the evolution of a system can be followed for short times only. In addition, the computational complexity of simulating quantum systems often scales super-linearly with the size of the problem[^scaling]. This poses significant challenges when attempting to model large-scale many-body systems accurately and for long times.

One of the fundamental distinctions between quantum and classical approaches lies in the treatment of electron contributions. In quantum calculations, the interactions and motions of electrons are computed explicitly, for instance with DFT methods, as discussed in the previous Chapter. This level of detail enables precise predictions of electronic structure and properties but comes at a high computational cost, particularly as the number of electrons and/or atoms increases.

In contrast, classical force fields adopt a simplified approach where the contributions of electrons are averaged out. Instead of explicitly modeling individual electron behaviors, classical force fields approximate the interactions between atoms using simplified mathematical models based on classical mechanics. Note that there exist also "mixed" methods, where some degrees of freedom are accounted for by using a quantum-mechanical treatment, while others are treated classically, *e.g.* the nuclear degrees of freedom in the CPMD approach.

In the rest of the Chapter (unless stated otherwise), we consider systems composed of $N$ interacting point particles following Newton's equations of motion set by a classical Hamiltonian:

$$
H = \sum_{i=1}^N \frac{p_i^2}{2m_i} + V(\{ \vec r_i \}),
$$

where $\vec p_i$ and $\vec r_i$ are the momentum and coordinates of particle $i$, $m_i$ is its mass and $V(\{ \vec r_i \})$ is the total potential energy. 

:::{seealso} Python implementation
Some of the algorithms presented in this section have been implemented in this [Jupyter notebook](./notebooks/MD.ipynb).
:::

[^scaling]: The complexity ranges from $\mathcal{O}(e^N)$ for brute-force implementations, to $\mathcal{N^3}$ for many DFT codes, but can be linear in some cases (see *e.g.* [](doi:10.1088/0034-4885/75/3/036503)).

# Molecular dynamics

In a molecular dynamics simulation, the atoms, or particles, follow Newton's equations of motions. As a result, the dynamics happens on a hypersurface defined by $E = \text{const}$, where $E$ is the energy of the system. In practice, the equations of motions are solved iteratively by discretising time: the quantities of interest (position $r$, velocity $v$ and force $F$) at time $t$ are used to obtain those at time $t + \Delta t$, where $\Delta t$ is the integration *time step*. The flow of a simple MD program is:

1. The simulation parameters (*e.g.* temperature, density, time step) are read and initialised.
2. The initial configuration (*i.e.* the initial positions and velocities) is read or generated.
3. The simulation runs, iterating the following steps:
   1. The forces (and possibly torques) on all particles are computed.
   2. The positions and velocities are updated according to Newton's equations.
4. The simulation stops when some condition is met (*e.g.* number of steps run).

The following snippet shows a pseudo-code implementation of the foregoing algorithm:

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

:::{important}
Nowadays there are sophisticated codes that can be used, in principle, to run MD simulations with no knowledge any of the algorithms employed. This is **very** dangerous, as it is very easy to produce wrongs results that appear correct if one does not know what they are doing. Therefore, I urge you to try to understand how (and why) these algorithms work.
:::

## Initialisation

Initialising a simulation is a trivial task when simulating simple systems, *e.g.* gases, most liquids and even many solids, but can become extremely difficult for highly heterogeneous systems interacting through complicated, and possibly steep, interactions. Biomacromolecules modelled at the all-atom level are definitely within this class of systems.

In general, the initial configuration should be such that there is no sensible overlap between the system's constituents (atoms, molecules, particles, *etc.*), in order to avoid large initial forces that could lead to numerical instabilities that could lead to crashes or, which is worse, unphysical behaviour, such as two chains crossing each other. For simple systems this can be done in several ways:

* Particle coordinates are chosen randomly one after the other, ensuring that the distance between the new particle and any other particle is larger than sum threhsold. This works great as long as the density is not too high.
* Particles are placed on a lattice, making it possible to go to generate highly dense systems.

For more complicated systems, the generation of the initial configuration is often done in steps. In some cases, and I will show some examples, it is also common to run short-ish simulations where the dynamics is constrained (*e.g.* some degrees of freedom such as bond lengths or angles are kept fixed) and the energy is minimised to remove (part of) the initial stress.

In the case of highly heterogeneous systems, it is very common to use external tools to build part of (or the whole) system. Take for instance the simulation of protein-membrane systems (*e.g.* [](#fig:GPCR)), where one has to build a membrane with a given composition and density, embed one or more proteins, add ions and solvate the system.

What about velocities (or, in general, momenta)? Here initialisation is more straightforward, but there are some subtleties that can lead to errors later on if one is not careful. First of all, recall that, in a system in thermal equilibrium,

$$
\langle v_{i,\alpha}^2 \rangle = \frac{k_B T}{m},
$$ (eq:equipartition)

where $v_{i,\alpha}$ is the $\alpha$-th component of the velocity of particle $i$. We can then define the instantaneous temperature of a system of $N$ particles[^only_velocities] as

$$
k_B T(t) = \sum_{i, \alpha}^N \frac{m v_{i,\alpha}^2}{d N},
$$ (eq:instantaneous_T)

where $d$ is the dimensionality of the system. We now use this relation in the following initialisation procedure:

1. Randomly extract each velocity component of each particle from a distribution with zero mean. A good choice is a Gaussian, but a uniform distribution will also do (but don't be lazy!).
2. Set the total momentum of the system to zero by computing $\langle v_\alpha \rangle = \frac{1}{N} \sum_i^N v_{i, \alpha}$ and subtracting it from each $v_{i, \alpha}$.
3. Rescale all velocities so that the initial temperature $T(0)$ matches the desired value $T$, *i.e.* multiply each component by $\sqrt{T/T(0)}$.

If you chose to use a Gaussian distribution for step 1., then you will end up with velocities distributed according to the Maxwell-Boltzmann probability distribution at temperature $T$, *viz.*

$$
P(v) = \left(\frac{\beta}{2 \pi m}\right)^{3/2} \exp\left(-\frac{\beta m v^2}{2}\right).
$$ (eq:maxwell-boltzmann)

However, since the initial configuration is most likely *not* an equilibrium configuration, the subsequent equilibration will change the average temperature, so that $\langle T(t) \rangle \neq T$. We will see later on how to couple the system to a [thermostat](#sec:thermostats) to enforce thermal equilibrium at the desired temperature.

[^only_velocities]: here I assume that we are dealing with point-like objects.

## Reduced units

In molecular simulations, reduced units are employed to simplify calculations by scaling physical quantities relative to characteristic properties of the system, such as particle size or interaction strength. This approach not only makes the equations dimensionless and more general, but it also prevents the numerical instability that can arise from dealing with values that are extremely small or large, which can occur when using standard physical units (*e.g.* SI units). For example, distances are often scaled by the size of the particles, *e.g.* the particle diameter $\sigma$, energies by a characteristic interaction energy, *e.g.* the well depth of the attraction $\epsilon$, and masses by the particle mass, $m$. Any other unit is then expressed as a combination of these basic quantities; for instance, given $\sigma$, $\epsilon$, and $m$, time is in units of $\sigma \sqrt{m / \epsilon}$, and pressure is in units of $\epsilon / \sigma^3$. This allows quantities like temperature to be expressed as $T^* = \frac{k_B T}{\epsilon}$, reducing the need to handle numbers with extreme magnitude, which can hinder numerical stability and therefore lead to inaccuracies in the simulation. As an alternative, it is also possible to use the characteristic constants directly as units of measurements. For instance, a distance may be written as $L = 10 \, \sigma$, or an energy as $U = 0.1 \, \epsilon$.

## The force calculation

This is the most time-consuming part of an MD simulation. It consists of evaluating the force (and, for rigid bodies, the torque) acting on each individual particle. Since, in principle, every particle can interact with every other particle, for a system of $N$ objects the number of interacting contributions will be equal to the number of pairs, $N (N - 1) / 2$. Therefore, the time required to evaluate all forces scales as $\mathcal{O}(N^2)$. Although, as we will see later, there are methods that can bring this complexity down to $\mathcal{O}(N)$, here I present the simplest case.

Consider, for simplicity, a system made of $N$ like atoms interacting through the Lennard-Jones potential modelling [Van der Walls forces](#sec:van-der-waals). The atom-atom interaction depends only on the inter-atomic distance $r$, and reads

$$
V_\text{LJ}(r) = 4 \epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right),
$$

where $\epsilon$ is the potential well depth and $\sigma$ is the atom diameter. The force associated to this potential is then

$$
\vec{f}(\vec r) = - \vec \nabla V_\text{LJ}.
$$

(sec:cut-off)=
### Interaction cut-off

:::{important} The range of the interaction
A "short-range potential" is a definition that applies to all those potentials that can be truncated at a certain distance $r_c$ (called *cut-off distance*), and whose tail contribution (*i.e.* the contribution for $r > r_c$) is finite. This can happen either because $V(r > r_c) = 0$, or if the potential decays to zero quickly enough. Indeed, the tail correction can be estimated by considering the continuous limit, *viz.*

$$
U_\text{tail} \approx \frac{N\rho}{2} \int_{r_c}^\infty V(r) d \vec{r},
$$ (eq:U_tail)

where $\rho$ is the number density of the system, and $d\vec{r} = 4 \pi r^2 dr$ in 3D and $d\vec{r} = 2 \pi r dr$ in 2D. In order to yield a finite value, the integrand in Eq. [](#eq:U_tail) should decay faster than $r^{-1}$. If $V(r)$ decays at large distances as $\sim r^{-\alpha}$, the convergence condition becomes $\alpha > d$, where $d$ is the dimensionality.
:::

Given the definition above, the LJ potential is short ranged[^long-range_interactions]. As a result, the total potential energy of a particle is dominated by the contributions of those particles that are closer than some cut-off distance $r_c$. Therefore, in order to save some computing time, it is common to truncate the interaction at $r_c$, so that the potential reads

(eq:truncation)=
\begin{align}
V_\text{tr}(r) = \begin{cases}
4 \epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right) & \; \text{if} \; r \leq r_c\\
0 & \; \text{otherwise.}
\end{cases}
\end{align}

In this way, only pairs of particles closer than $r_c$ feel a mutual interaction. However, the potential has a discontinuity, and therefore the force diverges, at $r = r_c$. The discontinuity introduces errors in the estimates of several quantities. For instance, the potential energy will lack the contributions of distance atoms, which can be estimated through Eq. [](#eq:U_tail) and for a 3D LJ potential is

$$
U_\text{tail} \approx 2 \pi N\rho \int_{r_c}^\infty V(r) r^2 dr = \frac{8}{3} \pi N \rho \epsilon \sigma^3 \left[ \frac{1}{3}\left(\frac{\sigma}{r_c} \right)^9 - \left(\frac{\sigma}{r_c} \right)^3 \right].
$$ (eq:U_tail_LJ)

Another important observable affected by the truncation is the pressure. The correction can be estimated by using the virial theorem and Eq. [](#eq:U_tail_LJ) as

$$
\Delta P_\text{tail} = 2 \pi \rho^2 \int_{r_c}^\infty \vec{r} \cdot \vec{f}(r) r^2 dr = \frac{16}{3} \pi \rho^2 \epsilon \sigma^3 \left[ \frac{2}{3}\left(\frac{\sigma}{r_c} \right)^9 - \left(\frac{\sigma}{r_c} \right)^3 \right].
$$

Note that this is the contribution that should be added to $P$ if one wishes to estimate the pressure of the true (untruncated) potential, rather than the true pressure of the truncated potential[^true_pressure]. This is a specific example of a more general lesson: one should always be aware of the consequences of any approximation that is introduced in the model, weighing in advantages and disadvantages and, if possible, finding a way to correct the results *a posteriori*, as in this case.

A common way of removing the divergence in Eq. [](#eq:truncation) is to truncate and shift the potential:

\begin{align}
V_\text{tr,sh}(r) = \begin{cases}
4 \epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right) - V(r_c) & \; \text{if} \; r \leq r_c\\
0 & \; \text{otherwise.}
\end{cases}
\end{align}

This is especially important in constant-energy MD simulations (*i.e.* simulations in the NVE ensemble), since the discontinuity of the truncated (but not shifted) potential would greatly deteriorate the energy conservation. In this case the pressure tail correction remains the same, while the energy requires an additional correction on top of Eq. [](#eq:U_tail_LJ), which accounts for the average number of particles that are closer than $r_c$ from a given particle, multiplied by $\frac{1}{2} V(r_c)$.

:::{warning} TODO
Add part about smoothing
:::

[^long-range_interactions]: But Coulomb and dipolar interactions are not, and has to be treated differently, as we will see later on.
[^true_pressure]: Check @frenkel2023understanding if you want an expression for the latter.

### Minimum image convention and periodic boundary conditions

The number of particles in a modern-day simulation ranges from hundreds to millions, thus being very far from the thermodynamic limit. Therefore, it is not surprising that finite-size effects are always present (at least to some extent). In particular, the smaller the system, the larger boundary effects are: since in 3D the volume scales as $N^3$ and the surface as $N^2$, the fraction of particles that are on the surface scales as $N^{-1/3}$, which is a rather slowly decreasing function of $N$. For instance, if $N = 1000$, in a cubic box of volume $V = L^3$ more than half of the particles are on the surface. One would have to simulate $\approx 10^6$ particles to see this fraction decrease below $10\%$!

Those pesky boundary effects can be decreased by using periodic boundary conditions (PBCs): we get rid of the surface by considering the volume containing the system, which is often but not always a cubic box, one cell of an infinite periodic lattice made of identical cells. Then, each particle $i$ iteracts with any other particle: not only with those in the original cell, but also with their "images", including its own images, contained in all the other cells. For instance, the total energy of a (pairwise interacting) cubic system of side length $L$ simulated with PBCs would be

$$
U_\text{tot} = \frac{1}{2} \sum_{i,j,\vec{n}}\phantom{}^{'} V(|\vec{r}_{ij} + \vec{n}L|),
$$ (eq:PBC_sum)

where $\vec r_{ij}$ is the distance between particles $i$ and $j$, $\vec{n}$ is a vector of three integer numbers $\in [-\infty, +\infty]$, and the prime over the sum indicates that the $i = j$ term should be excluded if $\vec n = (0, 0, 0)$. How do we handle such an infinite sum in a simulation? Keeping the focus on short-range interactions, it is clear that all terms with $|\vec{r}_{ij} + \vec{n}L| > r_c$ vanish, leaving only a finite number of non-zero interactions. In practice, Eq. [](#eq:PBC_sum) is evaluated by making sure that $r_c < L / 2$, so that each particle can interact **at most** with a single periodic image of another particle. Then, given any two particles $i$ and $j$, the vector distance $\vec r_{ij} = (x_{ij}, y_{ij}, z_{ij})$ between the closest pair of images is

(eq:minimum_image)=
\begin{align}
x_{ij} & = x_j - x_i - \text{round}\left(\frac{x_j - x_i}{L_x}\right) L_x\\
y_{ij} & = y_j - y_i - \text{round}\left(\frac{y_j - y_i}{L_y}\right) L_y\\
z_{ij} & = z_j - z_i - \text{round}\left(\frac{z_j - z_i}{L_z}\right) L_z,
\end{align}

where, for the sake of completeness, I'm considering a non-cubic box of side lengths $L_x$, $L_y$ and $L_z$, and $\text{round}(\cdot)$ is the function that rounds its argument to its closest integer. [](#fig:PBC) shows a 2D schematic of periodic-boundary conditions and of the minimum-image construction.

```{figure} figures/PBC.png
:name: fig:PBC
:align: center
:width: 500px

A schematic representation of periodic-boundary conditions and the minimum-image construction for a two-dimensional system. The system shown here contains 4 particles (labelled $i$, $j$, $k$, and $l$), with the black lines connecting particle $i$ with the closest periodic images of the other particles. Credits to [Abusaleh AA via Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Minimum_Image_Convention.png).
```

The total force acting on a particle $i$ can now be computed by summing up the contributions due to each particle $j \neq i$, where the distance is given by Eq. [](#eq:minimum_image).

The use of periodic boundary conditions gets rid of surface effects altogether, and more generally can greatly limit finite-size effects, even for small numbers of particles. However there are some important things that should always be kept in mind when using PBCs (which is by far the most common form of boundary conditions employed in molecular simulations):

1. If the system under study exhibit correlations that are of the same order of or exceed the box size, the system will "feel itself" through the PBCs, leading to spurious non-linear effects.
2. Only wavelengths compatible with the box size are allowed. This is true for any observable, and applies to a wide range of phenomena. For instance
   * If one wishes to simulate equilibrium crystals, the shape of the box and the length of its sides should be compatible with the lattice's unit cell.
   * Fluctuations with wavelengths longer than the box size are suppressed. This makes it hard to study critical phenomena, which have diverging correlation lengths.
3. Angular momentum is not conserved, since the system is **not** rotationally invariant (even if the Hamiltonian is).
4. When particles cross the box boundary to, for instance, exit the central cell, one of its images will enter it. It is important to keep track of the "original" particle to correctly compute the value of some observables (*e.g.* the mean-squared displacements, see below, or a chain connectivity).

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

This is the Verlet algorithm, which allows us to calculate the position $x(t + \Delta t)$ at the next time step using the current position $x(t)$, the previous position $ x(t - \Delta t)$, and the current acceleration $a(t)$. Although this basic Verlet method does not explicitly involve velocity, if we consider the expansion up to the second order it can be written as

$$
v(t) = \frac{x(t + \Delta t) - x(t - \Delta t)}{2\Delta t} + \mathcal{O}(\Delta t^2).
$$

Note the lower order of accuracy with respect to $x(t)$. We can modify the basic Verlet approach to include an explicit update for the velocity, leading to the Velocity Verlet method. Instead of relying on the positions from the previous and current time steps, the Velocity Verlet algorithm updates the position and velocity in a two-step process.

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
2. Update the position, $x(t + \Delta t) = x(t) + v(t + \Delta t / 2)\Delta t = x(t) + v(t) \Delta t + \frac{1}{2} a(t) \Delta t^2$ (*e.g.* Eq. [](#eq:velocity_verlet_x)).
3. Calculate the force (and therefore the acceleration) using the new position, $x(t + \Delta t) \to a(t + \Delta t) = F(t + \Delta t) / m$.
4. Update the velocity, second step, $v(t + \Delta t) = v(t + \Delta t / 2) + \frac{1}{2} a(t + \Delta t)\Delta t = v(t) + \frac{1}{2} \left[ a(t) + a(t + \Delta t) \right] \Delta t$ (*e.g.* Eq. [](#eq:velocity_verlet_v)).

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

## Energy conservation

The Velocity Verlet algorithm is [*symplectic*](https://en.wikipedia.org/wiki/Symplectic_integrator), which means that it conserves the energy (*i.e.* the value of the Hamiltonian)[^symplectic]. This is not a common property, as most schemes (such as the explicit Euler scheme or the famous Runge-Kutta method) are not symplectic, and therefore generate trajectories that do not live on the $H = \text{const}$ hypersurface. Since this is where the dynamics of systems evolving through Newton's equations move, MD integrators should always be symplectic. If this is the case, the resulting molecular dynamics simulations will conserve the total energy, which is a constant of motion. In turn, this means that MD simulations sample the microcanonical ensemble (as long as the system is ergodic, so that time and ensemble averages are equivalent).

Since we are discretising the equations, there are two caveats associated to the energy conservation:

1. If the time step $\Delta t$ is too large, then a deviation (drift) from the "correct" value of the energy will be observed. A good rule of thumb to avoid an energy drift, which would invalidate the simulation, is to calculate the highest vibrational frequency associated to the interaction potential employed, $\omega_c \equiv \sqrt{k_c / m}$, where $k_c$ is the highest curvature of the potential, *i.e.* the largest value of $d^2V(r)/dr^2$, and $m$ is the mass of the particle. Then, a sensible value for the integration time step is an order of magnitude less than the characteristic time associated to $\omega_c$, *i.e.* $\Delta t \approx \frac{1}{10} \frac{2 \pi }{\omega_c}$.
2. If the value of $\Delta t$ is appropriate, then the energy will be conserved *on average*, meaning that it will fluctuate around its average value with an amplitude, $\delta U \equiv \sqrt{\langle \Delta U^2 \rangle}$. If the integrator is of the Verlet family (*i.e.* Velocity Verlet), then $\delta U \sim \Delta t^2$. This behaviour can be exploited to ensure that codes and simulations are working correctly, at least from the point of view of the energy conservation. See [](#fig:energy_conservation) for an example.

```{figure} figures/energy_conservation.png
:name: fig:energy_conservation
:align: center
:width: 500px

The extent of the fluctuations of the total energy, $\delta U$, for a Lennard-Jones system simulated at $k_B T / \epsilon = 1.5$ and $\rho \sigma^3 = 0.36$ as a function of the time step $\Delta t$. The line is a quadratic fit. **Nota Bene:** this scaling requires that truncation errors are small, which is why I had to use a rather large cut-off, $r_c = 3.5$.
```

## Some observables

Since, in principle, we have access to the whole phase space, in MD simulations we can compute averages of any (classical) observable. Here I present a short list of some useful (and common) quantities.

(sec:compute_pressure)=
### Pressure

:::{tip}
The main source for this part is @allen2017computer.
:::

In the canonical ensemble, for any observable $A$ and generalised coordinate or momentum $h_k$, integrating by parts (and assuming reasonable boundary conditions) one finds

\begin{align}
\left\langle \frac{\partial A}{\partial h_k} \right\rangle &= \frac{1}{Q} \int \frac{\partial A}{\partial h_k} \exp(-\beta H) d\lbrace h_k \rbrace = \frac{1}{Q} \int \beta A \frac{\partial H}{\partial h_k} \exp(-\beta H) d\lbrace h_k \rbrace\\
&= \beta \left\langle A \frac{\partial H}{\partial h_k} \right\rangle,
\end{align}

which translates to the following generalised equipartition relation:

$$
\left\langle h_k \frac{\partial H}{\partial h_k} \right\rangle = k_B T.
$$ (eq:generalised_equipartition)

If we plug in a generalised momentum and sum over all momenta, Eq. [](#eq:generalised_equipartition) yields the usual equipartition principle:

$$
\left\langle \sum_{i=1}^N \frac{|\vec{p}_i|^2}{m_i} \right\rangle = 3 N k_B T,
$$

where $m_i$ is the mass of the $i$-th particle. By contrast, if we choose to use Cartesian coordinates as generalised coordinates in Eq. [](#eq:generalised_equipartition) and recall that the derivative of the Hamiltonian with respect to a particle coordinate is minus the total force acting on the particle along that coordinate, we find

$$
\left\langle \sum_{i=1}^N \vec{r}_i \cdot \vec{F}_i^{\rm tot} \right\rangle = -3 N k_B T,
$$

where $\vec{r}_i$ is the position of particle $i$. Note that here $\vec{F}_i^{\rm tot}$ is the sum of inter-molecular interactions, $\vec{F}_i^{\rm int}$, and external forces, $\vec{F}_i^{\rm ext}$. If the latter consist only of the forces exerted by the container walls (which keeps the system in a volume $V$ under a pressure $P$), we have

$$
\frac{1}{3} \left\langle \sum_{i=1}^N \vec{r}_i \cdot \vec{F}_i^{\rm ext} \right\rangle = -PV.
$$

Since $\vec{F}_i^{\rm tot} = \vec{F}_i^{\rm int} + \vec{F}_i^{\rm ext}$, we find

$$
P = \frac{N k_B T}{V} + \frac{1}{3V} \left\langle \sum_{i=1}^N \vec{F}_i^{\rm int} \cdot \vec{r}_i \right\rangle \equiv \langle P_{\rm inst} \rangle,
$$ (eq:pressure)

where I have defined the instantaneous pressure $P_{\rm inst}$. The second term in Eq. [](#eq:pressure) represents the contribution from interparticle forces, where $\vec{F}_i^{\rm int}$ is the force on particle $i$ due to all other particles. Applying Eq. [](#eq:pressure) makes it possible to evaluate the pressure in computer simulations.

### Radial distribution function

There are many ways to obtain structural information from simulation output. The simplest, and most straightforward, piece of such information is encoded in the pair correlation function, also known as the radial distribution function, $g(r)$.

Let us start by defining a function $n(r)$ such that $n(r)dr$ represents the number of particles found between $r$ and $r+dr$, assuming the origin is placed at a random particle. With this definition, the total number of particles contained within a sphere of radius $R$ centered on a random particle is:

$$
N(R) = \int_0^R n(r)dr.
$$

For an ideal gas, $n(r) = 4\pi r^2 \rho$, where $\rho$ is the number density $(N-1)/V$.

The function $g(r)$ is defined as:

$$
g(r) \equiv \frac{n(r)}{4\pi r^2 \rho},
$$

and it indicates how much denser the system is at a distance $r$ compared to an ideal gas. The function $g(r)$ is always positive and, for large values of $r$, $g(r) \to 1$. The $g(r)$ is a central observable in molecular simulations, and can be used to estimate many other observables (see *e.g.* @hansen2013theory).

```{figure} figures/rdf.png
:name: fig:rdf
:align: center
:width: 550px

The radial distribution function of a Lennard-Jones fluid at $k_B T / \epsilon = 0.71$ and $\rho \sigma^3 = 0.844$. Superimposed on the plot is a cartoon that shows the interpretation of the features of the $g(r)$, highlighting the first and second shells around the central particle (coloured in red). Adapted from [Wikipedia](https://en.wikipedia.org/wiki/Radial_distribution_function).
```

[](#fig:rdf) shows a typical shape of $g(r)$, together with a 2D cartoon that gives an idea of how the structural features of the system affects the resulting radial distribution function. Near the origin, $g(r) = 0$ due to the excluded volume. This is always observed in simple liquids and, more generally, in systems with short-range repulsions (such as atom-based systems). For soft (or ultrasoft) potentials, which are typically effective, $g(r)$ may not be zero even at the origin.

[](#fig:rdf) also conveys a second message: the packing of objects with excluded volume inevitably generates oscillations in the $g(r)$ function, with a periodicity determined by the diameter of the particles themselves.

```{figure} figures/rdf_water.png
:name: fig:rdf_water
:align: center
:width: 400px

The (a) $O-O$ and (b) $O-H$ radial distribution functions of SPC/E liquid water at $-10^\circ$ C. Taken from [](doi:10.1063/1.465158
).
```

In a system composed by atoms (or particles) of different types, it is possible to define radial distribution functions between pairs of atom types. For instance, when simulating water it is common to define $g_{OO}(r)$, $g_{HH}(r)$, and $g_{HO}(r)$, which provide information on the oxygen-oxygen, hydrogen-hydrogen and hydrogen-oxygen pair correlations, respectively. See [](#fig:rdf_water) for an example.

### Mean-squared displacement

The mean squared displacement (MSD) quantifies the extent of particle diffusion within a system over time. It provides insights into the mobility of particles by tracking how far each particle moves from its original position. Formally, for a given particle $i$, the MSD at time $t$ is calculated as the average of the squared differences between its position $\vec r_i(t)$ at time $t$ and its initial position $\vec r_i(0)$, represented as $\text{MSD}_i(t) = \langle |\vec r_i(t) - \vec r_i(0)|^2 \rangle$, so that the overall MSD is

$$
\text{MSD}(t) = \frac{1}{N} \sum_{i=1}^N \langle |\vec r_i(t) - \vec r_i(0)|^2 \rangle
$$

In MD simulations, the MSD can be computed by first recording the position of each particle at each timestep, then calculating the squared displacement for each particle relative to its starting position, and finally averaging this value over all particles. If the system is ergodic, as it is often the case, we can also average over different time origins to decrease the statistical error. This is done by averaging contributions of the type $|\vec r_i(t_0 + t) - \vec r_i(t_0)|^2$ for multiple values of $t_0$.

```{figure} figures/msd.png
:name: fig:msd
:align: center
:width: 550px

The mean-squared displacement of a coarse-grained system (a tetravalent [patchy particle](https://en.wikipedia.org/wiki/Patchy_particles) model) computed at fixed density and varying ratio between the attraction strength $\epsilon$ (kept constant) and the thermal energy, $k_B T$. Taken from [](doi:10.1063/5.0209151).
```

The MSD provides essential information on diffusive behavior, as its time dependence can help distinguish between different diffusion regimes, such as ballistic ($\text{MSD}(t) \sim t^2$), brownian ($\text{MSD}(t) \sim t$), or subdiffusive (*e.g.* $\text{MSD}(t) \sim t^\alpha$, with $\alpha$ < 1), which are characteristic of different types of material and dynamical properties. [](#fig:msd) provides an example of a system that, depending on temperature, experience all these regimes.

### Root mean-squared deviation

The root mean-squared deviation (RMSD) is a measure commonly used to assess the structural deviation of a molecule or a set of particles from a reference structure, typically the initial or an experimentally determined structure. It quantifies the average distance between corresponding atoms in two structures, giving insight into conformational changes over time. The RMSD is particularly useful in studying the stability and flexibility of molecular structures, such as proteins, over the course of a simulation.

Formally, for a system with $N$ atoms, the RMSD at time $t$ relative to a reference structure is given by:

$$
\text{RMSD}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} |\vec r_i(t) - \vec r_i^{\text{ref}}|^2}
$$ (eq:rmsd)

where $\vec r_i^{\text{ref}}$ is the position of the $i$-th atom in the reference structure. To compute the RMSD in MD simulations Eq. [](#eq:rmsd) is applied after that the structure at each time point is aligned with the reference structure (typically minimizing rotational and translational differences). Formally, this can be written as

$$
\text{RMSD}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} |\hat R \cdot \vec r_i(t) + \vec t - \vec r_i^{\text{ref}}|^2},
$$ (eq:aligned_rmsd)

where $\hat R$ and $\vec t$ are the rotation matrix and translation vector that minimise the resulting RMSD, respectively. The most common use of the RMSD is to plot it as a function of time to analyse the stability of the molecule: lower RMSD values indicate that the structure remains close to the reference, while higher values suggest significant conformational changes.

```{figure} figures/rmsd_amber.png
:name: fig:rmsd_amber
:align: center
:width: 550px

The RMSD of a Alanine dipeptide simulated with [Amber](https://ambermd.org/). Taken from [this Amber tutorial](https://ambermd.org/tutorials/basic/tutorial0/index.php).
```

[](#fig:rmsd_amber) shows the RMSD of a simulation of a small molecule (an Alanine dipeptide). In this example the small value of the RMSD throughout the trajectory shows that there is no significant conformational change in the positions of the atoms relative to the starting structure.

## Tricks of the trade

(sec:neighbour_lists)=
### Neighbour lists

If the interaction potential is short-ranged, in the sense that it goes to zero at some distance $r_c$, it is clear that particles that are further away than $r_c$ will not feel any reciprocal force. If this is the case, calculating distances between all pairs of particles would be wasteful, as only a
fraction of pairs will feel a mutual interaction. There are several techniques that can be used to optimise the force calculation step (and therefore the simulation performance) by performing some kind of bookkeeping that makes it possible to evaluate only those contributions due to pairs that are "close enough" to each other. Here I will present the most common ones: cells and Verlet lists.

The idea behind cell lists is to partition the simulation box into smaller boxes, called cells, with side lengths $\geq r_c$[^cubic_cells]. This partitioning is done by using a data structure that, for each cell, stores the list of particles that are inside it. Then, during the force calculation loop, we consider that particle $i$, which is in cell $c$, can interact only with particles that are either inside $c$ or in one of its neighbouring cells (8 in 2D and 26 in 3D). The cell data structure can either be built every step, or updated after each integration step. In this latter case, the code checks whether a particle has crossed a cell boundary, and in this case it removes the particle from the old cell and adds it to the new one. This can be done efficiently with [linked lists](https://en.wikipedia.org/wiki/Linked_list). With this technique, each particle has a number of possibly-interacting neighbours that depends only on particle density and cell size, and therefore is independent on $N$. As a result, the algorithmic complexity of the simulation is $\mathcal{O}(N)$ rather than $\mathcal{O}(N^2)$. 

Differently from cell lists, in Verlet lists for each particle the code stores a list of particles that are within a certain distance $r_v = r_c + r_s$, where $r_s$ is a free parameter called "Verlet skin". Every time the lists are updated, the current position of each particle $i$, $\vec r_{i, 0}$ is also stored. During the force calculation step of particle $i$, only those particles that are in $i$'s Verlet list are considered. For a homogeneous system of density $\rho$, the average number if neighbours is

$$
N_v = \frac{4}{3} \pi r_v^3 \rho,
$$

which should be compared to the average number of neighbours if cell lists are used, which in three dimensions is

$$
N_c = 27 r_c^3 \rho.
$$

In the limit $r_s \ll r_c$, which is rather common, $r_v \approx r_c$, so that $N_c / N_v = 81 / 4 \pi \approx 6$, which means that if Verlet lists are used, the number of distances to be checked is six times smaller than with cell lists.

Verlet lists do not have to be updated at every step, or we would go back to checking all pairs, but only when when any particle has moved a distance larger than $r_s / 2$ from its original position $\vec r_{i, 0}$. Note that the overall performance will depend on the value of $r_s$, as small values will result in very frequent updates, while large values will generate large lists, which means useless distance checks on particles that are too far away to interact. Although the average number of neighbours of each particle $i$ is independent on $N$ and therefore the actual force calculation is $\mathcal{O}(N)$, the overall algorithmic complexity is still $\mathcal{O}(N^2)$, since list updating, even if not done at each time step, requires looping over all pairs. To overcome this problem it is common to use cell lists to build Verlet lists, bringing the complexity of the list update step, and therefore of the whole simulation, down to $\mathcal{O}(N)$, while at the same time retaining the smaller average number of neighbours of Verlet lists.

[^cubic_cells]: In principle, cells don't have to be cubic.

(sec:ewald)=
### Long-range interactions

:::{tip}
Most of the text of this part comes from @frenkel2023understanding, while some bits have been adapted from @schlick2010molecular.
:::

In molecular simulations, long-range interactions refer to forces that decay slowly with distance. The definition of "long-range" can be made unambiguous if we consider the general form of a pairwise interaction potential $V(r)$ between two particles separated by a distance $r$. The energy contribution of these interactions beyond a certain cut-off distance $r_c$ is given by the "tail" correction of [](#eq:U_tail). For a potential that decays as $1/r^\alpha$, the tail correction is, asymptotically,

$$
U_{\text{tail}} \sim \int_{r_c}^\infty \frac{r^2}{r^\alpha} \, dr = \int_{r_c}^\infty \frac{1}{r^{\alpha - 2}} \, dr \sim \left. \frac{1}{r_c^{\alpha - 3}}\right|_{r_c}^\infty
$$

which converges for any $r_c$ as long as $\alpha > 3$ (*e.g.* the Lennard-Jones potential, for which $\alpha = 6$)[^dipole_convergence]. By contrast, if $\alpha < 3$, the integral diverges, indicating that the energy contribution from the interaction tail remains significant when any cut-off is applied: 
the tail of the interaction potential beyond $r_c$ contributes non-negligibly to the total energy of the system, making a direct cut-off inaccurate.

I'll now quickly introduce the most popular methods used to handle long-range interactions, together with two techniques used to improve its efficiency. First, we start by considering a system consisting of $N$ charged particles in a box of volume $V = L^3$ and periodic boundary conditions. We assume that particles cannot overlap (*i.e.* that there is at least an additional short-range repulsion), and that the system is electrically neutral, *e.g.* that $\sum_i q_i = 0$. The total electrostatic energy of the system (in Gaussian units, which make the notation lighter) is

$$
U_\text{el} = \frac{1}{2} \sum_{i, j, \vec{n}}\phantom{}^{'} \frac{q_i q_j}{|\vec r_{ij} + \vec n L|},
$$ (eq:U_el)

where the prime on the summation indicates that the sum is over all periodic images $\vec n$ and over all particle pairs $(i, j)$, except $i = j$ if $\vec n = (0, 0, 0)$, *i.e.* particle $i$ interacts with all its periodic images but not with itself. Unfortunately, Eq. [](#eq:U_el) cannot be used to compute the electrostatic energy in a simulation because it contains a conditionally convergent sum.

To improve the convergence of the expression for the electrostatic potential energy, we follow [](doi:10.1002/andp.19213690304) and rewrite the expression for the charge density. In equation [](#eq:U_el) we have represented the charge density as a sum of $\delta$-functions. The contribution to the electrostatic potential due to these point charges decays as $1 / r$. Now consider what happens if we assume that every particle $i$ with charge $q_i$ is surrounded by a diffuse charge distribution of the opposite sign, such that the total charge of this cloud exactly cancels $q_i$ . In that case only the fraction of $q_i$ that is not screened contributed to the electrostatic potential due to particle $i$. At large distances, this fraction goes to $0$ in a way that depends on the functional form of the screening charge distribution, which we will take as Gaussian in the following.

```{figure} figures/ewald.png
:name: fig:ewald
:align: center
:width: 500px

The elecrostatic effect due to point charges can be seen as a sum of screend charges, minus the smoothly varying screening background. Taken from @frenkel2023understanding.
```

The contribution to the electrostatic potential at a point $r$ due to a set of screened charges can be easily computed by direct summation, because the
electrostatic potential due to a screened charge is a rapidly decaying function of $r$. However, it was not our aim to evaluate the potential due to a
set of screened charges but due to point charges. Hence, we must correct for the fact that we have added a screening charge cloud to every particle. This is shown schematically in [](#fig:ewald). This compensating charge density varies smoothly in space (because the screening charge distribution is a smoothly varying function!). We wish to compute the electrostatic energy at the site of ion $i$. Of course, we should exclude the electrostatic interaction of the ion with itself. We have three contributions to the electrostatic potential: first of all, the one due to the point charge $q_i$, secondly, the one due to the Gaussian screening charge cloud with charge $-q_i$ , and finally the one due to the compensating charge cloud with charge $q_i$. In order to exclude Coulomb self-interactions, we should not include any of these three contributions to the electrostatic potential at the position of ion $i$. However, it turns out that it is convenient to retain the contribution due to the compensating charge distribution and correct for the resulting spurious interaction afterwards. The reason we retain the compensating charge cloud for ion $i$ is that, if we do so, the compensating charge distribution is not only a smoothly varying function, but it is also periodic. Such a function can be represented by a (rapidly converging) Fourier series, and this turns out to be essential for the numerical implementation. Of course, in the end we should correct for the inclusion of a spurious "self" interaction between ion and the compensating charge cloud.

Considering screening Gaussians with width $\sqrt{2/\alpha}$, the splitting can be written as follows:

$$
\frac{1}{|\vec{r}_i - \vec{r}_j|} = \frac{\text{erfc}(\alpha |\vec{r}_i - \vec{r}_j|)}{|\vec{r}_i - \vec{r}_j|} + \frac{\text{erf}(\alpha |\vec{r}_i - \vec{r}_j|)}{|\vec{r}_i - \vec{r}_j|},
$$ (eq:ewald)

where $\alpha$ is a parameter that controls the width of the Gaussian distribution used to split the potential, and $\text{erfc}(x)$ and $\text{erf}(x)$ are the complementary error function and error function, respectively, and they come out from integrating the Gaussian screening function. In Eq. [](#eq:ewald), the first term decays quickly as $|\vec{r}_i - \vec{r}_j|$ increases, so it is computed only for nearby particles within a cut-off $ r_c$. By contrast, the second term decays slowly in real space, but can be computed efficiently in Fourier space, using a sum over reciprocal lattice vectors $\vec{k}$.

The total electrostatic energy using Ewald summation is the sum of the real-space term, the reciprocal-space term, and the self-interaction correction:

$$
E_{\text{tot}} = E_{\text{real}} + E_{\text{rec}} + E_{\text{self}}.
$$

* The real-space contribution to the total energy, $E_{\text{real}}$ is given by:

  $$
  E_{\text{real}} = \frac{1}{2} \sum_{i=1}^N \sum_{j \neq i}^N \frac{q_i q_j \, \text{erfc}(\alpha |\vec{r}_i - \vec{r}_j|)}{|\vec{r}_i - \vec{r}_j|},
  $$

  where the sum is truncated at a cut-off distance $r_c$.

* The reciprocal-space contribution $E_{\text{rec}}$ is computed using the Fourier transform of the charges and involves a sum over the reciprocal lattice vectors $\vec{k}$, which can be safely truncated at some cut-off wave vector $k_c$, since the exponential term ensures that the sum converges rapidly:

$$
E_{\text{rec}} = \frac{1}{2 V} \sum_{\vec{k} \neq 0} \frac{4 \pi}{k^2} \exp\left( -\frac{k^2}{4 \alpha^2} \right) \left| \sum_{j=1}^N q_j \exp(i \vec{k} \cdot \vec{r}_j) \right|^2.
$$ (eq:ewald_reciprocal)


* The self-interaction term, $E_{\text{self}}$, is:

$$
E_{\text{self}} = -\frac{\alpha}{\sqrt{\pi}} \sum_{i=1}^N q_i^2.
$$

Note that there are many subtleties that are linked to the boundary conditions that are applied to the system (even though the latter is infinite!). Have a look at @frenkel2023understanding and @schlick2010molecular (and references therein) for additional details.

For fixed values of $k_c$ and $r_c$, the algorithmic time scales as $\mathcal{O}(N^2)$. However, this behaviour can be improved by realising that there are values of the cut-offs that minimise the error due to the truncations. Using these values, that depend on $N$, the complexity can be brought down to $\mathcal{O}(N^{3/2})$.

A further scaling improvement can be obtained by using the so-called Particle Mesh Ewald method ([](doi:10.1063/1.464397)). In the PME method, the reciprocal space contribution to electrostatic interactions is computed through the following series of steps that involve transforming particle charges into a charge density on a grid, applying Fourier transforms, and then transforming back to obtain the forces and energies:

1. Assign each particle's charge to a regular 3D grid that spans the simulation box. This is achieved through an interpolation scheme, where each particle's charge is "spread" over multiple neighboring grid points. The most common method is B-spline interpolation (see [](doi:10.1063/1.470117) for details).
2. Once the charges are mapped onto the grid, the Fourier transform of the grid-based charge density is computed. This is done by using the Fast Fourier Transform (FFT), an efficient algorithm for performing discrete Fourier transforms, with a computational complexity of $\mathcal{O}(N \log N)$.
3. In Fourier space, the electrostatic potential for each grid point at wave vector is computed using the Ewald screening function (*i.e.* Eq. [](#eq:ewald_reciprocal)).
4. Once the reciprocal-space potential has been calculated, an inverse Fourier transform is applied to convert the potential back to real space on the grid, yielding the electrostatic potential on each grid point in the simulation box.
5. The final step is to interpolate the grid-based potential back to the positions of the actual particles, which allows the forces and energies to be computed for each particle based on the smoothed potential field. This is essentially the reverse of the initial charge assignment process, with each particle's force interpolated from the neighboring grid points.

If the cut-off of the real-space contribution is chosen so that $E_\text{real}$ can be computed in $\mathcal{O}(N)$, the overall complexity of the PME method is $\mathcal{O}(N \log N)$, which is **much** better than $\mathcal{O}(N^{3/2})$. The method has some constant overhead that makes it slower than the original Ewald sums only for small systems, which is why PME is the *de-facto* standard of biomolecular simulations.

Note that it is possible to achieve linear scaling ($\mathcal{O}(N)$) using sophisticated techniques such as the Fast Multipole Method (see *e.g.* [](doi:10.1016/0021-9991(87)90140-9), but those are rather difficult to implement, and not widely used yet.

[^dipole_convergence]: Note that if $\alpha = 3$ (*e.g.* a dipole-dipole interaction), the tail correction is $U_\text{tail} \sim [\log r]_{r_c}^\infty$, which diverges.

# Other ensembles

We know from statistical mechanics that, in the thermodynamic limit, all ensembles are equivalent. However, in simulations it can be convenient to use different ensembles, depending on the phenomena one wishes to study. I will now briefly present some methods that can be used to fix the temperature (rather than energy) and the pressure (rather than volume) in MD simulations.

[^symplectic]: it is also time-invariant and conserves volumes in phase space.

(sec:thermostats)=
## Thermostats

Before introducing some of the schemes that are used to perform Molecular Dynamics simulations at constant temperature, it is important to understand what "constant temperature" even means. From a statistical mechanical point of view, there is no ambiguity: we can impose a temperature on a system by bringing it into thermal contact with a large heat bath. Under those conditions, the distribution of the velocities is given by the Maxwell-Boltzmann distribution (Eq. [](#eq:maxwell-boltzmann)), whose second moment is connected to the temperature *via* Eq. [](#eq:equipartition). Given in this form, it is clear that the temperature also fluctuates, with the fluctuations being linked to the second and fourth moment of the Maxwell-Boltzmann distribution (*e.g.* $\propto \langle v^4 \rangle - \langle v^2 \rangle^2$). Therefore, any algorithm devised to fix the temperature of a MD simulation should reproduce not only the average $T$, but also its fluctuations. Here I will present three such algorithms.

:::{important} No more energy conservation!
As soon as a thermostat is coupled to the system, the energy will not be conserved any more. Therefore, we have one fewer way of testing the software, or the simulation parameters (*e.g.* the time step) we have chosen to use. My advice is to always turn off the thermostat and run a short simulation to see whether everything looks in order or not.
:::

(sec:andersen_thermostat)=
### Andersen

The Andersen thermostat is a widely used method in molecular dynamics simulations for controlling temperature, introduced by [Hans C. Andersen in 1980](doi:10.1063/1.439486). Its fundamental idea is to maintain a system's temperature by coupling the particles to an external heat bath through random collisions. In this approach, particles in the simulation periodically undergo stochastic collisions with a fictitious heat bath, resulting in velocity reassignment according to a Maxwell-Boltzmann distribution that corresponds to the desired temperature. This random reassignment of velocities, which can be considered as a Monte Carlo move that transports the system from one constant-energy shell to another, mimics the effect of a thermal reservoir, ensuring that the system reaches and maintains thermal equilibrium.

The thermostat has a parameter $\nu$ that is the frequency of stochastic collisions, which represents the strength of the coupling to the heat bath: by adjusting this collision frequency, the user can control how frequently the system interacts with the heat bath, thus influencing the rate at which the system equilibrates. In practice, at each time step each particle has a probability $\nu \Delta t$ of undergoing a collision, *i.e.* of being reassigned a velocity from a Maxwell-Boltzmann distribution corresponding to the target temperature. This scheme reproduced the canonical ensemble, since temperature fluctuations align with those expected from statistical mechanics. 

One drawback of the Andersen thermostat is connected to its stochastic nature: the random collisions "disturb" the dynamics in a way that is not realistic. In particular, they enhance the decorrelation of the particle velocities, thereby affecting quantities such as the diffusion constant. This effect depends on the value of $\nu$: the larger the collision frequency, the faster the decorrelation, the smaller the diffusion constant. Therefore, the Andersen thermostat should be used when we are interested in static properties only, and dynamics observables are not important.

### Nose-Hoover

The Nosé-Hoover thermostat is a more sophisticated method for controlling temperature in molecular dynamics simulations, designed to generate the canonical ensemble while preserving the continuous evolution of the system's dynamics. Unlike the Andersen thermostat, which introduces stochastic velocity reassignment, the Nosé-Hoover thermostat operates deterministically by coupling the system to a fictitious heat bath via an extended Lagrangian formalism. This allows the system to exchange energy with the thermostat in a smooth and continuous manner, maintaining temperature without introducing artificial randomness.

The extended Lagrangian formalism starts by introducing a scaling factor, $s$, that modifies the velocities of particles in the system ([](doi:10.1080/00268978400101201)). The extended Lagrangian for the Nosé-Hoover thermostat is expressed as:

$$
L = \sum_{i=1}^{N} \frac{m_i}{2} \left( \frac{ \dot{\vec r}_i }{s} \right)^2 - V(\{ \vec{r}_i \}) - g k_B T \log s,
$$

where $m_i$ and $\dot{\vec r}_i = \vec v_i$ are the mass and velocity particle $i$, $V(\{r\})$ is the potential energy of the system, and the term $ g k_B T \ln(s) $, where $ g $ is the number of degrees of freedom, introduces the necessary coupling between the system and the heat bath. The auxiliary variable $s$ is responsible for controlling the temperature by scaling the velocities of the particles so that the system's kinetic energy corresponds to the target temperature.

The dynamics of the system are then derived from this Lagrangian. The equations of motion for the positions and velocities of the particles are modified by the thermostat, resulting in:

$$
\ddot{\vec r}_i = \frac{\vec F_i}{m_i} - \zeta \dot{\vec r}_i
$$

where $\vec F_i$ is the force acting on particle $i$, and the term $\zeta \equiv \dot{s} / s$ is a friction-like coefficient that emerges from the thermostat variable. This friction term adjusts the particle velocities in response to deviations from the target temperature, ensuring that the system remains at the correct thermal equilibrium, and evolves according to:

$$
\dot{\zeta} = \frac{1}{Q} \left( \sum_{i=1}^{N} \frac{\vec p_i^2}{m_i} - g k_B T \right),
$$

where $Q$ is a parameter that controls the strength of the coupling between the system and the thermostat, $\vec p_i = m \vec r_i$ is the momentum conjugate to $\vec r_i$, so that $\sum_{i=1}^{N} \frac{\vec p_i^2}{m_i}$ is the total kinetic energy of the system, and $g k_B T$ represents the target thermal energy. If the system's kinetic energy exceeds the desired value, the variable $\zeta$ increases, effectively damping the particle velocities to bring the temperature back in line. Conversely, if the kinetic energy is too low, $\zeta$ decreases, allowing the velocities to rise and the temperature to stabilize at the target value. The "inertia" of this process is controlled by the value of $Q$, which therefore plays the role of a fictitious mass. This deterministic feedback mechanism distinguishes the Nosé-Hoover thermostat from stochastic approaches. By continuously adjusting the velocities of all particles, this method preserves the natural evolution of the system’s dynamics while still achieving temperature control. The benefit of using this extended Lagrangian framework is that it allows for smooth, continuous temperature regulation without disrupting important dynamical properties, such as diffusion coefficients or time-dependent correlations.

Note that the above equations of motion conserve the following quantity

$$
U_\text{NH} = \sum_{i=1}^N \frac{\vec p_i^2}{m_i} + V(\{ r \}) + \frac{\zeta^2 Q}{2} + g k_B T \log s,
$$

which can be used as a check when implementing the thermostat, or when testing the simulation parameters (*e.g.* the time step). [](doi:10.1103/PhysRevA.34.2499) demonstrated that the above equations of motion are unique, in the sense that are different equations of the same form (*i.e.* containing an additional friction-like parameter) cannot lead to a canonical distribution. However, the simulation samples from the canonical distribution *only* if there is a single non-zero constant of motion, the energy. If there are more conservation laws, which can happen, for instance, if there are no external forces and the total momentum of the system is different from zero, than the Nosé-Hoover method does not work any more.

```{figure} figures/nose_hoover.png
:name: fig:nose_hoover
:align: center
:width: 800px

The trajectory of the harmonic oscillator for a single initial condition simulated (from left to right) without a thermostat, with the Andersen thermostat, with the Nosé-Hoover thermostat. Adapted from @frenkel2023understanding.
```

[](#fig:nose_hoover) shows the failure of this thermostat for a cases that is deceptively simple: that of a harmonic oscillator. It is obvious that the dynamics simulated in the microcanical ensemble or with the Nosé-Hoover thermostat becomes non-ergodic. Although in practice one can always set the initial velocity of the centre of mass of the system to zero to remain with a single non-trivial conserved quantity and therefore avoid this issue, it would be better to have an algorithm that works well in the general case[^why_not_andersen]. 

As demonstrated in [](doi:10.1063/1.463940), this issue can be overcome by coupling the Nosé-Hoover thermostat to another thermostat or, if necessary, to a whole chain of thermostats, which take into account additional conservation laws. Here I provide the equations of motion for a system coupled to an $M$-link thermostat chain, where each link $i$ is an additional degree of freedom associated to a "mass" $Q_i$:

$$
\begin{align}
\ddot{\vec r}_i &= \frac{\vec F_i}{m_i} - \zeta \dot{\vec r}_i\\
\dot \zeta_k & = \frac{p_{\zeta_k}}{Q_k}\\
\dot p_{\zeta_1} &= \left( \sum_{i=1}^{N} \frac{\vec p_i^2}{m_i} - g k_B T \right) -  \frac{p_{\zeta_2}}{Q_2} p_{\zeta_1}\\
\dot p_{\zeta_k} &= \left[ \frac{p^2_{\zeta_{k-1}}}{Q_{k-1}} - k_BT \right] - \frac{p_{\zeta_{k+1}}}{Q_{k+1}} p_{\zeta_k}\\
\dot p_{\zeta_M} &= \left[ \frac{p^2_{\zeta_{M-1}}}{Q_{M-1}} - k_BT \right].
\end{align}
$$

As discussed in @frenkel2023understanding (where all these arguments are presented more in depth) of the $M$ additional degrees of freedom only two (the first link and the thermostat "centre", $\zeta_c \equiv \sum_{k=2}^M \zeta_k$) are independently coupled to the dynamics, which is what is needed when there are two conservation laws.

:::{important}
The chained Nosé-Hoover thermostat is the most common thermostat used in all-atom simulations. However, the default parameters change from package to package (at the time of writing $M = 3$ for LAMMPS and 10 for GROMACS, for instance). The choice of $Q$ (or of $\{ Q_i \}$ in the case of chains) is important and should be made with care. If $Q$ is very large, the energy exchange between the system and the reservoir is very slow, as in the $Q \to \infty$ limit we recover the Hamiltonian dynamics. By contrast, in the small-$Q$ limit the high coupling can be give raise to unphysical energy (and therefore temperature) fluctuations.
:::

[^why_not_andersen]: provided that the dynamics is of interest; otherwise, you can rely on the good old Andersen thermostat.

### Bussi-Donadio-Parrinello

We know from the equipartition theorem, Eq. [](#eq:equipartition), that the instantaneous temperature $T(t)$ is related to the kinetic energy of the system, $K(t)$, through:

$$
T(t) = \frac{2 K(t)}{3 N k_B}.
$$

Since on a computer we always have a discrete dynamics, in the following I will slightly abuse the notation by using $A(k)$ to mean $A(k \Delta t)$, where $k$ is an index that keeps track of the time step and $A$ is a time-dependent function.

The simplest way of fixing the temperature is to scale the velocities to achieve the target temperature instantaneously: the new velocities $\vec{v}_i(k + 1) $ are obtained from the old velocities $ \vec{v}_i(k) $ by multiplying them by a scaling factor $\alpha$, *e.g.* $\vec{v}_i(k + 1) = \alpha \vec{v}_i(k)$, where $\alpha$ is given by

$$
\alpha = \sqrt{\frac{T}{T(k)}} = \sqrt{\frac{K}{K(k)}}
$$ (eq:v_rescaling)

where $T$ and $K$ are the target temperature and kinetic energy. This method has two drawbacks:

1. It does not yield the correct fluctuations for the kinetic energy.
2. The rescaling affects all particles at the same time in a abrupt way, greatly disturbing the dynamics, as was the case with the Andersen thermostat.

The second issue can be mitigated by weakly coupling the system to a heat bath with a characteristic relaxation time $\tau$. With this method, called the [Berendsen thermostat](doi:10.1063/1.448118), the velocities are gradually adjusted to bring the temperature toward the target value. The velocity scaling factor in the Berendsen thermostat is given by:

$$
\alpha = \sqrt{1 + \frac{\Delta t}{\tau} \left( \frac{K}{K(k)} - 1 \right)}.
$$ (eq:berendsen)

This factor ensures that the system's temperature approaches $T$ over time, with the relaxation time $\tau$ controlling the speed of this adjustment. The Berendsen thermostat achieves smooth temperature control, but it still suppresses the natural energy fluctuations required for proper sampling of the canonical ensemble.

[](doi:10.1063/1.2408420) introduced a thermostat that is based on the idea of velocity rescaling, but it samples the canonical ensemble. Note that in the continous limit, *i.e.* for $\Delta t \to 0$, Eq. [](#eq:berendsen) becomes

$$
dK = (K - K(t)) \frac{dt}{\tau}.
$$

This equation makes it obvious that $K(t) \to K$ exponentially, in a deterministic way (*i.e.* without fluctuations). A way of introducing fluctuations is to add a Weiner noise $dW$, weighted by a factor that ensures that the fluctuations of the kinetic energy are canonical[^wiener_factor]:

$$
dK = (K - K(t)) \frac{dt}{\tau} + 2 \sqrt{\frac{K(t)K}{N_f}} \frac{dW}{\sqrt{\tau}},
$$ (eq:bussi)

where $N_f$ is the number of degrees of freedom in the system. Since the total momentum is conserved, $N_f = 3N - 1$ in a 3D simulation. With some non-trivial math [](doi:10.1063/1.2408420) demonstrated that Eq. [](#eq:bussi) yields the following scaling factor:

$$
\begin{align}
\alpha^2 &= e^{-\Delta t / \tau} + \frac{K}{N_f K(k)} \left( 1 - e^{-\Delta t / \tau} \right) \left( R_1^2 + \sum_{i = 2}^{N_f} R_i^2 \right)\\
& + 2 R_1 e^{-\Delta t / \tau} \sqrt{\frac{K}{N_f K(k)} \left( 1 - e^{-\Delta t / \tau} \right) },
\end{align}
$$

where the $\{ R_i \}$ are independent random numbers extracted from a Gaussian distribution with zero mean and unit variance.

The interpretation of [](#eq:bussi) is that, in order to sample the canonical ensemble by rescaling the velocities, the target kinetic energy (*i.e.* the target temperature) should be chosen randomly from the associated canonical probability distribution function[^kinetic_PDF]:

$$
P(K(t)) \propto K(t)^{N_f / 2 - 1} e^{-\beta K(t)}.
$$ (eq:K_PDF)

One way of doing this would be to randomly extract a value $\bar K$ from Eq. [](#eq:K_PDF) each time a rescale should be performed, and use it as the target kinetic energy in Eq. [](#eq:v_rescaling). However, this would generate abrupt changes in the dynamics, as for the simple velocity rescaling and Andersen thermostats. By contrast, the stochastic process defined in Eq. [](#eq:bussi) is smoother, since the rescaling procedure happens on a timescale given by $\tau$, in the same vein as with the Berendsen thermostat.

As for the Nosé-Hoover thermostat, this thermostat, which is sometimes called the stochastic velocity rescaling thermostat, has an associated conserved quantity that can be used to check the choice of the time step and of the other simulation parameters:

$$
U_\text{BDP} = \sum_{i=1}^N \frac{\vec p_i^2}{m_i} + V(\{ r \}) - \int_0^t (K - K(t')) \frac{dt'}{\tau} - 2 \int_0^t \sqrt{ \frac{K(t')K}{N_f}}. \frac{dW(t')}{\sqrt{\tau}}
$$

[^wiener_factor]: It turns out that there is some freedom in choosing this factor, but only if we don't want to keep the Berendsen part untouched.
[^kinetic_PDF]: This is a chi-squared PDF, which appears when dealing with sums of the squares of independent Gaussian random variables. Since the kinetic energy is the sum of squared velocities, each of which is distributed normally, its PDF is a chi-square.

### Other thermostats

There are many other thermostats, which I will not introduce here because of time constraints. Here I list a few, just for reference:

* Langevin thermostats for systems where the dynamics is overdamped. Examples are particles dispersed in a viscous fluid.
* Brownian thermostats for systems where inertia is negligible. Examples are concentrated suspensions of colloids or bacteria.
* Stochastic rotational dynamics to model hydrodynamic interaction. Examples are diluted solutions of colloids or bacteria where long-range hydrodynamic effects are important and of interest (*e.g.* sedimentation or flocking).

:::{important}
I reiterate that in equilibrium simulations the static properties of a system are independent on the thermostat used (provided that the thermostat reproduces the correct temperature fluctuations). However, if the dynamics is of interest (which is always the case in out-of-equilibrium systems, but it is often the case also for equilibrium systems), then choosing the right thermostat is [very important](doi:10.1140/epje/i2018-11689-4).
:::

## Barostats

Experiments are usually performed at constant pressure (*e.g.* ambient pressure). Moreover, while statistical mechanics results ensure that ensembles are equivalent in equilibrium, there are many situations in which it is more convenient to simulate in a specific ensemble. For instance, equilibrating a crystal structure usually requires that the edge lengths of the box can relax and fluctuate, or the coexisting region between a gas and a liquid is extended in $(T, \rho)$, but it is just a line in the $(T, P)$ projection. In the same as controlling the temperature requires coupling the system to a thermostat, fixing $P$ means coupling the system to a barostat. When simulating at constant pressure it is very common, but by no means the only possibility, to also fix the temperature, thereby simulating in the so-called isothermal-isobaric ensemble.

As for thermostats, there exist many barostats, each having its own unique advantages and disadvantages. Here I briefly present three barostats, two based on the extended-Lagrangian formalism and a stochastic one.

### Nosè-Hoover barostat

The Nosè-Hoover formalism can be extended to fix the pressure. In this case the Equations of motion, as presented in [](doi:10.1063/1.467468), which improves the original scheme by [Hoover](doi:10.1103/PhysRevA.34.2499) that only yields approximated constant-$P$ distributions, are:

$$
\begin{align}
\dot{\vec r_i} &= \frac{\vec p_i}{m_i} + \frac{p_\epsilon}{W} \vec r_i\\
\dot{\vec p_i} &= \vec F_i - \left( 1 + \frac{d}{dN} \right) \frac{p_\epsilon}{W} \vec p_i - \frac{p_\zeta}{Q} \vec p_i\\
\dot{V} &= \frac{dVp_\epsilon}{W}\\
\dot{p_\epsilon} &= dV (P(t) - P) + \frac{1}{N} \sum_{i=1}^N \frac{\vec p_i^2}{m_i} - \frac{p_\zeta}{Q}\vec p_\epsilon,
\end{align}
$$ (eq:nose-hoover_barostat)

where $d$ is the space dimensionality, $\epsilon = \log(V / V(0)$, where $V(0)$ is the volume at $t = 0$, $W$ is the inertia parameter associated to $\epsilon$, $p_\epsilon$ is the conjugate momenum of $\epsilon$ (here I'm using the same notation as @frenkel2023understanding), and $P$ and $P(t)$ are the target and instantaneous pressures, respectively. In particular, the instantaneous pressure is

$$
P(t) = \frac{1}{dV} \left[ \sum_{i=1}^N \left( \frac{\vec p_i^2}{m_i} + \vec r_i \cdot \vec F_i \right) - dV \frac{\partial U(V)}{\partial V} \right],
$$ (eq:pressure_NPT)

which is the virial pressure of a system with a non-constant volume. Note that the second term in [](#eq:pressure_NPT) is non-zero if the energy depends explicitly on volume, which is always the case for interaction potentials that are long-range, or for which long-range corrections due to cut-offs have to be evaluated.

Note that equations [](#eq:nose-hoover_barostat) also contain a coupling to a Nosé-Hoover thermostat, which can be generalised to a chain thermostat if $Q \to Q_1$, $p_\zeta \to p_{\zeta_1}$ the equations of motions of the other $M - 1$ links are added.

### Parrinello-Rahman

The Parrinello-Rahman barostat, introduced in [](doi:10.1103/PhysRevLett.45.1196), uses an extended Hamiltonian to account for both particle motions and cell shape fluctuations, enabling the simulation of anisotropic pressure conditions, which is very common, for instance, when dealing with crystalline structures. This Hamiltonian includes original Hamiltonian of the system, which accounts for the kinetic and potential energy of the particles, plus additional terms that describe the kinetic and potential contributions of the simulation box itself.

The additional terms can be written in terms of the matrix

$$
\hat H = \begin{bmatrix}\vec a & \vec b & \vec c\end{bmatrix} = \begin{bmatrix} a_x & b_x & c_x \\ a_y & b_y & c_y \\ a_z & b_z & c_z \end{bmatrix},
$$

where $\vec a$, $\vec b$ and $\vec c$ are the basis vectors of a generic (possibly anisotropic) simulation box, and coincide with its edges. The position of particle $i$ in the box is just

$$
\vec r_i = \hat H \vec s_i = x_i \vec a + y_i \vec b + z_i \vec c,
$$

where $0 < x_i, y_i, z_i < 1$ are the scaled (fractional) coordinates of particle $i$. With this formalism, the distance between particles $i$ and $j$ is

$$
r^2_{ij} = \vec s_{ij} \hat H^T \hat H \vec s_{ij} = \vec s_{ij} \hat G \vec s_{ij},
$$

where $G \equiv \hat H^T \hat H$.

We can now write down the expression for the extended Hamiltonian:

$$
H = \sum_{i=1}^N \frac{1}{2} m_i \dot{\vec{s}}_i^T \hat G \dot{\vec{s}}_i + V(\{ \vec r_i \}) + \frac{1}{2} W \text{Tr}(\dot{\hat H}^T \cdot \dot{\hat H}) + P(t) \det(\hat H),
$$ (eq:H_parrinello_rahman)

where $\det(\hat H) = V$ is the box volume, $P(t)$ is the istantaneous pressure, and the parameter $W$ is the strength of the barostat coupling, and can be interpreted as the mass of the fictitious piston exerting the external pressure on the system.

By deriving Eq. [](#eq:H_parrinello_rahman), the following equations of motion are obtained:

\begin{align}
m_i \ddot{s}_i &= \hat H^{-1} \vec f_i - m_i \hat G^{-1} \dot{\hat G} \dot{\vec s}_i\\
W \ddot{\hat H} & = (P(t) - P)V(\hat H^{-1})^T,
\end{align}

where $\vec f_i = - \vec \nabla_{\vec r_i} V(\{ \vec r_j \})$ is the force acting on particle $i$.

### Stochastic cell rescaling

```{warning}
TODO
```

# Classical force fields

:::{tip}
The main references for this part are @schlick2010molecular and @leach2001molecular.
:::

Now that we know the algorithms used to run MD codes, we shift our attention to the interaction potentials acting between the components that make up the system we want to simulate. We start from the so-called classical (or empirical) force fields, which describe the interactions between atoms using simplified mathematical models based on quantum mechanics calculations, empirical observations, and physical principles from classical mechanics. Unlike quantum methods (introduced in [](./quantum.md)), which rigorously account for the wave-like nature of particles and their interactions through principles like superposition and entanglement, classical force fields offer a pragmatic approach that balances computational feasibility with accuracy.

The potential energy functions used in classical force fields typically consist of terms that describe various types of interactions, including covalent bonds, van der Waals (dispersion) forces and electrostatic interactions. Despite their simplicity, classical force fields can exhibit remarkable predictive power and have been successfully applied across diverse scientific disciplines. Common force fields used in biomolecular simulations are [AMBER](https://ambermd.org/AmberModels.php) and [CHARMM](https://mackerell.umaryland.edu/charmm_ff.shtml).

In a classical force field (FF from now on), the main assumption is that of *additivity*, which makes it possible to write down the total energy of a system as a sum of several contributions. These are classified into terms that account for atoms that are linked by covalent bonds (*bonded* or *local* terms), and terms that account for noncovalent interactions (*non-bonded* or *non-local* terms), so that the total energy of the system is written as

$$
E_{\rm tot}(\dofs) = E_{\rm bonded}(\dofs) + E_{\rm non-bonded}(\dofs),
$$ (eq:tot-energy)

where I made explicit the fact that the energy (and its contributions) are functions of the set of positions of the atoms, $\dofs$.

```{note}
The separation of contributions into bonded and non-bonded terms is not only a conceptual one, but it has also practical implications. For instance, bonded interactions act on pairs or groups of atoms that do not change during the course of the simulation, and therefore the data structures that keep track of these interacting atoms do not need to be updated. Moreover, non-bonded interactions are usually much less steep than bonded ones, which makes it possible to implement multiple-timestep protocols to speed-up simulations (see *e.g.* Chapter 14 of @schlick2010molecular).
```

The functional forms of the functions (implicitly) defined in Eq. [](#eq:tot-energy) depend on the force field employed, and are parametrised so as to reproduce the experimental structure and dynamics of target molecular systems, but also leveraging *ab initio* calculations, which are often used to complement the experimental results.

```{figure} figures/ff_spectra.png
:name: fig:ff_spectra
:align: center
:width: 700px

Characteristic (left) stretching and (right) bending and torsional vibrational frequencies extracted from experimental spectra. Taken from @schlick2010molecular.
```

The traditional sources of parameters are provided by vibrational spectra (see [](#fig:ff_spectra) for examples), which have been routinely used to derive the force constants of many of the functional forms discussed below. The idea is that the frequencies extracted from experimental spectra can be compared with the characteristic frequencies of the normal modes as evaluated in simulations in order to find the values of the force constants that minimise the differences. 

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

(sec:bond_stretching)=
### Bond stretching

This interaction occurs when atoms connected by a covalent bond move closer together or farther apart, resulting in the stretching or compression of the bond, and can be considered as "excitation" terms that accounts for small deviations from reference values, which are usually taken from experimental measurements.

The simplest type of interaction that accounts for bond deformations is the harmonic potential, which is based on Hooke's law and takes a quadratic form:

$$
E_\text{harmonic}(r) = \frac{1}{2}k (r - r_\text{ref}),
$$ (eq:harmonic_bond)

where $k$ is the force constant (which can be estimated by the, possibly reduced, mass and frequency of a bond vibration through $k = m \omega^2$), and $r_\text{ref}$ is the reference value. Eq. [](#eq:harmonic_bond) works only for rather small deformations ($\approx 0.1 \angstrom$).

Going beyond small deformations requires more complicated functional forms. An example is the Morse potential:

$$
E_\text{Morse}(r) = D [1 - \exp(-S_m(r - r_\text{ref})]^2,
$$ (eq:morse)

where the constants $S_m$ and $D$ controls the width and the depth of the potential well. This potential goes to infinity for $r \to 0$, while it tends to $D$ as $r \to \infty$, modelling the bond dissociation. Evaluating exponential functions is a rather slow business on a computer, which is why the Morse potential is often expanded in series, and only the first terms are retained. The expansion up to the fourth order reads

$$
E_\text{Morse}(r) \approx DS_m^2(r - r_\text{ref})^2 - DS_m^3(r - r_\text{ref})^3 + \frac{7}{12}DS_m^4(r - r_\text{ref})^4.
$$ (eq:morse_expansion)

I note in passing that the previous equation implies the force constant of Eq. [](#eq:harmonic_bond) can be related to the Morse parameters by $k = DS_m^2$.

```{figure} figures/bond_stretching.png
:name: fig:bond_stretching
:align: center
:width: 500px

Morse, harmonic, cubic and two quartic bond potentials for H-Br. Taken from @schlick2010molecular.
```

A comparison between the Morse potential and its expansions at the second, third and fourth order is shown in [](#fig:bond_stretching). The change of curvature of the cubic potential yields unphysical behaviour for large deformations, but even the other approximated forms cannot model bond dissociation. Finally, some force fields use the following different (special) quartic form to avoid having to compute odd powers of $r$ (which require the computation of a square root, which is an expensive function to compute):

$$
E_\text{quartic} = S_q(r^2 - r_\text{ref}^2)^2,
$$

where the constant can be linked to the Morse parameters by matching the two second-order expansions, yielding $S_q = DS_m^2 / 4 r_\text{ref}^2$. This function is also shown in [](#fig:bond_stretching).

### Angle bending

When three atoms are connected by two consecutive covalent bonds, the angle between these bonds can change, causing the atoms to deviate from a linear configuration. Angle bending interactions describe the energy associated with such deformations and are often modeled using harmonic potentials, where the energy increases quadratically with the deviation of the angle from its equilibrium value. At first approximation, the reference (equilibrium) value is given by the type of orbital hybridisation due to the bonds (*e.g.* $180^\circ$ for *sp*, $120^\circ$ for *sp^2* and $109.47^\circ$ for *sp^3*). However, the orbitals are often deformed, and real values differ from ideal ones. For instance, in propane the $C-C-C$ and $H-C-H$ bond angles are $\approx 112.5^\circ$, and $\approx 107.5^\circ$, respectively.

The most common potentials used have a harmonic form involving either angles or cosines:

(eq:bond_angle)=
\begin{align}
E^\theta_\text{harmonic}(\theta) & = K_h(\theta - \theta_\text{ref})^2\\
E^\theta_\text{trig}(\theta) & = K_t(\cos \theta - \cos \theta_\text{ref})^2.
\end{align}

If the second equation is expanded in series around $\theta_\text{ref}$ and compared to the first equation, it is possible to obtain the following relation, which connects $K_h$ and $K_t$ so that the small-angle fluctuations of the two forms are similar:

$$
K_t = K_h \sin^2 \theta_\text{ref}.
$$ (eq:bond_angle_constants)

Since it does not require the computation of inverse trigonometric functions, the trigonometric form is more efficient from the computational point of view. [](#fig:bond_angle) shows the two forms defined in Eq. [](#eq:bond_angle), where the value $K_t$ has been chosen according to Eq. [](#eq:bond_angle_constants).

```{figure} figures/bond_angle.png
:name: fig:bond_angle
:align: center
:width: 500px

Harmonic bond-angle potentials of the forms given in Eq. [](#eq:bond_angle) for an aromatic $C-C-C$ bond angle ($CA-CA-CA$ atomic sequence in CHARMM) with parameters $K_h = 40$ kcal/(mol rad$^2$) and $\theta_\text{ref} = 2.1$ rad ($120^\circ$). The $K_t$ force constant is evaluated by using Eq. [](#eq:bond_angle_constants). Taken from @schlick2010molecular.
```

### Dihedral rotation

Dihedral or torsional interactions arise when four atoms are connected by three consecutive covalent bonds, forming a torsion angle or dihedral angle. Rotating one group of atoms around the axis defined by the central bond changes the conformation of the molecule, leading to different energy minima corresponding to different dihedral angles.

As noted [before](./proteins.md#sec:molecular_vibrations), in proteins the most important dihedral (torsional) angles are those associated to rotations around the $N - C^\alpha$ and $C^\alpha - C$ bonds, $\phi$ and $\psi$, which feature free-energy barriers of the order of the thermal energy. However, rotations around the peptide bond (dihedral $\omega$) and around sp$^3$-sp$^3$ bonds such as those found in aliphatic said chains (dihedral $\chi$) can also play a role in dictating the flexibility of proteins.

Dihedral interactions are typically described using periodic potentials that capture the periodicity of the energy as a function of the dihedral angle. A generic torsional interaction associated to a dihedral angle $\varphi$ takes the form

$$
E_t(\varphi) = \sum_{n} \frac{V_n}{2}[1 + \cos(n\varphi - \varphi_0)],
$$ (eq:dihedral_potential)

where $n$ is an integer that determines the periodicity of the barrier of height $V_n$, and $\varphi_0$ is a reference angle, which is often set to $0$ or $\pi$.

```{figure} figures/ff_dihedrals.png
:name: fig:dihedral_potential
:align: center
:width: 700px

Twofold and threefold torsion-angle potentials and their sums for an $O-C-C-O$ rotational sequence in nucleic-acid riboses ($V_2 = 1.0$ and $V_3 = 2.8$ kcal/mol) and a rotation about the phosphodiester ($P-O$) bond in nucleic acids ($V_2 = 1.9$ and $V_3 = 1.0$ kcal/mol).  Taken from @schlick2010molecular.
```

The most common potentials of the form of [](#eq:dihedral_potential)) are those comprising twofold and threefold terms, which are enough to reproduce common energy differences (*e.g.* cis/trans and trans/gauce). Two examples from the CHARMM force field are presented in [](#fig:dihedral_potential):

* The rotational interaction in the $O-C-C-O$ sequence in nucleic acids (*e.g.*, $O3'-C3'-C2'-O2'$ in ribose) shows a minimum at the trans state.
* The rotation about the phosphodiester bond ($P-O$) in nucleic acids shows a very shallow minimum at the trans state. 

For smaller molecules is often necessary to include higher-order terms ($n = 1, 2, 3, 4$ and $6$) to reproduce the experimental torsional frequencies accurately.

### Improper rotation

Improper torsion interactions occur when four atoms are bonded in a way that one atom is not directly in line with the other three. This configuration is often necessary to maintain the correct geometry or chirality of a molecule. Improper torsions are used to model the energy associated with deviations from the ideal geometry. The potential energy associated with improper torsions is typically described by a harmonic potential or a cosine function, depending on the force field parameters and the desired behavior of the molecule. A common form is

$$
E_\text{imp}(\chi) = \frac{V_\text{imp}}{2} \chi^2,
$$

where $\chi$ is the improper Wilson angle, which has the following definition: given four atoms $i, j, k$ and $l$, where $j$ is bonded to the other three, $\chi$ is the angle between bond $j-l$ and the plane defined by $i, j$ and $k$, as shown in [](#fig:improper).

```{figure} figures/improper.png
:name: fig:improper
:align: center
:width: 500px

The definition of the Wilson angle $\chi$, as seen (left) from the top and (right) from the side.
```

(sec:cross_terms)=
### Cross terms

```{figure} figures/ff_crossterms.png
:name: fig:ff_crossterms
:align: center
:width: 500px

Schematic illustrations for various cross terms involving bond stretching, angle bending, and torsional rotations. Here UB stands for Urey-Bradley, which is a way of modelling the stretch/bend coupling by taking into account only the distance between the "far" atoms in the triplet. Taken from @schlick2010molecular.
```

It is sometimes necessary to include terms that couple different degrees of freedom in order to improve the accuracy, especially in force fields used to model small molecules. These *cross terms* are usually simple (quadratic) functions of the quantities in play (*e.g.* atom-atom distances or bond angles) and are fitted to the experimental vibrational spectra. [](#fig:ff_crossterms) pictorially shows some of these contributions.

## Non-bonded interactions

*Non-bonded interactions* refer to the energy contributions that arise between atoms or molecules that are not directly connected by chemical bonds, and comprise several terms. Here I present the most common ones. Note that non-electrostatic non-bonded interactions are usually short-ranged and therefore are cut-off at some distance to improve performance, as discussed in [](#sec:cut-off).

### Van der Waals Interactions

Van der Waals interactions are weak forces that arise due to fluctuations in the electron density of atoms or molecules. These interactions include both attractive forces, arising from dipole-dipole interactions and induced dipole-induced dipole interactions (van der Waals dispersion forces, see @israelachvili2011intermolecular for a derivation of these terms), and repulsive forces, resulting from the overlap of electron clouds at close distances. Van der Waals interactions are described by empirical potential energy functions. The most common forms are the Lennard-Jones and Morse potentials, which we have already encountered when discussing [interactions in proteins](#sec:van-der-waals) and [](#sec:bond_stretching), respectively. Here I report their functional forms, both of which account for attractive and repulsive components, for your convenience:

\begin{align}
V_\text{LJ}(r) &= 4 \epsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right)\\
V_\text{Morse}(r) &= \epsilon [1 - \exp(-a(r - r_0)]^2,
\end{align}

where $\epsilon$ sets the depth of the attractive well, $\sigma$ is the LJ diameter, $r_0$ is the position of the minimum of the Morse potential, and $a$ is linked to the curvature close to such a minimum.

### Electrostatic Interactions

Electrostatic interactions arise from the attraction or repulsion between charged particles, such as ions or polar molecules. These interactions are governed by Coulomb's law, which states that the force between two charged particles is proportional to the product of their charges and inversely proportional to the square of the distance between them. In molecular systems, electrostatic interactions include interactions between charged atoms or molecules (ion-ion interactions), as well as interactions between charged and neutral atoms or molecules (ion-dipole interactions or dipole-dipole interactions). Since electrostatic interactions are long-ranged, they have to be handled with care, as discussed in [](#sec:ewald).

### Hydrogen Bonding

While sometimes hydrogen bonding interactions are modeled using empirical *ad-hoc* potentials that account for the directionality and strength of hydrogen bonds, in modern force fields they arise spontaneously as a result of the combination of the Van der Waals and electrostatic interactions acting between the atoms.

## More on force fields

Modern force fields contain tens, hundreds or even more parameters, depending on what they have been designed to model. The value of every single parameter has been optimised to reproduce some target quantity, either generated with higher-fidelity numerical methods, or measured experimentally, or a mixture of both. However, no force field is perfect, as each FF has its own advantaged and disadvantages.

The force field should be chosen carefully, based on the problem at hand. For proteins and nucleic acids, common FFs are [AMBER](https://ambermd.org/AmberModels.php) and [CHARMM](https://mackerell.umaryland.edu/charmm_ff.shtml), but there are other possibilities. Note that it is often necessary to simulate molecules or functional groups that are not supported by the specific FF we wish to use. In this case, it is possible to use multiple force fields, provided the two FFs are compatible, *i.e.* someone has parametrised the "cross interactions" between the atoms or molecules modelled with distinct force fields. As a golden rule, you should **never** mix parameters from different force fields.

In general, you should always follow the advice and good-practices given in the documentation of the FF(s) you are using, unless you **really** know what you are doing. For instance, the CHARMM force field recommends to use the TIP3P model for water rather than more sophisticated and realistic ones, since it has been parametrised using that particular model. Using another model can, in principle, be made to work, but would also make your results somehow questionable, which is **not** what you want from a scientific point of view.

## GROMACS

GROMACS, an acronym for Groningen Machine for Chemical Simulations, is an open-source software package designed primarily for molecular dynamics simulations of biomolecular systems, providing insights into the structural dynamics of proteins, lipids, nucleic acids, and other complex molecular assemblies. GROMACS can run simulations efficiently on a wide range of hardware platforms, from single processors to large parallel computing clusters, and implements advanced algorithms and techniques, such as domain decomposition, particle-mesh Ewald summation for long-range electrostatics, and multiple time-stepping schemes. 

Interestingly (and differently from other simulation engines) GROMACS also offers a comprehensive suite of analysis tools, making it possible to perform tasks such as trajectory analysis, free energy calculations, and visualization of simulation results.

:::{warning} 
The text that follows is just a collection of notes. It will be updated in the future.
:::

* Use `mdp` files to create a `tpr` (which is a binary, non-human-readable file): at this step we need a topology file (`topol.top`)[^topol_file]
* Index files (`.ndx`) are used to assign atoms to categories that can be used to easily find them. Use `gmx make_ndx -f file.tpr -o output` to make a new index file.
* The `-deffnm` flag tells Gromacs to use the input file as a template for the filenames of output quantities
* The default trajectory file (`.trr`) contains all the info (coordinates, velocities, *etc.*). By contrast, `xtc` files are lighter since they store only the atom coordinates.

```{warning} Berendsen thermostat
On newer versions, Gromacs will spit out a warning if you use the Berendsen thermostat (`tcoupl = berendsen`). However, by default `gmx grompp` considers a single warning as a fatal error. Use the `-maxwarn` flag to raise the number of acceptable warnings (*e.g.* `-maxwarn `).
```

[^topol_file]: in Gromacs lingo, a topology file contains the details of the force field

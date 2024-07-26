---
title: A simple molecular dynamics code
---

:::{attention}
Remember to read the [homeworks](#sec:homeworks) section to understand what you need to submit for this assignment.
:::

# The main assignment

Molecular dynamics code integrate Newton's equations to follow the dynamics of a systems made of components interacting with each other through a well-defined set of potentials. The mandatory part of this assignment is to write a molecular dynamics code that can be used to simulate a 3D Lennard-Jones system at a given temperature $T$ and density $\rho = N / V$, where $N$ is the number of particles and $V = L^3$ is the volume of the simulation box. The code should 

1. Take as input at least the number of particles, temperature, density and the integration time step, $\delta t$.
2. Use the Velocity-Verlet integrator for the equations of motion.
3. Take into account periodic boundary conditions (see [below](#sec:pbc) for details).
4. Make it possible to optionally couple a thermostat to the system. The simplest to implement is the [Andersen thermostat](#sec:andersen_thermostat).

Make use of reduced units: the particle diameter is $\sigma$, the depth of the LJ well is $\epsilon$, and the mass of a particle is $m$, so that time is measured in units of $t_0 = \sigma \sqrt{m / \epsilon}$.

The documentation accompanying the code should contain a plot that shows that the fluctuations of the total energy is proportional to $\delta t$[^range_dt].

[^range_dt]: This will be true only if $\delta t$ is not too large: be careful.

# Possible extensions

1. Add [neighbour lists](#sec:neighbour_lists) to speed up computation.
2. Add a routine to calculate the [pressure](#sec:compute_pressure) of the system. See [below](#sec:compare_pressure) for values you can use to check your code.
3. Extend your code to also support Lennard-Jones *dumbbells*. Now your basic particle is made of two spheres connected by a *bonded* interaction. You can check that your code works by comparing the pressure as a function of temperature and density to Fig. 2 of [this paper](doi:10.1103/PhysRevE.107.034607)[^packing_fraction], where the bonded interaction acting between the two spheres that compose a dumbbell is $V_{\rm bond} = k(r - \sigma)^2$, with $k = 3000 \epsilon / \sigma^3$.

[^packing_fraction]: In the figure, the packing fraction $\eta$ is the density multiplied by the volume of a dumbbell, *e.g.* $\eta = \rho \pi \sigma^3 / 3$, where $\sigma$ is the diameter of a Lennard-Jones sphere.

# Additional details

(sec:pbc)=
## Periodic boundary conditions

Applying the periodic boundary conditions in a cubic box with side $L$ means that the distance between two particles with coordinates $\vec{r}_i$ and $\vec{r}_j$ is

$$
\vec{r} = (\vec{r_2} - \vec{r_1}) - {\rm round}\left(\frac{\vec{r_2} - \vec{r_1}}{L}\right)L,
$$

where $\rm round()$ is a function that rounds its argument to the nearest integer.

(sec:compare_pressure)=
## Pressure in a Lennard-Jones system

There are many (mostly old) papers dealing with estimating the equation of state of a Lennard-Jones fluid. In most of them, the pressure computed in an MD simulation is corrected by the following term, which takes into account the truncated tail of the potential in a continuous way:

$$
\Delta P = \frac{32}{9} \rho^2 \left( \left(\frac{\sigma}{r_c}\right)^{9} - \frac{3}{2} \left( \frac{\sigma}{r_c}\right)^{3} \right).
$$ (eq:P_correction)

With this in mind, you can try to evaluate the pressure, correct it by applying the shift of Eq. [](#eq:P_correction), and compare the results with those reported in [](doi:10.1080/00268979300100411). However, note that in the paper a large cut-off ($r_c = 4.0\, \sigma$) was used. 

If you feel that your code is not fast enough, you can use the following values, obtained with $\rho\sigma^3 = 0.60$, $r_c = 2.5 \, \sigma$, $\delta t = 0.001 \, t_0$ and simulations of $10^6$ time steps:

| $T$ [$\epsilon / k_B$] | $P$ [$k_B T / \sigma^3$] | $P_{\rm corr}$[^P_corrected] [$k_B T / \sigma^3$]|
| :---: | :---:| :---:|
| 2.0 | 2.2 | 1.8 |
| 3.0 | 4.0 | 3.6 |
| 4.0 | 5.7 | 5.3 |

:::{warning}
Here I run short simulations to generate numbers that can be comparable with your. Such short runs do not allow to compute the pressure with high precision, nor to provide a sensible estimate of the error, which nevertheless should be of the order of the last reported digit.
:::

[^P_corrected]: The pressure corrected with Eq. [](#eq:P_correction).

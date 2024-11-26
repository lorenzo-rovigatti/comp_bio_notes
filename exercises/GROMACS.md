---
title: Simulating a polypeptide with GROMACS
---

:::{attention}
Remember to read the [homeworks](#sec:homeworks) section to understand what you need to submit for this assignment.
:::

# The main assignment

Perform all-atom molecular dynamics simulations of an aqueous solution of a poly-alanine with an a-helix conformation at $T = 290$ K using the CHARMM force field. In the work you should:

1. Use the minimized structure (`mini.gro`) of an aqueous solution of poly-alanine equilibrate the system in temperature ($T = 290$ K) and pressure ($P = 1$ bar), following a procedure similar to what was done at $T = 350$ K during the tutorial. You should plot the time evolution of the temperature during the $NVT$ simulation, and as well as the time evolution of the density during the $NPT$ simulation to check the system equilibration.
2. Run 10 ns of simulation at $290$ K and check the temperature and density during the data acquisition.
3. Calculate the root mean-squared deviation and the radius of gyration of the poly-alanine as a function of time and compare it to the behavior observed at 350 K and discuss the plot.

<!-- # Possible extensions

Perform all-atom molecular dynamics simulations of an aqueous solution of a poly-alanine at three different temperature values in the range between $290$ K < $T$ < $350$ K. Calculate the root mean-squared deviation and the radius of gyration of the poly-alanine as a function of time and compare the behavior observed at different temperatures. Plot the values of the radius of gyration averaged over the 10 ns of simulation time as a function of temperature.-->

# Additional details

The files required to run the simulation, together with some instructions about how to install VMD and GROMACS, can be found at [this link](https://elearning.uniroma1.it/pluginfile.php/1413675/mod_label/intro/tutorial_poly-alanine_x_stud.zip). Note that you have to be logged in as a Sapienza user to be able to download the file.

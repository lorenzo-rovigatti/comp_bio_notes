---
title: The hydrophobic-polar protein lattice model
---

:::{attention}
Remember to read the [homeworks](#sec:homeworks) section to understand what you need to submit for this assignment.
:::

# The main assignment

# Possible extensions

# Additional details

The Hydrophobic-Polar (HP) protein lattice model is a simplified computational approach used to investigate the principles of protein folding and structure. In this model, proteins are abstracted as sequences composed solely of hydrophobic (H) and polar (P) amino acids. This binary representation simplifies the complex nature of real proteins, allowing researchers to focus on the primary forces driving the folding process.

In the HP model, the protein sequence is mapped onto a grid, or lattice, which can be two-dimensional (2D) or three-dimensional (3D). Each amino acid in the sequence occupies a point on this grid, and consecutive amino acids are connected by edges representing peptide bonds. Common lattice types used include square lattices for 2D models and cubic lattices for 3D models, but other choices are possible. This lattice representation constrains the possible configurations of the protein, mimicking the spatial restrictions that occur in real proteins.

The HP model is widely used in computational biology for studying protein folding mechanisms, developing algorithms for protein structure prediction, and understanding the impact of mutations on protein stability and function. Despite its simplicity, the HP model effectively captures essential features of protein folding, providing valuable insights into the roles of hydrophobic and polar interactions in determining protein structure and stability. By reducing the complexity of protein structures to fundamental interactions, the HP model serves as a powerful tool for exploring the basic principles underlying protein behavior.

The folding of the protein is simulated by arranging the sequence on the lattice to achieve a conformation that minimizes the system's energy, which quantifies the thermodynamic stability of the specific conformation. Typically, the energy is proportional to the number of contacts between non-adjacent (*i.e.* non-consecutive) hydrophobic monomers, capturing the essence of the hydrophobic collapse, which is a key driving force in protein folding. Indeed, hydrophobic amino acids tend to cluster together to avoid exposure to the solvent, reflecting the natural behavior of real hydrophobic residues in aqueous environments. In contrast, polar residues are generally indifferent to being near other residues, and their interactions do not significantly affect the energy calculation.

The primary objective is to find the conformation with the lowest possible energy, representing the protein's native folded state. This involves a computational search through possible conformations to identify the optimal arrangement of hydrophobic and polar residues.

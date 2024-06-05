---
title: Proteins
---

```{tip}
The main references for this part are @finkelstein2002protein and @schlick2010molecular.
```

Proteins are macromolecules composed by amino acids linked by peptide bonds. Let's analyse what these words mean.

# Macromolecules

A macromolecule is a molecule composed by a great number of covalently bonded atoms[^macromolecule]. The most common class of macromolecules is that of biopolymers, which comprise three of the four main macromolecular components of life: proteins, nucleic acids and carbohydrates. Biopolymers are, in turn, a subclass of polymers, which are defined as molecules composed by smaller subunits, the *monomers*, covalently linked together. While there exist polymeric substances with more complex architectures, the main macromolecules of life have a linear (chain) structure that makes it possible to assign a one dimensional *sequence* to each molecule. This sequence is just the list of monomers composing the chain, spelt from one chain end to the other.

```{figure} figures/polymer_sequence.png
:name: fig:polymer-sequence
:align: center
:width: 500px

A cartoon of a polymer composed by 9 monomers (the coloured spheres) connected by covalent bonds (black lines). Panel (A) shows two 2D conformations, panel (B) shows the two 1D sequences built by listing the monomers that make up the chain, starting from either end.
```

[](#fig:polymer-sequence) shows an imaginary short polymer composed by 9 monomers of different nature (coloured differently). In general, there are many different spatial arrangements that the same (bio)polymer can take in solution. By contrast, its sequence is fixed, being given by the list of covalently-linked monomers. However, as shown in the figure, in absence of any convention, the sequence can be read from either end, giving raise to an ambiguity. As we will see, there exist conventions for proteins and nucleic acids that get rid of this ambiguity.

[^macromolecule]: "great number" is a purposedly vague qualifier: there is no strict definition about the number of atoms required for a molecule to be dubbed a macromolecule.

# Amino acids

An amino acid (AA) is a molecule that consists of

* an amino group: A functional group containing nitrogen, written as $-NH_2$. At neutral pH (pH $\approx 7$) this group is protonated, *i.e.* it becomes $-NH_3^+$.
* a carboxyl group: A functional group consisting of a carbon atom double-bonded to an oxygen atom and bonded to a hydroxyl group, written as $-COOH$. At neutral pH (pH $\approx 7$) this group is negatively charged by donating a proton, becoming $-COO^-$. Note that the carbon of this group is often named $C'$ to distinguish it from the carbons of the side chains.
* a side chain (R group): A variable group that differs among amino acids and determines the characteristics and properties of each amino acid. The side chain can be as simple as a hydrogen atom (as in glycine) or more complex like a ring structure (as in tryptophan).
* a central carbon atom ($C^\alpha$) that links together the three foregoing chemical groups plus an additional hydrogen atom.

```{figure} figures/amino_acids.png
:name: fig:amino-acids
:align: center

(A) The structure of a generic amino acid with side chain R. (B) A cartoon showing the spatial difference between the left-handed (L) and right-handed (D) amino-acid enantiomers. (C) Another representation showing that looking down the $H-C^\alpha$ towards the latter, the first letter(s) of the chemical groups read in clockwise order spell out a proper word ($CORN$) for L-amino acids, but a non-existing word ($CONR$) for D-amino acids. In (B) and (C) the chemical groups are coloured as in panel (A).
```

The four bonds of the $C^\alpha$ are arranged in a tetrahedral fashion, which means that all amino acids but glycine[^glycine] are chiral molecules: each AA can, in principle, exist into two distinct forms that are one the mirror image of the other (*i.e.* they are [enantiomers](https://en.wikipedia.org/wiki/Enantiomer)). [](#fig:amino-acids) shows the chemical structure of amino acids in panel (A), and the two enantiomers in panel (B) with two different representations. For reasons that are not yet understood, most proteins found in nature are made of L-amino acids[^D-proteins]. By contrast, D-amino acids, which are synthesised by special enzymes, are rare but present in living beings, being involved in some specific biological processes (see *e.g.* [](doi:10.1007/s00018-010-0571-8) and references therein).

## List of amino acids

Given the generic nature of the side chain $R$, there exist countless different amino acid molecules. However, there are 20 standard amino acids that are commonly found in proteins and are encoded by the genetic code. These 20 amino acids are the building blocks of proteins in all known forms of life. 

```{figure} figures/structures_amino_acids.png
:name: fig:structures_AA
:align: center

The chemical structure of the 20 standard amino acids. As discussed in the text, the AAs are separated into categories according to the properties of their side chains ($R$ groups). Adapted from [here](https://commons.wikimedia.org/wiki/File:Amino_Acids002.svg).
```

[](#fig:structures_AA) shows the chemical structure of the 20 standard amino acids, together with their three- and one-letter abbreviations. In the figure the AAs are grouped in the following classes:

1. Charged amino acids: AAs with side chains that are charged at physiological pH, making them highly hydrophilic. These are:
    * Lysine (Lys, K) and Arginine (Arg, R) - positively charged
    * Aspartic acid (Asp, D) and Glutamic acid (Glu, E) - negatively charged
    * Histidine (His, H) - being [amphoteric](https://en.wikipedia.org/wiki/Amphoterism), can be positively charged, negatively charged or neutral depending on the pH, playing a role in enzyme active sites.
2. Polar unchaged armino acids: AAs with side chains that can form hydrogen bonds, making them hydrophilic but not charged at physiological pH. These are:
    * Serine (Ser, S), Threonine (Thr, T), Asparagine (Asn, N), Glutamine (Gln, Q)
    * Cysteine (Cys, C) - has been traditionally considered a polar (hydrophilic) AA, but this view has been challenged, which is why here it is listed among the hydrophobic AAs. See for instance [](doi:10.1016/S0014-5793(99)01122-9).
3. Aromatic amino acids: AAs with side chains with aromatic rings, contributing to protein structure and function through stacking interactions. These are:
    * Phenylalanine (Phe, F), Tryptophan (Trp, W) - the side chains have a strong hydrophobic character
    * Tyrosine (Tyr, Y) - the presence of the hydroxil group $-OH$ endowes this AA with both hydrophobic and hydrophilic features
4. Hydrophobic amino acids: AAs with side chains that are hydrophobic and likely to be found in the interior of proteins. These are:
    * Glycine (Gly, G), Alanine (Ala, A), Valine (Val, V), Isoleucine (Ile, I), Leucine (Leu, L), Methionine (Met, M), Proline (Pro, P)
    * Cysteine (Cys, C) - has been traditionally considered a polar (hydrophilic) AA, but this view has been challenged, which is why here it is listed among the hydrophobic AAs. See for instance [](doi:10.1016/S0014-5793(99)01122-9).

```{note} Non-standard amino acids
:class: dropdown
In addition to the 20 standard amino acids, there are a few non-standard amino acids that are found in some proteins or are used in specialized biological processes. Two notable ones are:

* Selenocysteine (Sec, U): Sometimes referred to as the 21st amino acid, it is incorporated into proteins by a unique mechanism that involves a specific tRNA and a specific sequence in the mRNA.
* Pyrrolysine (Pyl, O): Sometimes referred to as the 22nd amino acid, it is found in some archaeal and bacterial proteins and is also incorporated by a specific tRNA and sequence in the mRNA.

There are also many other amino acids that are not incorporated into proteins but have important roles in metabolism, such as ornithine and citrulline. Additionally, post-translational modifications can lead to the formation of amino acid derivatives within proteins, such as phosphorylated serine or hydroxyproline.
```

Note that in amino acids, the carbon atoms in the side chains are named systematically based on their position relative to the alpha carbon $C^\alpha$. The naming convention is to use successive Greek letters (and possibly numerals for branched side chains) to denote each carbon atom, where the beta carbon, $C^\beta$, is attached to $C^\alpha$, $C^\gamma$ follows the beta carbon, and so on and so forth.

Here are some examples:

* In alanine the side chain has only one carbon, the beta carbon $C^\beta$.
* In leucine the side chain has four carbons and branches at the gamma carbon, so that the names are $C^\beta, C^\gamma, C^{\delta_1},$ and $C^{\delta_2}$.
* In lysine the side chain has four carbons: $C^\beta, C^\gamma, C^\delta$, and $C^\epsilon$.

[^glycine]: To convince yourself that glycine is not chiral substitute R with H in the rightmost subpanel of [](#fig:amino-acids)(C) and look down the $R-C^\alpha$ bond towards $C^\alpha$: you will see that the resulting view is the same as that of the leftmost panel, which makes the two "enantiomers" identical.
[^D-proteins]: Incorporation of D-amino acids in proteins has been observed to occur only outside of ribosomes. [](doi:10.1007/s00018-010-0571-8) reports some examples.

(peptide-bond)=
# The peptide bond

Two amino acids can be linked together through a *peptide bond*, which is a covalent bond that connects the carboxyl group ($-COOH$) of one AA to the amino group ($-NH_2$) of the other. The reaction through which a peptide bond is formed is called "dehydration synthesis" since the carboxyl group loses a hydroxyl group ($-OH$), and the amino group loses a hydrogen atom ($-H$), which combine to form a water molecule ($H_2O$) that is released in solution. Once an AA has been incorporated into the growing chain and the water molecule has been removed, what remains of the molecule is called an *amino acid residue*, or just *residue*. Note that under physiological conditions this reaction is not spontaneous, and therefore requires catalysis (usually performed by ribosomes).

```{figure} figures/peptide_bond.png
:name: fig:peptide-bond
:align: center

The formation of a dipeptide through a dehydration synthesis: two amino acids are joined together, and a water molecule is subsequently released. The peptide bond is indicated by the arrow.
```

[](#fig:peptide-bond) shows the chemical reaction that takes place when two amino acids are joined together, highlighting the peptide bond that links them. The generic term for a chain of amino acids connected by peptide bonds is *polypeptide*, which is a class of biopolymers that comprises (but it is not strictly equivalent to) proteins. While there is some ambiguity in the definition, a polypeptide is considered to be a protein when it takes a specific three-dimensional structure[^intrinsically_disordered] and a well-defined biological function. Note that by this definition a protein is not necessarily formed by a single polypeptide chain. Indeed, there are many proteins, like hemoglobin, that are made of multiple polypeptide units linked by non-covalent bonds (see [](#sec:quaternary_structure) for a more comprehensive discussion on this matter).

```{important}
The peptide bond has a partial double-bond character arising from resonance: the lone pair of electrons on the nitrogen delocalizes to the carbonyl oxygen, creating a double-bond character between the carbon and nitrogen that restricts rotation around the peptide bond, making it planar and rigid.
```

The linear sequence of amino acids that are covalently linked together by peptide bonds to form a polypeptide chain is called *the primary sequence* of a protein. The primary sequence is the most basic level of protein structure and dictates the specific order in which amino acids are arranged, and it determines the protein's ultimate shape and function through interactions that lead to higher levels of structural organization which we will discuss later on. As a result, any changes or mutations in the primary sequence can significantly impact the protein's overall structure and function.

[^intrinsically_disordered]: Contrary to this definition, [not all proteins have a specific three dimensional structure](https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins), either because they are fully intrinsically-disordered, or because they contain regions that are intrinsically disordered.

(sec:molecular_vibrations)=
# Molecular vibrations

The length of chemical (covalent) bonds is of the order of an angstrom, with $C-H$, $O-H$ and $N-H$ being almost exactly $1 \,\AA$, $C-C$ being $\approx 1.5 \,\AA$, and $C=O$ and the peptide bond $C-N$ being in between ($\approx 1.3\, \AA$).

Regarding covalent bond angles, these depend on the nature of the hybridization. We are primarily concerned with sp$^2$- and sp$^3$-hybridized atoms, which result in planar (approximately $120^\circ$) and tetrahedral (approximately $109.5^\circ$) structures, respectively. In polypeptides, sp$^2$ hybridization is observed in the carbon and nitrogen atoms of the peptide bond, while sp$^3$ hybridization is seen in the alpha carbon, which forms four bonds, and in oxygen and sulfur atoms, which typically form two bonds.

As we will discuss [later](./all_atom.md), vibrations of bond lengths and angles have characteristic frequencies associated to the infrared (IR):

* The bond lengths of interest have typical frequencies that go from $3 \cdot 10^{13}$ Hz for $C-C$ to $9 \cdot 10^{13}$ Hz for $C-H$ (corresponding to 1000 cm$^{-1}$ and 3000 cm$^{-1}$, respectively)
* Bond angles have typical frequencies of $1.8 \cdot 10^{13}$ to $2.7 \cdot 10^{13}$ Hz for $C^\alpha$, $1.5 \cdot 10^{13}$ to $2.4 \cdot 10^{13}$ Hz for $C$, and $2.1 \cdot 10^{13}$ to $3.6 \cdot 10^{13}$ Hz for $N$ (500 to 1200 cm$^{-1}$ for the total range).

The values of these typical frequencies can be compared to the thermal energy, $k_B T$, by using Planck's relation, $U = h \nu$, which yields

$$
\nu_{\rm th} = \frac{k_B T}{h} \approx 6.2 \cdot 10^{12} \, {\rm Hz}.
$$

Noting that this value is roughly half of the lowest vibrational or bending frequency reported above, we can conclude that thermal energy does not significantly impact bond stretching and angle bending vibrations[^bond_angles]. As a result, covalent bond and angle vibrations contribute very little to the conformational flexibility of polypeptides, at least under normal conditions.

[^bond_angles]: Temperature-induced fluctuations of bond angles are somewhat larger than those of bond lengths, and therefore cannot be disregarded completely. However the effect is rather minor, as typical amplitudes ar of the order of 5$^\circ$ (@finkelstein2002protein).

By contrast, the typical frequencies of rotations around single bonds are generally much lower than those of bond stretching or bond angle bending vibrations. In fact, some of these frequencies are in the range that is accessible by thermal energy at room temperature, endowing polypeptides with a rotational flexibility that is essential for their conformational dynamics, allowing them to adopt various functional states.

```{figure} figures/dihedral.png
:name: fig:dihedral
:align: center
:width: 400px

The dihedral angle $\varphi$ is defined as the angle formed by the planes determined by the $ABC$ and $BCD$ atoms, which are connected by covalent bonds represented by blue lines. Note that in this picture only one of the two possible choices for the  angle, *i.e.* the dihedral angle $\varphi$, is shown explicitly (see text for details).
```

The rotation angle around a bond is called torsional or *dihedral* angle, which, as shown in [](#fig:dihedral), is the angle formed by two intersecting planes determined by the positions of four atoms. Given the four atoms in figure, $ABCD$, let $\overrightarrow{XY}$ be the vector connecting two covalently linked atoms $X$ and $Y$. Given the ambiguity of defining the normal to a plane, the definition of the torsional angle requires a convention. If we choose to define the normal vector to the $ABC$ plane as $\vec{n}_1 = \overrightarrow{AB} \times \overrightarrow{BC}$, the two normal vectors to the $BCD$ plane are

$$
\begin{aligned}
\vec{n}_2 & = \overrightarrow{DC} \times \overrightarrow{CB}\\
\vec{n}'_2 & = \overrightarrow{BC} \times \overrightarrow{CD} = -\vec{n}_2.
\end{aligned}
$$

We can now define two torsional angles, $\varphi$, $\hat{\varphi}$, as the angles between the normal vectors, *viz.*

$$
\begin{aligned}
\cos \varphi &= \frac{\vec{n}_1 \cdot \vec{n}_2}{\lVert \vec{n}_1 \rVert \lVert \vec{n}_2 \rVert}\\
\cos \hat{\varphi} &= \frac{\vec{n}_1 \cdot \vec{n}'_2}{\lVert \vec{n}_1 \rVert \lVert \vec{n}'_2 \rVert} = -\cos \varphi.
\end{aligned}
$$

The two angles are connected by the relation $\varphi + \hat{\varphi} = \pi$. Following @schlick2010molecular, I will call $\varphi$ the dihedral angle, and $\hat{\varphi}$ the torsional angle, although the two terms are often used interchangeably.

As noted [above](#peptide-bond), the peptide bond, whose associated dihedral angle is called $\omega$, has a partial double-bond character that restricts rotations around it. The peptide bond is considered to be "planar", *i.e.* that the dihedral angle takes values $\omega = \pi$ or $\omega = 0$, with the latter being somewhat less common. Deviations from these values are considered to be rare, but this view has been challenged (see *e.g.* [](doi:10.1073/pnas.1107115108)).

By contrast, rotations around bonds that connect sp$^2$- and sp$^3$-hybridized atoms are associated to energy barriers that are of the order of $k_BT$ and therefore are the main contributors to the flexibility of the macromolecules. In the main chains of peptides the dihedral angles involved in these rotations are those associated to the $N - C^\alpha$ and $C^\alpha - C$ bonds, which are called $\phi$ and $\psi$. 

## Ramachandran Plots

The small energy barriers associated to rotations around $\phi$ and $\psi$ only take into account the atom themselves, not the rest of the molecule to which they are attached. Steric hindrances, resulting from the repulsive interactions between atoms and groups attached to the $C$, $C^\alpha$ and $N$ atoms, disfavor or even prohibit certain combinations of $(\phi, \psi)$ values.

A graphical representation of the allowed and disallowed dihedral angles of amino acid residues in protein structures is the Ramachandran plot, introduced for the first time in [](doi:10.1016/S0022-2836(63)80023-6). By plotting $\phi$ and $\psi$ on the $x$ and $y$ axes, respectively, the plot reveals regions where the angles are sterically favorable, corresponding to common secondary structures and motifs like $\alpha$-helices and $\beta$-sheets which will be introduced below. This visualization is crucial for validating protein structures, as it highlights conformational possibilities and identifies potential structural anomalies (*i.e.* AA conformations that lie in disallowed regions of the plot).

```{warning}
In many papers and books on the subject of protein structure, there is an important fact that is often not stressed or even omitted (perhaps because it should be obvious to people who, unlike us, know some biochemistry): in the original explanation of the Ramachandran plot, the steric clashes that contribute to excluding some $(\phi, \psi)$ combinations are those involving at least third neighbours (often called 1-4 interactions).
```

Given a protein structure, it is possible to extract and plot the $\phi$ and $\psi$ values obtained for each residue on the same figure. However, different amino acids display different flexibilities depending on their associated side chains. Therefore, it is common to produce multiple Ramachandran plots, each serving specific purposes and providing insights into various aspects of protein structure and conformation. The most common ones are:

* **Glycine Ramachandran Plot:** glycine residues have more conformational freedom due to the absence of a side chain, and therefore a wider range of allowed $\phi$ and $\psi$ angles compared to the general case, reflecting the famous enhanced flexibility of this AA.
* **Proline Ramachandran Plot:** proline residues have restricted $\phi$ and $\psi$ angles due to the cyclic nature of its side chain, resulting in a limited range of conformations, highlighting the unique structural constraints of this AA.
* **Pre-Proline Ramachandran Plot:** residues that precede proline in the primary sequence often exhibit distinct conformational preferences given by the steric influence of proline.
* **Ile-Val Ramachandran Plot:** the branched carbons of isoleucine (Ile) and valine (Val) give them a distinct shape of disallowed $\phi$-$\psi$ regions.
* **General Ramachandran Plot:** $\phi$ and $\psi$ angles for all residues that are not part of one of the foregoing categories.

```{figure}
:name: fig:ramachandran_plot
:align: center

![The Ramachandran plot of glycines.](figures/1a3n_Glycine.png)
![The Ramachandran plot of residues that precede a proline ("pre-proline").](figures/1a3n_Pre-proline.png)
![The Ramachandran plot of prolines.](figures/1a3n_Proline.png)
![The Ramachandran plot of of leucines and valines.](figures/1a3n_Ile-Val.png)
![The Ramachandran plot of all amino acids that are not part of any of the above classes (also known as "general plot").](figures/1a3n_General.png)
![The Ramachandran plot of all amino acids.](figures/1a3n_All.png)

Ramachandran plots of amino acids that compose the human hemoglobin protein [1A3N](https://www.rcsb.org/structure/1A3N) (points), plotted on background  contour lines generated by analysing the Top8000 PDB data set. I have used [this software](https://github.com/Joseph-Ellaway/Ramachandran_Plotter) to make the plots.
```

[](#fig:ramachandran_plot) shows examples of the plot types described above, with the empty (non-coloured) "space" being associated to angle combinations that are sterically forbidden. The points in the panels refer to the $(\phi, \psi)$ combinations of the residues of the human hemoglobin protein, while the background contour plots have been calculated by extracting the values of $\phi$ and $\psi$ from 8000 reference protein structures listed in a database[^ramachandran_plot_reference].

A comparison between [](#fig:ramachandran_plot-a) and the other plots shows that glycine has a much greater conformational flexibility than any other AA. By contrast, [pre-prolines](#fig:ramachandran_plot-b) and [prolines](#fig:ramachandran_plot-c) (and, to a smaller extent, [leucines and valines](#fig:ramachandran_plot-d)) are much more constrained than an [average AA](#fig:ramachandran_plot-e).

[^ramachandran_plot_reference]: While it is not explicitly stated in the repository of [the software](https://github.com/Joseph-Ellaway/Ramachandran_Plotter) I have used to make the figure, I believe that the reference Top8000 database used is [this one](http://kinemage.biochem.duke.edu/research/top8000/). The updated version of the same database can instead be downloaded [here](https://zenodo.org/records/5777651).

## Rotamers

Almost all amino acids can adopt different orientations of side chains around single bonds. Exceptions are 

* glycine, which has no rotational freedom in terms of side chain orientation,
* alanine, whose side chain is a methyl group, $-CH3$, which is small and symmetrical, offering no significant variation in rotameric states,
* proline, whose side chain forms a ring by bonding back to the backbone nitrogen, severely restricting the rotation around its side chain.

The distinct conformations are called *rotamers* and arise due to the rotation around the bonds connecting the side chains to the main backbone of the polypeptide chain. The dihedral angles connected to the rotamers are referred to as $\chi$ angles, which describe the rotations around the bonds within the side chains of amino acids. Each rotatable bond in a side chain is associated with a specific $\chi$ angle, labeled sequentially from the backbone outward. For instance, $\chi_1$ is the dihedral angle around the bond between the alpha carbon and the beta carbon, $\chi_2$ is the dihedral angle around the bond between the beta carbon and the gamma carbon, and so on for longer side chains.

Each amino acid side chain can adopt multiple rotameric states, influenced by steric interactions, hydrogen bonding, and other intramolecular forces. The different rotamers contribute to the overall flexibility and diversity of protein structures, allowing proteins to achieve their functional conformations and interact effectively with other molecules.



# Secondary structure

# Tertiary structure

(sec:quaternary_structure)=
# Quaternary structure

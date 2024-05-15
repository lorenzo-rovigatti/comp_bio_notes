---
title: Proteins
---

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
* a carboxyl group: A functional group consisting of a carbon atom double-bonded to an oxygen atom and bonded to a hydroxyl group, written as $-COOH$. At neutral pH (pH $\approx 7$) this group is negatively charged by donating a proton, becoming $-COO^-$.
* a side chain (R group): A variable group that differs among amino acids and determines the characteristics and properties of each amino acid. The side chain can be as simple as a hydrogen atom (as in glycine) or more complex like a ring structure (as in tryptophan).
* a central carbon atom ($C^\alpha$) that links together the three foregoing chemical groups plus an additional hydrogen atom.

```{figure} figures/amino_acids.png
:name: fig:amino-acids
:align: center

(A) The structure of a generic amino acid with side chain R. (B) A cartoon showing the spatial difference between the left-handed (L) and right-handed (D) amino-acid enantiomers. (C) Another representation showing that looking down the $H-C^\alpha$ towards the latter, the first letter(s) of the chemical groups read in clockwise order spell out a proper word ($CORN$) for L-amino acids, but a non-existing word ($CONR$) for D-amino acids. In (B) and (C) the chemical groups are coloured as in panel (A).
```

The four bonds of the $C^\alpha$ are arranged in a tetrahedral fashion, which means that all amino acids but glycine[^glycine] are chiral molecules: each AA can, in principle, exist into two distinct forms that are one the mirror image of the other (*i.e.* they are [enantiomers](https://en.wikipedia.org/wiki/Enantiomer)). [](#fig:amino-acids) shows the chemical structure of amino acids in panel (A), and the two enantiomers in panel (B) with two different representations. For reasons that are not yet understood, most proteins found in nature are made of L-amino acids[^D-proteins], while D-amino acids are rare but present in living beings, being involved in some specific biological processes (see *e.g.* [](doi:10.1007/s00018-010-0571-8) and references therein).

## List of amino acids

Given the generic nature of the side chain $R$, there exist countless different amino acid molecules. However, there are 20 standard amino acids that are commonly found in proteins and are encoded by the genetic code. These 20 amino acids are the building blocks of proteins in all known forms of life. 

```{figure} figures/structures_amino_acids.png
:name: fig:structures_AA
:align: center

The chemical structure of the 20 standard amino acids. As discussed in the text, the AAs are separated into categories according to the properties of their side chains ($R$ groups). Adapted from [here](https://commons.wikimedia.org/wiki/File:Amino_Acids002.svg).
```

[](#fig:structures_AA) shows the chemical structure of the 20 standard amino acids, together with their three- and one-letter abbreviations. In the figure the AA are grouped in the following classes:

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

[^glycine]: To convince yourself that glycine is not chiral substitute R with H in the rightmost subpanel of [](#fig:amino-acids)(C) and look down the $R-C^\alpha$ bond towards $C^\alpha$: you will see that the resulting view is the same as that of the leftmost panel, which makes the two "enantiomers" identical.
[^D-proteins]: Incorporation of D-amino acids in proteins has been observed to occur only outside of ribosomes. [](doi:10.1007/s00018-010-0571-8) reports some examples.

# The peptide bond

Two amino acids can be linked together through a *peptide bond*, which is a covalent bond that connects the carboxyl group ($-COOH$) of one AA to the amino group ($-NH_2$) of the other. The reaction through which a peptide bond is formed is called "dehydration synthesis" since the carboxyl group loses a hydroxyl group ($-OH$), and the amino group loses a hydrogen atom ($-H$), which combine to form a water molecule ($H_2O$) that is released in solution. Once an AA has been incorporated into the growing chain and the water molecule has been removed, what remains of the molecule is called an *amino acid residue*, or just *residue*. Note that under physiological conditions this reaction is not spontaneous, and therefore requires catalysis (usually performed by ribosomes).

```{figure} figures/peptide_bond.png
:name: fig:peptide-bond
:align: center

The formation of a dipeptide through a dehydration synthesis: two amino acids are joined together, and a water molecule is subsequently released. The peptide bond is indicated by the arrow.
```

[](#fig:peptide-bond) shows the chemical reaction that takes place when two amino acids are joined together, highlighting the peptide bond that links them. The generic term for a chain of amino acids connected by peptide bonds is *polypeptide*, which is a class of biopolymers that comprises (but it is not strictly equivalent to) proteins. While there is some ambiguity in the definition, a polypeptide is considered to be a protein when it takes a specific three-dimensional structure[^intrinsically_disordered] and a well-defined biological function.



```{important}
The peptide bond has a partial double-bond character arising from resonance: the lone pair of electrons on the nitrogen delocalizes to the carbonyl oxygen, creating a double-bond character between the carbon and nitrogen. This resonance restricts rotation around the peptide bond, making it planar and rigid.
```

[^intrinsically_disordered]: Contrary to this definition, [not all proteins have a specific three dimensional structure](https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins), either because they are fully intrinsically-disordered, or because they contain regions that are intrinsically disordered.

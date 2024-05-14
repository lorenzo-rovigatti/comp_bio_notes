---
title: Basic notions
---

# Introduction to the course

As with most computational courses, this course is supposed to have a practical side that should not be overlooked. However, as some of you may have noticed, there will be frontal lessons only. While this may seem contradictory (and in some sense it is), it also means that you are strongly advised to practice on your (or someone else's) computer what you'll hear about (and be shown) during the lectures. Moreover, I will also set up some hands-on (bring-your-own-laptop) lectures to guide you through the most important technical hurdles we will be encountering.

```{tip}
Most of the acronyms that are scattered throughout the text show tooltips when hovered with the mouse cursor, like HPC.
```

## Prerequisites

1. This is a computational course, and as such it requires some proficiency with (or at least having the right attitude towards) computers. The most important skill you should have is using aterminal, since most of the computations (running simulations, analysing results, *etc.*) will be launched from there. Linux and macOS come with pre-installed terminals, while Windows does not[^windows_terminal]. However, it is possible to install a very handy ["Windows Subsystem for Linux"](https://learn.microsoft.com/en-us/windows/wsl/) that makes it possible to have a good Linux-like experience within Windows. Install instructions can be found [here](https://learn.microsoft.com/en-us/windows/wsl/install).
2. Knowing to code is important, for several reasons. First of all, implementing algorithms and techniques introduced in class is the best way of understanding how they work and what are their advantages and disadvantages, even though most of research-grade results during one's career will be obtained with production-ready HPC codes written by others. Secondly, reading other people's code can come in very handy (i) to understand what it does and (ii) to extend it to suit your needs. Finally, in computational physics it is common to need to perform custom analyses for which no libraries or codes are available, leaving no other option than writing your own script or program. For the sake of the course, the programming language you are more familiar with is not important, although Python is the *de facto* standard language used to write analysis scripts, while C and C++ (and more rarely FORTRAN) are used to write performance-critical software (or part thereof).
3. By definition, life is an out-of-equilibrium process, since it continuously uses up energy. However, many biophysical processes are in (or close to) equilibrium. Therefore, they can be analysed and understood using the language of thermodynamics and statistical mechanics. You should be familiar with concepts such as free energy, entropy, ensembles, partition function, Boltzmann distribution, as well as with the mathematical tools used to work with these quantities.
4. *Any* knowledge of biology is welcome (and will be useful to its holder). However, I will introduce most of what we need at the beginning, and then some more as the need arises.

[^windows_terminal]: Or at least it does not come with a Unix-like (POSIX) terminal, which is what we will be using.

# Introduction to computational physics

# Introduction to biophysics

## The central dogma of biology

The single most important fact about biology is the so-called "central dogma of biology", which is a framework that describes the flow of genetic information in biological systems. It was formulated by Francis Crick[^dogma] and presented to the public in a famous lecture given on September 19th in 1957 [](doi:10.1371/journal.pbio.2003243). The first formulation of the central dogma is shown in [](#central-dogma), where arrows indicate the *flow of information*: what is used to build what. Note that DNA $\to$ protein and RNA $\to$ DNA information transfers were highly speculative at that time (and in fact DNA $\to$ protein transfers have not been observed *in vivo*).

```{figure} figures/central_dogma.png
:name: central-dogma
:align: center

The first outline of the central dogma, from un unpublished note made by Crick in 1956. Credit: [Wellcome Library](https://wellcomecollection.org/works/xmscu3g4), London. Taken from [](doi:10.1371/journal.pbio.2003243).
```

In later formulations Crick expanded and clarified his view on this matter, writing that 

> the central dogma could be stated in the form "once (sequential) information has passed into protein it cannot get out again".
>
> [...]
>
> the central dogma was a negative statement, saying that transfers from protein did not exist.
> 
> -- [](doi:10.1038/227561a0)

**Nota Bene:** this version of the central dogma is not only the original one, but also the only one that stood the test of time: popularised versions such as "DNA $\to$ RNA $\to$ protein" or "DNA makes RNA, RNA makes protein" are oversimplified and ambiguous statements[^ambiguous_dogma].

```{attention}
Try asking ChatGPT what is the central dogma of biology: Do you see anything wrong with the answer?
```

[^dogma]: Apparently Crick was not familiar with the meaning of "dogma", which is why he used that word instead of other options that would have perhaps avoided later misunderstandings (see *e.g.* [](doi:10.1371/journal.pbio.2003243)).
[^ambiguous_dogma]: Oversimplified because they leave out important information transfers such as RNA $\to$ RNA and RNA $\to$ DNA, but also ambiguous since the meaning of arrows and "makes" are not spelt out.

## Proteins

Proteins are macromolecules composed by amino acids linked by peptide bonds. Let's analyse what these words mean.

### Macromolecules

A macromolecule is a molecule composed by a great number of covalently bonded atoms[^macromolecule]. The most common class of macromolecules is that of biopolymers, which comprise three of the four main macromolecular components of life: proteins, nucleic acids and carbohydrates. Biopolymers are, in turn, a subclass of polymers, which are defined as molecules composed by smaller subunits, the *monomers*, covalently linked together. While there exist polymeric substances with more complex architectures, the main macromolecules of life have a linear (chain) structure that makes it possible to assign a one dimensional *sequence* to each molecule. This sequence is just the list of monomers composing the chain, spelt from one chain end to the other.

```{figure} figures/polymer_sequence.png
:name: fig:polymer-sequence
:align: center
:width: 500px

A cartoon of a polymer composed by 9 monomers (the coloured spheres) connected by covalent bonds (black lines). Panel (A) shows two 2D conformations, panel (B) shows the two 1D sequences built by listing the monomers that make up the chain, starting from either end.
```

[](#fig:polymer-sequence) shows an imaginary short polymer composed by 9 monomers of different nature (coloured differently). In general, there are many different spatial arrangements that the same (bio)polymer can take in solution. By contrast, its sequence is fixed, being given by the list of covalently-linked monomers. However, as shown in the figure, in absence of any convention, the sequence can be read from either end, giving raise to an ambiguity. As we will see, there exist conventions for proteins and nucleic acids that get rid of this ambiguity.

[^macromolecule]: "great number" is a purposedly vague qualifier: there is no strict definition about the number of atoms required for a molecule to be dubbed a macromolecule.

### Amino acids

An amino acid (AA) is a molecule that consists of

* an amino group: A functional group containing nitrogen, written as -NH$_2$. At neutral pH (pH $\approx 7$) this group is protonated, *i.e.* it becomes -NH$_3^+$.
* a carboxyl group: A functional group consisting of a carbon atom double-bonded to an oxygen atom and bonded to a hydroxyl group, written as -COOH. At neutral pH (pH $\approx 7$) this group is negatively charged by donating a proton, becoming -COO$^-$.
* a side chain (R group): A variable group that differs among amino acids and determines the characteristics and properties of each amino acid. The side chain can be as simple as a hydrogen atom (as in glycine) or more complex like a ring structure (as in tryptophan).
* a central carbon atom (C$^\alpha$) that links together the three foregoing chemical groups plus an additional hydrogen atom.

```{figure} figures/amino_acids.png
:name: fig:amino-acids
:align: center

(A) The structure of a generic amino acid with side chain R. (B) A cartoon showing the spatial difference between the left-handed (L) and right-handed (D) amino-acid enantiomers. (C) Another representation showing that looking down the H-C$^\alpha$ towards the latter, the first letter(s) of the chemical groups read in clockwise order spell out a proper word (CORN) for L-amino acids, but a non-existing word (CONR) for D-amino acids. In (B) and (C) the chemical groups are coloured as in panel (A).
```

The four bonds of the C$^\alpha$ are arranged in a tetrahedral fashion, which means that all amino acids but glycine[^glycine] are chiral molecules: each AA can, in principle, exist into two distinct forms that are one the mirror image of the other (*i.e.* they are [enantiomers](https://en.wikipedia.org/wiki/Enantiomer)). [](#fig:amino-acids) shows the chemical structure of amino acids in panel (A), and the two enantiomers in panel (B) with two different representations. For reasons that are not yet understood, most proteins found in nature are made of L-amino acids[^D-proteins], while D-amino acids can be found in other biological processes (see *e.g.* [](doi:10.1007/s00018-010-0571-8) and references therein).

[^glycine]: To convince yourself that glycine is not chiral substitute R with H in the rightmost subpanel of [](#fig:amino-acids)(C) and look down the R-C$^\alpha$ bond towards C$^\alpha$: you will see that the resulting view is the same as that of the leftmost panel, which makes the two "enantiomers" identical.
[^D-proteins]: Incorporation of D-amino acids in proteins has been observed to occur only outside of ribosomes. [](doi:10.1007/s00018-010-0571-8) reports some examples.

### The peptide bond

```{figure} figures/codon-flower.svg
:name: fig:codon-flower
:align: center

The meaning of each sequence made of three nucleotides (also known as a *codon* or *codon triplet*), as read by the ribosome. Each codon should be spelt starting from the inner circle and going outwards. Out of 64 combinations, three are recognised as "stop signals", while the rest encode amino acids, with the exception of the AUG codon which, in addition to represent methionine, is also typically used to initiate protein syntesis. Most amino acids (excluding methionine and tryptophan) are encoded by more than a single codon.
```

## Nucleic acids

## Towards larger length-scales

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
4. *Any* knowledge of biology and biochemistry is welcome (and will be useful to its holder). However, I will introduce most of what we need at the beginning, and then some more as the need arises.

[^windows_terminal]: Or at least it does not come with a Unix-like (POSIX) terminal, which is what we will be using.

# About these notes

* Internal links are styled [like this](#sec:intro_comp_phys). Hover a link to see a preview.
* Wikipedia links are styled [like this](https://en.wikipedia.org/wiki/Biophysics). Hover a wikipedia link to see a preview.
* External links are styled [like this](https://www.youtube.com/watch?v=dQw4w9WgXcQ).
* Hover over a footnote reference to show the associated text[^footnote].
* Hover over a reference to show the associated citation and a clickable DOI (see *e.g.* [](doi:10.1126/science.177.4047.393)).

[^footnote]: like this.

(sec:intro_comp_phys)=
# Introduction to computational physics

(sec:intro_biophys)=
# Introduction to biophysics

## Some important numbers and quantities

| Quantity | Value |
| :--- | :---: |
| $k_B T$ at 298 K| $4.11 \times 10^{-21}$ J|
| " | $4.11 \, {\rm pn \cdot nm}$|
| RT at 298 K | $0.59 \, {\rm kcal / mol}$|
| " | $2.48 \, {\rm kJ / mol}$|
| Average nucleotide mass|$\approx 330 \to 300$ Da|
| Average amino acid mass|$\approx 110 \to 100$ Da|

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

## The genetic code and the information flow

```{figure} figures/codon-flower.svg
:name: fig:codon-flower
:align: center

The meaning of each sequence made of three nucleotides (also known as a *codon* or *codon triplet*), as read by the ribosome. Each codon should be spelt starting from the inner circle and going outwards. The outer circle shows 
```

[](#fig:codon-flower) shows the meaning of each codon (a sequence of three nucleotides). Since each nucleotide in a codon can have one out of 4 possible values (A, C, G or T), there exist $4^3 = 64$ possible codons. Out of these 64 combinations, three are recognised as "stop signals", while the rest encode amino acids, with the exception of the AUG codon which, in addition to represent methionine, is also typically used to initiate protein syntesis. Most amino acids (excluding methionine and tryptophan) are encoded by more than a single codon. The outer circle of [](#fig:codon-flower) contains the list of standard amino acids, together with their three- and one-letter abbreviations, which are shown in bold and after each AA name, respectively.

## The experimental characterisation of the structure of biomolecules

### Three-dimensional structure

The 3D structure of biomacromolecules such as proteins and nucleic acids can be experimentally determined using several techniques. The most commonly used methods are:

1. **X-ray Crystallography**:
   - **Procedure**: The biomolecule is crystallized, and then X-rays are directed at the crystal. The X-rays diffract as they pass through the crystal, creating a diffraction pattern that can be recorded. The diffraction pattern is analyzed to determine the electron density map, which is used to model the atomic structure of the biomolecule.
   - **Advantages**: High resolution, can be used for large biomolecules.
   - **Limitations**: Requires high-quality crystals, which can be difficult to obtain for some biomolecules (*e.g.* molecules that are not soluble in water).

2. **Nuclear Magnetic Resonance (NMR) Spectroscopy**:
   - **Procedure**: The biomolecule is dissolved in solution, and NMR spectra are obtained by applying a magnetic field and measuring the resulting interactions between nuclear spins. The spectra provide information about the distances and angles between atoms, which are used to calculate the 3D structure.
   - **Advantages**: Provides information about the dynamics and flexibility of biomolecules in solution.
   - **Limitations**: Generally limited to smaller biomolecules (up to about 50 kDa), although recent advances are pushing these limits.

3. **Cryo-Electron Microscopy (cryo-EM)**:
   - **Procedure**: The biomolecule is rapidly frozen in a thin layer of ice, and then imaged using an electron microscope. Thousands of 2D images are collected and computationally combined to reconstruct a 3D model of the biomolecule.
   - **Advantages**: Does not require crystallization, suitable for large and complex biomolecules, including those that are difficult to crystallize.
   - **Limitations**: Typically lower resolution than X-ray crystallography, although recent advances have significantly improved the resolution.

4. **Small-Angle X-ray Scattering (SAXS)**:
   - **Procedure**: The biomolecule is dissolved in solution, and X-rays are scattered at small angles by the sample. The scattering pattern provides low-resolution information about the overall shape and size of the biomolecule.
   - **Advantages**: Can study biomolecules in near-physiological conditions, provides information about conformational changes and flexibility.
   - **Limitations**: Lower resolution compared to X-ray crystallography and cryo-EM, providing only general shape information.

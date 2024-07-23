---
title: Basic notions
---

# Introduction to the course

As with most computational courses, this course is supposed to have a practical side that should not be overlooked. However, as some of you may have noticed, there will be frontal lessons only. While this may seem contradictory (and in some sense it is), it also means that you are strongly advised to practice on your (or someone else's) computer what you'll hear about (and be shown) during the lectures. Moreover, I will also set up some hands-on (bring-your-own-laptop) lectures to guide you through the most important technical hurdles we will be encountering.

# About these notes

* Internal links are styled [like this](#sec:intro_comp_phys). Hover a link to see a preview.
* Wikipedia links are styled [like this](https://en.wikipedia.org/wiki/Biophysics). Hover a wikipedia link to see a preview.
* External links are styled [like this](https://www.youtube.com/watch?v=dQw4w9WgXcQ).
* Hover over a footnote reference to show the associated text[^footnote].
* Hover over a reference to show the associated citation and a clickable DOI (see *e.g.* [](doi:10.1126/science.177.4047.393)).
* Most of the acronyms that are scattered throughout the text show tooltips when hovered with the mouse cursor, like HPC.

:::{tip} A box
There are some boxes scattered throughout the text. Their colour weakly correlates with the content: pay particular attention to the yellow and red ones!
:::

[^footnote]: like this.

## Prerequisites

1. This is a computational course, and as such it requires some proficiency with (or at least having the right attitude towards) computers. The most important skill you will need is to use a terminal, since most of the computations (running simulations, analysing results, *etc.*) will be launched from there. Linux and macOS come with pre-installed terminals, while Windows does not[^windows_terminal]. However, it is possible to install a very handy ["Windows Subsystem for Linux"](https://learn.microsoft.com/en-us/windows/wsl/) that makes it possible to have a good Linux-like experience within Windows. Install instructions can be found [here](https://learn.microsoft.com/en-us/windows/wsl/install).
2. Knowing to code is important, for several reasons. First of all, implementing algorithms and techniques introduced in class is the best way of understanding how they work and what are their advantages and disadvantages, even though most of research-grade results during one's career will be obtained with production-ready HPC codes written by others. Secondly, reading other people's code can come in very handy (i) to understand what it does and (ii) to extend it to suit your needs. Finally, in computational physics it is common to need to perform custom analyses for which no libraries or codes are available, leaving no other option than writing your own script or program. For the sake of the course, the programming language you are more familiar with is not important, although Python is the *de facto* standard language used to write analysis scripts, while C and C++ (and more rarely FORTRAN) are used to write performance-critical software (or part thereof).
3. By definition, life is an out-of-equilibrium process, since it continuously uses up energy. However, many biophysical processes are in (or close to) equilibrium. Therefore, they can be analysed and understood using the language of thermodynamics and statistical mechanics. You should be familiar with concepts such as free energy, entropy, ensembles, partition function, Boltzmann distribution, as well as with the mathematical tools used to work with these quantities.
4. *Any* knowledge of biology and biochemistry is welcome (and will be useful to its holder). However, I will introduce most of what we need at the beginning, and then some more as the need arises.

[^windows_terminal]: Or at least it does not come with a Unix-like (POSIX) terminal, which is what we will be using.

## Evaluation

Your final grade will be based on:

* the final group project ($40\%$)
* the part of the individual oral exam devoted to discussing the final project ($30\%$)
* the part of the individual oral exam where you will be asked questions about the lectures OR the homework assignements carried out during the course ($30\%$)

About the last bullet point, during the course you will be asked to complete three homework assignments. Those are not mandatory, but carrying them out satisfactorily will allow you to skip the questions about the course program during the oral exam.

### Homework assignments

More or less every couple of weeks I will ask you to do a homework assignment in which you will have to write a code that does something "useful", connected to what explained in the classroom. These assignments are not mandatory, but if you carry out enough of them in a satisfactorily manner, you will be able to "skip" part of the oral exam where you would be asked questions about the frontal lessons.

The output of each assignment should be:

1. Your code, written in any (sensible) language (*e.g.* C, C++, Python, Fortran, Java, but other languages may be also fine, just let me know beforehand).
2. A short (but not too short!) document containing your code documentation, and some results obtained with the code.

The code documentation should be clear and concise, but also list features, limitations, and known bugs of your code. In addition, it should present a comprehensive guide about how to you use the code (*e.g.* how to call it, how to specify the input, how to read the output, *etc.*). I will provide examples of "good" documentation during the course.

(sec:intro_biophys)=
# Introduction to biophysics

:::{tip} Source
Some parts of this section have been adapted from [here](https://bio.libretexts.org/Bookshelves/Computational_Biology/Book:_Computational_Biology_-_Genomes_Networks_and_Evolution_(Kellis_et_al.)/01:_Introduction_to_the_Course/1.04:_Crash_Course_in_Molecular_Biology).
:::

## Some definitions

Procaryotes
: Unicellular organisms that lack a true nucleus and membrane-bound organelles. Their genetic material is located in a nucleoid, and they typically have a single, circular chromosome along with plasmids. Prokaryotic cells have a rather small size ($0.1-5.0 \, {\rm \mu m}$) and reproduce through binary fission[^binary_fission]. Most have a rigid cell wall made of peptidoglycan (in bacteria). Examples of prokaryotes include bacteria and archaea.

Eukaryotes
: Organisms whose cells contain a true nucleus enclosed by a nuclear membrane and various membrane-bound organelles such as mitochondria, endoplasmic reticulum, and Golgi apparatus. Their DNA is linear and organized into chromosomes within the nucleus. Eukaryotic cells are generally larger than prokaryotes ($10-100 \, {\rm \mu m}$) and reproduce through mitosis for somatic cells and meiosis for gametes. While some eukaryotes have cell walls (*e.g.*, plants and fungi), others do not. Examples of eukaryotes include plants, animals, fungi, and protists.

[^binary_fission]: For what concerns us, this is a simpler version of mitosis.

## DNA

The DNA molecule stores the genetic information of an organism. DNA contains regions called genes, which encode for proteins to be produced. Other regions of the DNA contain regulatory elements, which partially influence the level of expression of each gene. Within the genetic code of DNA lies both the data about the proteins that need to be encoded, and the control circuitry, in the form of regulatory motifs.

As we will see in depth, in cells DNA is usually found in a "double-stranded" helical form, where each strand is a long chain of repeating units, called nucleotides, of just four types: A(adenine), C(cytosine), T (thymine), and G (guanine). The list of nucleotides is called *sequence* or *primary structure*. In the double strand, A pairs with T and G with C, with the A-T pairing being weaker than the C-G pairing[^GC_pairing].

The two DNA strands in the double helix are complementary, meaning that if there is an A on one strand, it will be bonded to a T on the other, and if there is a C on one strand, it will be bonded to a G on the other. The DNA strands have a chemical directionality, or polarity, with the convention being to list a strand's sequence with the direction along which enzymes[^enzyme] called polymerases[^polymerase] synthesise DNA and RNA, called the 5' to 3' direction. With this in mind, we can say that that the DNA strands are anti-parallel, as the 5' end of one strand is adjacent to the 3' end of the other. As a result, DNA can be read both in the 3' to 5' direction and the 5' to 3' direction, and genes and other functional elements can be found in each.

```{figure} figures/chromatin.svg
:name: fig:chromatin
:align: center

The multi-scale organisation of DNA in a eukaryotic cell. Adapted from a figure [made by Phrood, via Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Chromatine.svg).
```

Base pairing between nucleotides of DNA constitutes its *secondary structure*. In addition to DNA's secondary structure, there are several extra levels of structure that allow biological DNA to be tightly compacted and influence gene expression, as schematically shown in [](#fig:chromatin). The tertiary structure describes the twist in the DNA ladder that forms a helical shape. In the quaternary structure, DNA is tightly wound around small proteins called histones. These DNA-histone complexes are further wound into tighter structures seen in chromatin.

Before DNA can be replicated or transcribed into RNA, the chromatin structure must be locally "unpacked". Thus, gene expression may be regulated by modifications to the chromatin structure, which make it easier or harder for the DNA to be unpacked. This regulation of gene expression via chromatin modification is an example of [epigenetics](https://en.wikipedia.org/wiki/Epigenetics).

The structure of DNA allows the strands to be easily separated for the purpose of DNA replication, as well as transcription, translation, recombination, and DNA repair, among others. This was noted by [Watson and Crick](doi:10.1038/171737a0) as "It has not escaped our notice that the specific pairing that we have postulated immediately suggests a possible copying mechanism for the genetic material". In the replication of DNA, the two complementary strands are separated, and each of the strands are used as templates for the construction of a new strand.

DNA polymerases attach to each of the strands at the origin of replication, reading each existing strand from the 3' to 5' direction and placing down complementary bases such that the new strand grows in the 5' to 3' direction. Because polymerases build the chains from 5' to 3', one strand (the leading strand) can be copied continuously, while the other (the lagging strand) grows in pieces which are later glued together by another enzyme, DNA ligase. The end result is 2 double-stranded pieces of DNA, where each is composed of one old strand, and one new strand; for this reason, DNA replication is a semiconservative process.

Many organisms have their DNA split into several chromosomes. Each chromosome contains two strands of DNA, which are complementary to each other but are read in opposite directions. Genes can occur on either strand of DNA. The DNA before a gene (in the 5' region) is considered "upstream", whereas the DNA after a gene (in the 3' region) is considered "downstream".

[^GC_pairing]: For this reason, the genetic composition of bacteria that live in hot springs is $80\%$ G-C.
[^enzyme]: An enzyme is a protein that speeds up a specific biochemical reaction.
[^polymerase]: DNA and RNA polymerases are particular proteins (or rather enzymes[^enzyme]) that synthesise nucleic acids.

## Transcription and RNA

```{figure} figures/transcription.svg
:name: fig:transcription
:align: center
:width: 600

A sketch showing three steps that lead to the synthesis of an mRNA: RNA polymerase (RNAp) attaches to the transcription start site (transcription initiation) and start building an RNA chain that is fully complementary to the DNA template strand (transcription elongation). When the gene has been fully transcribed, the RNAp and mRNA both detach and leave (transcription termination). Adapted from [here](https://commons.wikimedia.org/wiki/File:Simple_transcription_initiation1-fr.svg), [here](https://commons.wikimedia.org/wiki/File:Simple_transcription_elongation1-fr.svg), and [here](https://commons.wikimedia.org/wiki/File:Simple_transcription_termination1-fr.svg?uselang=fr).
```

When a gene product is to be expressed, the associated DNA is partially unwound to form a "bubble", where the two strands are separated. At this point transcription, which is the process by which RNA is produced using a DNA template, can commence: RNA polymerase is recruited to the transcription start site by regulatory protein complexes, and then reads the DNA from the 3' to 5' direction, placing down complementary bases to form a messenger RNA (mRNA) strand[^mRNA_direction]. RNA uses the same nucleotides as DNA, except Uracil is used instead of Thymine.

mRNA in eukaryotes experience post-translational modifications, or processes that edit the mRNA strand further, starting from the "pre-mRNA" molecule that is complementary to the template DNA strand of a gene. Most notably, a process called splicing removes introns, regions of the gene which don't code for protein, so that only the coding regions, the exons, remain. Different regions of the primary transcript may be spliced out to lead to different protein products (alternative splicing). In this way, an enormous number of different molecules may be generated based on different splicing permutations. In addition to splicing, both ends of the mRNA molecule are processed. The 5' end is capped with a modified guanine nucleotide, while at the 3' end, roughly 250 adenine residues are added to form a poly(A) tail. In most cases, the final molecule will be used by the cell machinery (and, in particular, by ribosomes) to synthesise proteins. As we will see in a moment, the base unit of a protein is an amino acid, and the information about what amino acid to add to the nascent protein is stored in the mRNA as a three-nucleotide sequence, also known as a *codon*.

RNA is not only the intermediary step to code a protein, but RNA molecules also have catalytic and regulatory functions. Though proteins are generally thought to carry out the main cellular functions, it has been shown that RNA molecules can have complex three-dimensional structures and perform diverse functions in the cell. In 1989 the Chemistry Nobel Prize was awarded to Thomas R. Cech and Sidney Altman for their "discovery of catalytic properties of RNA", which greatly contributed to the establishment of the concept of *ribozymes* (**ribo**nucleic acid en**zymes**), which are RNA molecules that have the ability to catalyze specific biochemical reactions.

The most common RNA types are:

1. mRNA (messenger RNA) contains the information to make a protein and is translated into protein sequence.
2. tRNA (transfer RNA) specifies codon-to-amino-acid translation. It contains a 3 base pair anti-codon complementary to a codon on the mRNA, and carries the amino acid corresponding to its anticodon attached to its 3' end.
3. rRNA (ribosomal RBA) forms the core of the ribosome, the organelle responsible for the translation of mRNA to protein.
4. snRNA (small nuclear RNA) is involved in splicing (removing introns from) pre-mRNA, as well as other functions.

The discovery of rybozymes contributed to the "RNA world" hypothesis, for which early life was based entirely on RNA. RNA served as both the information repository (like DNA today) and the functional workhorse (like protein today) in early organisms. Protein is thought to have arisen afterwards via ribosomes, and DNA is thought to have arisen last, via reverse transcription.

[^mRNA_direction]: Since the template strand is read 3' $\to$ 5', the nascent mRNA is built in the 5' $\to$ 3' direction. 

## Translation and the genetic code

Translation is the process by which the information encoded in mRNA is used to synthesize a protein. Unlike transcription, where the nucleotide sequence in DNA is transcribed into a corresponding nucleotide sequence in RNA, translation involves converting this nucleotide sequence into an amino acid sequence, forming the primary structure of the protein. Since there are 20 different amino acids and only 4 nucleotides, the mRNA is read in sets of three nucleotides, known as codons, each of which specifies a particular amino acid.

```{figure} figures/codon-flower.svg
:name: fig:codon-flower
:align: center

The meaning of each sequence made of three nucleotides (also known as a *codon* or *codon triplet*), as read by the ribosome. Each codon should be spelt starting from the inner circle and going outwards. The outer circle shows what corresponds to each codon.
```

[](#fig:codon-flower) shows the meaning of each codon. Since each nucleotide in a codon can have one out of 4 possible values (A, C, G or T), there exist $4^3 = 64$ possible codons. Out of these 64 combinations, three are recognised as "stop signals", while the rest encode amino acids, with the exception of the AUG codon which, in addition to represent methionine, is also typically used to initiate protein syntesis. Since there are 64 possible codon sequences, the code is degenerate, and most amino acids (excluding methionine and tryptophan) are encoded by more than a single codon. Most of the degeneracy occurs in the 3rd codon position.

The outer circle of [](#fig:codon-flower) contains the list of standard amino acids, together with their three- and one-letter abbreviations, which are shown in bold and after each AA name, respectively.

## Information flow: the central dogma of biology

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

## The experimental characterisation of the structure of biomolecules

I will refer to this part when talking about the secondary and tertiary structure of proteins and nucleic acids.

### Secondary structure

Several experimental methods can be used to understand the secondary structure of a particular biomolecule, such as proteins or nucleic acids. Here are some common techniques:

1. **Circular Dichroism (CD) Spectroscopy:**: CD spectroscopy is widely used to study the secondary structure of proteins. It measures the difference in the absorption of left-handed versus right-handed circularly polarized light, providing information about the overall content of $\alpha$-helices, $\beta$-sheets, nucleic acid helices, and random coils.

2. **Fourier Transform Infrared (FTIR) Spectroscopy:**: FTIR spectroscopy measures the absorption of infrared light by the biomolecule, which provides information about the types of chemical bonds and functional groups present. The amide I and II bands are particularly informative for analyzing protein secondary structures.
   
3. **Ultraviolet (UV) Absorption Spectroscopy:**: UV absorption spectroscopy can provide some information about the secondary structure, especially in nucleic acids. It is less informative for proteins compared to CD spectroscopy but can still be useful when combined with other techniques.

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

## Some useful numbers and quantities

| Quantity | Value |
| :--- | :---: |
| $k_B T$ at 298 K| $4.11 \times 10^{-21}$ J|
| " | $4.11 \, {\rm pn \cdot nm}$|
| RT at 298 K | $0.59 \, {\rm kcal / mol}$|
| " | $2.48 \, {\rm kJ / mol}$|
| Average nucleotide mass|$\approx 330 \to 300$ Da|
| Average amino acid mass|$\approx 110 \to 100$ Da|

:::{warning} Abuse of notation
When dealing with free energies and Boltzmann factors, it is practical to use the quantity $\beta \equiv 1 / k_B T$, which is expressed in units of measure of $J^{-1}$. However, in chemical physics and biophysics it is more common to use calories. As a result, I will frequently use $\beta$ to mean $1 / RT$.
:::

(sec:intro_comp_phys)=
# Introduction to computational physics

:::{warning} Using ChatGPT
You are welcome to use ChatGPT or similar tools during the course. However, harsh penalties will be handed out to students who use AI-assisted tools without actually understanding their answers. This applies to both code and text submitted for the homework assignments and for the final project.
:::

## Dynamic programming

:::{warning}
TODO
:::

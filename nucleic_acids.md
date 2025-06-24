---
title: Nucleic acids
exports:
   - format: pdf
---

```{tip}
The main references for this part are @lehninger2005lehninger, and @schlick2010molecular.
```

Here I take for granted that you are familiar with the general concepts introduced in the previous chapter regarding the interactions between atoms and molecules.

# The basic structure of DNA and RNA

Deoxyribonucleic acid (DNA) and ribonucleic acid (RNA) are the fundamental molecules that carry and execute the genetic instructions essential for all forms of life. Like proteins, they are macromolecules. However, the repeating monomer is not an amino acid, but a nucleotide, which is a molecule consisting of a a cyclic sugar (a [pentose](https://en.wikipedia.org/wiki/Pentose)), a phosphate group, and a nitrogenous base. A nucleotide without the phosphate group is called a nucleoside.

```{figure} figures/DNA_schematics.png
:name: fig:DNA_schematics
:align: center
:width: 500px

(a) The sugar and nitrogenous bases that make up DNA and RNA. Non-hydrogen atoms are labelled, and the broken lines indicate links that bind the compound to the other building blocks. Adapted from @schlick2010molecular. (b) A sketch showing two bonded DNA dinucleotides. The presence of a $-CH_3$ group connected to the pyrimidine ring highlighted in red makes the yellow-shaded base a thymine instead of a uracil. Adapted from [here](https://commons.wikimedia.org/wiki/File:0322_DNA_Nucleotides_Numbered.jpg).
```

In the classic DNA double helix (or *duplex*) described by [Watson and Crick](doi:10.1038/171737a0), a flexible ladder-like structure is formed by two chains (strands) that wrap around a virtual central axis. The two chains' backbones consist of alternating sugar and phosphate units, while the rungs that connect them consist of nitrogenous bases held together by hydrogen bonds and stabilised by stacking interactions. A schematic of a single nucleotide and of the Watson-Crick mechanism are shown in [](#fig:DNA_schematics).

:::{seealso} The discovery of the double helix

```{figure} figures/photo_51.jpg
:name: fig:photo_51
:align: center
:width: 300

The famous "Photo 51", taken 1952 by Raymond Gosling during his work with Rosalind Franklin. Credits to [Raymond Gosling/King's College London](http://www-project.slac.stanford.edu/wis/images/photo_51.jpg).
```

The discovery of the DNA double helix in 1953 marked a revolutionary advancement in molecular biology, primarily credited to James Watson and Francis Crick. Their success was built on the pivotal contributions of Rosalind Franklin and Maurice Wilkins, whose X-ray diffraction images of DNA provided critical evidence for the helical structure. In particular, these images, among which there was the famous "Photo 51", shown in [](#fig:photo_51), revealed the density and helical form of DNA. Wilkins shared these images, without her direct permission, with Watson and Crick who, thanks to their understanding of chemical bonding, enabled them to propose the accurate double helix model, fundamentally transforming our comprehension of genetic information storage and transmission.

Use [this page](https://www.roma1.infn.it/~dileorob/content/apps/photo51.php), prepared by prof. Roberto di Leonardo, to understand how the geometry of the double helix controls its X-ray diffraction pattern.

:::

Let's briefly discuss the main building blocks of a nucleotide.

## The pentose

A pentose is a simple sugar (*i.e.* a monosaccharide) containing five carbon atoms. In nucleic acids, pentose is present in its cyclic form, where the ring contains four carbons and one oxygen. In DNA, one of the hydroxyl ($-OH$) groups is replaced by a hydrogen atom. This is not a small change since, as we already know, a hydroxyl group can be involved in a hydrogen bond.

The three-dimensional conformation of the five-carbon sugar ring in nucleotides is non-planar, and the ring can adopt various forms, with the the C2'-endo and C3'-endo puckers being the most common forms for B-DNA and A-RNA (and A-DNA as well), which are helical structure that will be introduced in a moment. These puckering conformations, where the C2' or C3' atoms, respectively, are above the plane formed by the other four atoms, affect the overall structure and flexibility of DNA and RNA molecules and are one of the main source of conformational flexibility in nucleic acids (@schlick2010molecular).

## The phosphate group

A phosphate group is, by definition, the same for DNA and RNA and is composed by a phosphate atom attached to four oxygens. When in its anion form, a phosphate can be written as $[PO_4]^{-3}$: three of the oxygens are charged, and the fourth one is connected to the phosphate through a double bond. In nucleic acids, the phosphate group form a [phosphodiester bond](https://en.wikipedia.org/wiki/Phosphodiester_bond) with the 3' and 5' carbon atoms of the pentoses of consecutive nucleosides, joining them together to form the strand's backbone. I note here that you should not imagine the backbone as a rigid, fixed object, but as composed by segments around which the molecule can, to some extent, rotate. Indeed, there are six dihedral angles in the backbone (named, rather unoriginally, $\alpha, \beta, \gamma, \delta, \epsilon,$ and $\zeta$) that take different values depending on the local conformation (see @schlick2010molecular for more details).

Finally, as shown in [](#fig:DNA_schematics)(b), a phosphate that connects two nucleotides has a single residual negative charge. By contrast, phosphates at the beginning of a strand (also known as the 5' end) are doubly charged.

## The nitrogenous base

The nitrogenous bases endow the monomer with the exquisite selectivity that underlie the Watson and Crick double helix. Nitrogenous bases found in nucleic acids are 

* **Purines**, molecules with five- and a six-membered rings fused together: guanine (G) and adenine (A), present in both DNA and RNA.
* **Pyrimidines**, molecules with a single six-membered rings: cytosine (C), present in both DNA and RNA, thymine (T, DNA-only), uracil (U, RNA-only)[^RNA_only].

The chemical structure of the bases are shown in [](#fig:DNA_schematics)(b). Note that uracil has the same chemical structure as thymine, minus the $-CH_3$ group connected highlighted in red in the figure. Each base is connected to the C1' atom of the corresponding pentose through a nitrogen, which results in a *glycosyl* bond around which the base can rotate, contributing to the flexibility of the polymer. The dihedral angle associated to this rotation is called $\chi$.

In the classic pairing scheme unveiled by Watson and Crick, the geometry of the bases is such that G binds only to C, and A binds only to T (or U in RNA) by means of three and two hydrogen bonds, respectively. Additional non-canonical base pairings that are not only possible, but also biologically relevant exist. Here I will mention wobble base pairs, where G can bind with U[^wobble] with a thermodynamic stability that is similar to that of canonical base pairs, and Hoogsteen base pairs, where nucleotides can bind with an alternative geometry that can drive the formation of uncommon secondary structures such as [triple-stranded helices](https://en.wikipedia.org/wiki/Triple-stranded_DNA) or [G-quadruplexes](https://en.wikipedia.org/wiki/G-quadruplex).

[^wobble]: and also U, A and C can bind to the non-standard nucleotide [hypoxanthine](https://en.wikipedia.org/wiki/Hypoxanthine).
[^RNA_only]: Here I simplify, since chemical modifications are possible (*e.g.* methylation), and sometimes U and T can end up in DNA and RNA, respectively.

## Chain polarity

As mentioned earlier, the phosphodiester bond links the C5' atom of a nucleotide with the C3' atom of the one the follows it. As a result, DNA and RNA strands have a polarity, *i.e.* they are inherently directional. By convention, a strand starts at the $C5'-OH$ group (termed 5' end or terminal) and ends at the $C3'-OH$ group (the 3' end or terminal), and its sequence is also specified in this way. This is important since, as we will discuss in a moment, only anti-parallel strands (or anti-parallel sections of the same strand) can pair and form secondary structures.

## Hydrophobicity and hydrophilicity

The nitrogenous bases have both hydrophobic and hydrophilic regions: the aromatic ring structures have hydrophobic properties, but they also contain functional groups that can form hydrogen bonds, contributing to some degree of hydrophilicity. However, the partial hydrophobic character of the bases is more than counterbalanced by the backbone, which is highly hydrophilic due to the negatively charged phosphate groups and to the hydroxyl groups ($-OH$) of the sugar that can form hydrogen bonds with water. As a result, nucleotides readily dissolve in water.

The strong hydrophilic character of DNA and RNA strands makes their tertiary structure somewhat simpler compared to proteins, since once the [secondary structure elements](#sec:NA_secondary_structure) are formed, there is no strong interaction driving the formation of tightly packed super-secondary or tertiary structures[^NA_tertiary].

[^NA_tertiary]: Of course, this argument applies to solutions of nucleic acids only. In the cell, RNA, DNA, and proteins can and do interact together to form higher order structures.

(sec:hybridisation)=
## Hybridisation and denaturation

Two nucleotides coming together and binding to each other form a base pair (BP). Since hydrogen bonds are [highly directional](#sec:hydrogen-bonds), HB formation requires that the two nucleotides approach each other with the right mutual orientation. However, note that bases on their own (*i.e.* not part of a strand) can pair with any other base, as well as with themselves, often with more than one arrangement ([](doi:10.1016/S0079-6603(08)60565-6)). It is only in an extended paired region that the correct orientation is obtained only if the two nucleotides are compatible (*i.e.* if they can form a canonical or a non-canonical base pair, as mentioned above), and the strands run anti-parallel to each other. In this geometry, which is adopted by the most common DNA and RNA helical conformations,

* the bases on one strand can align perfectly with their complementary bases on the opposite strand, allowing stable hydrogen bonding;
* the sugar-phosphate backbones run in opposite directions, allowing the bases to stack neatly on top of each other, contributing to the overall stability of the double helix through van der Waals forces and hydrophobic interactions;
* negatively charged phosphate groups on the sugar-phosphate backbone are positioned in a way that minimizes repulsion.

Interestingly the spontaneous process through which single strands pair and form a duplex, which is called *hybridisation*, is now taken for granted, but was initially met with skepticism[^without_enzyme], since researchers at the time believed that duplex formation required enzymes to overcome electrostatic and entropic penalties ([](doi:10.1017/S0033583509004776)).

```{figure} figures/HB_stacking.png
:name: fig:HB_stacking
:align: center
:width: 600

A base pair step where hydrogen bonds and stacking interactions, the two main mechanisms responsible for hybridisation, are highlighted by arrows.
```

In addition to hydrogen bonding, duplexes are stabilised by interactions acting between the aromatic rings of consecutive nitrogenous bases (see [](#fig:HB_stacking)). These so-called base stacking interactions have a hydrophobic nature that is complemented by van der Waals and dipole-dipole termsthat contribute significantly to the stability and structural integrity of the DNA molecule. Rather counterintuitively, there is now ample evidence that base stacking is at least as important for duplex stability as base pairing, if not more (see *e.g.* [](doi:10.1093/nar/gkj454), [](doi:10.1016/j.plrev.2017.11.012), and [](doi:10.1038/s41565-023-01485-1)).

(box:stacking)=
:::{hint} Measuring the stacking strength

The classic study of [](doi:10.1016/j.jmb.2004.07.075) used a clever experimental setup based on the observation that a duplex where one of the composing strands is broken in two (*i.e.* it is *nicked*) has a mobility that is lower compared to an intact one. The reason is that a nicked strand can take two conformations: one straight, where all the stacking interactions that resembles that of an intact duplex are present, and one kinked, where all base pairs are formed but the stacking interaction at the nicked site is broken.

Assuming that

1. the stacking strength between two consecutive bases depends only weakly on whether the two bases are part of the same strand or not;
2. the nicked duplex is in the straight conformation with probability $p_s$ and it is kinked with probability $p_k = 1 - p_s$;
3. the mobility of the straight and kinked conformations are $\mu_s$ and $\mu_k$, respectively;

it is possible to write an expression for the mobility of a nicked duplex $\mu$ as a weighted average over the two conformations:

$$
\mu = p_s \mu_s + p_k \mu_k,
$$

which can be recast as

$$
\frac{p_k}{p_s} = \frac{\mu -\mu_s}{\mu_k - \mu}.
$$

Therefore, if $\mu_s$ and $\mu_k$ are known, a measure of $\mu$ can be linked to an estimate of the ratio between the equilibrium populations of the two states. However, in equilibrium the latter is connected to the free-energy difference between the two states, $\Delta G^{\rm ST}$ through the relation

$$
\frac{p_k}{p_s} = \exp\left( -\frac{\Delta G^{\rm ST}}{k_BT} \right).
$$

By changing the sequence of the nucleotides that flank the nicking site it is possible to obtain the values of $\Delta G^{\rm ST}$ for each dinucleotide step.

The same technique has been used to understand how the stacking parameters depend on temperature and ionic strength in [](doi:10.1093/nar/gkj454).
:::

The hybridisation of complementary strands is an enthalpy-driven process, and therefore it is promoted by a temperature decrease: at high temperature, where entropy dominates, the two strands are separated. As $T$ decreases, the energetic gain of forming BPs progressively leads to the stabilisation of the duplex. Strand hybridisation is a cooperative transition characterised by melting temperature and width that are controlled by the sequence, as well as by the buffer conditions (pH, ionic strength, *etc.*).

The transition is fully reversible: a duplex can be separated into its two composing strands (*i.e.* melted or denatured) by raising the temperature sufficiently.

[^without_enzyme]: "You mean [that a double helix can form] without an enzyme?" [](doi:10.1017/S0033583509004776)

(sec:NA_secondary_structure)=
# Secondary structure

Hydrogen bonding and base stacking drives are the main drivers for the formation of secondary structure in nucleic acids which, akin to the secondary structure of proteins, refers to the local conformation taken by polynucleotide chains. In the biological setting there is a marked difference between the possible secondary structures of DNA and RNA. Indeed, the secondary structure of biological DNA is rather boring, as it tends to form long and stable double-stranded helices storing the genetic code. By contrast, the secondary structure of biological RNA is more diverse, owing to the additional hydroxyl group present in the sugar, and to the fact that RNA is nearly always found as a single strand. Its complex secondary structure allows RNA to fulfill its many roles in processes such as catalysis, regulation, and protein synthesis.

:::{note}
The last two decades has seen the steady rise of DNA nanotechnology, which uses synthetic DNA to mimic, at least to some degree, the complexity of RNA secondary structure in order to exploit it for applications.
:::

Assuming that each base can either be unbound or involved in a single base pair[^no_triplex], the secondary structure of nucleic acids is defined by the pairing information of the nucleotides. By assigning a unique index to each nucleotide, this information can be represented as a list of index pairs, where each entry represents a single base pair. While, as I mentioned already, there are key differences between DNA and RNA, the surge of DNA nanotechnology has blurred the boundaries between the two, and many of the secondary structures that were unique to RNA can also be found in synthetic DNA systems. As a result, in this part I will refer to generic "nucleic acids", if not explicitly stated otherwise.

```{figure} figures/secondary_structure.png
:name: fig:secondary_structure
:align: center
:width: 500

The secondary structure of an RNA molecule. Dashed rectangles and labels have been added to the different elements. The double-stranded parts are coloured in green, while the other colours are used to highlight the other main secondary structures: 5' and 3' unpaired regions in orange, internal loops and bulges in yellow, internal loops of size one (also known as *mismatches*) in violet, hairpin loops in blue, multibranched loops in red. The graph representation has been generated with [forna](http://rna.tbi.univie.ac.at/forna/).
```

[](#fig:secondary_structure) shows the secondary structure of an RNA molecule drawn as a graph, which is a common representation. The backbone runs from the 5' to the 3' end and it is drawn with grey sticks representing the graph's edges, while the single nucleotides are the vertices, drawn as circles. Base pairs are depicted with red connections. The RNA molecule has been designed to feature the main secondary structure motifs, which are highlighted by coloured dashed rectangles. These are:

* **Double-stranded part**: Sections made of consecutive base pairs.
* **Hairpin** (also known as a *stem-loop* structure): A single strand that folds back on itself to form a stem-loop structure. The stem is a double-stranded region where bases pair with complementary bases, and the loop is a single-stranded region at the tip.
* **Bulge**: Occurs when one or more unpaired nucleotides bulge out from one side of a double-stranded stem.
* **Internal loop**: Formed when there are unpaired nucleotides on both sides of the double-stranded stem, creating a loop in the middle of the stem. An internal loop made of one nucleotide on each strand is often called a **mismatch**.
* **Multibranched loop**: A loop where three or more double-stranded stems converge.
* **Pseudoknot**: Formed when bases in a loop pair with complementary bases that are separated by at least one double-stranded section. This results in a crossover of strands, forming a knot-like structure. A common way a pseudoknot forms is when a single-stranded loop of a hairpin pairs with a complementary sequence outside the hairpin.

:::{hint} The structure of a biological RNA molecule

Transfer RNA (tRNA) is a type of short (76 to 90 nucleotides) RNA molecule that plays a crucial role in the process of translating genetic information from messenger RNA (mRNA) into proteins. tRNA molecules function as adaptors that match specific amino acids to corresponding codons in the mRNA sequence during protein synthesis.

In particular, the primary function of tRNA is to translate the genetic code from mRNA into a sequence of amino acids, ultimately forming a protein. Each tRNA molecule is "loaded" with its corresponding amino acid by an enzyme called aminoacyl-tRNA synthetase[^aaRS], forming an aminoacyl-tRNA complex. The anticodon of the tRNA pairs with the complementary codon on the mRNA strand in the ribosome, whose machinery catalyses the formation of the peptide bond through which the amino acid carried by the tRNA is added to the growing polypeptide chain.

```{figure} figures/tRNA.png
:name: fig:tRNA
:align: center
:width: 600

The (a) secondary and (b) three-dimensional structure of a tRNA ([yeast phenylalanine tRNA, 1EHZ](https://www.rcsb.org/structure/1EHZ)). The specific secondary motifs are coloured differently in (b), and highlighted by coloured dashed rectangles in (a). The nucleotides that bear chemical modifications are shown in blue[^modified_nucleotides].
```

[](#fig:tRNA) presents the secondary and three-dimensional structure of a tRNA molecule, highlighting its [cloverleaf](https://en.wikipedia.org/wiki/Four-leaf_clover) structure. It is worth mentioning that, as shown in the figure, many of the nucleotides of a tRNA molecule are chemically modified in order to tune their interaction with the ribosome or with the mRNA codons. The four main parts of a tRNA molecule are:

1. **The Acceptor Stem:** This is the site where a specific amino acid is attached. The 3' end of the tRNA has a conserved sequence (CCA) where the amino acid binds.
2. **The Anticodon loop:** This contains a sequence of three nucleotides (the anticodon) that is complementary to a specific mRNA codon, ensuring the correct amino acid is added to the growing polypeptide chain.
3. **The D loop:** Named for the presence of [dihydrouridine](https://en.wikipedia.org/wiki/Dihydrouridine), this arm is involved in tRNA folding.
4. **The T$\psi$C loop:** Named for the conserved sequence of thymine, [pseudouridine](https://en.wikipedia.org/wiki/Pseudouridine), and cytosine, this arm is important for tRNA recognition by the ribosome.

[^aaRS]: There is one aminoacyl-tRNA synthetase for each amino acid of the genetic code: in humans, there are 20.
[^modified_nucleotides]: **m$_2^G$**: 2-methyl-guanosine, **D**: 5,6-Dihydrouridine, **m$_2^2$G**: N2-dimethylguanosine, **C$_m$**: O2'-methyl-cytdine, **G$_m$**: O2'-methyl-guanosine, **T**: 5-Methyluridine (Ribothymidine), **Y**: wybutosine (Y-base), $\psi$: pseudouridine, **m$^5_C$**: 5-methyl-cytidine, **m$^7_G$**: 7-methyl-guanosine, **m$^1_A$**: 1-methyl-adenosine.
:::

(sec:dot-paren)=
## Dot-paren notation

Dot-paren notation provides a clear and concise way to represent the secondary structure of nucleic acids, showing which nucleotides are paired or unpaired, regardless of their identity. It is a useful tool for storing information about complex secondary structures, and for exchanging this information to and from computational tools.

:::{warning}
Dot-paren notation is primarily designed for representing the secondary structure of single-stranded nucleic acids, where it is very effective in illustrating the pattern of base pairing. However, for multi-strand systems, the dot-paren notation can become limited and cumbersome.
:::

The dot-paren notation uses dots and parentheses to indicate unpaired and paired nucleotides, respectively.

1. **Unpaired Nucleotides:**
   - Dots (.) represent nucleotides that are not involved in base pairing.
   - Example: `....` represents four unpaired nucleotides.
2. **Paired Nucleotides:**
   - Parentheses indicate base pairs, with matching opening and closing parentheses representing the paired nucleotides.
   - An opening parenthesis `(` signifies the 5' end of a base pair, while a closing parenthesis `)` signifies the 3' end.
   - A well-formed dot-paren representation should have the same number of opening and closing parentheses.
3. **Hairpins:**
   - Hairpin loops are represented with a series of dots inside parentheses.
   - Example: `(((....)))` indicates a stem of three base pairs with a loop of four unpaired nucleotides.
4. **Bulges and Internal Loops:**
   - Bulges are unpaired nucleotides on one side of a stem.
   - Internal loops have unpaired nucleotides on both sides.
   - Example: `((..((...))..))` represents a structure where there is a small internal loop of size 2 within a larger stem, while `((.(...)))` represents a stem-loop with a single bulge.
5. **Multibranch Loops and Pseudoknots:**
   - Multibranch loops involve multiple stems and are, in general, complex to represent.
   - Pseudoknots, which involve base pairs crossing one another, are not well-represented by simple dot-paren notation but can be indicated with additional notation if necessary (see the complex example in the [box below](#box:dot_paren_complex)).

:::{hint} Two simple examples
Consider the following RNA sequence, which forms a hairpin with stem and loop of size 4:

```
5'- GGCAUCGUUGCC -3'
```

The dot-paren notation for this structure is:
```
((((....))))
```

If a nucleotide is added in third position the sequence will also feature a bulge. An example of such a sequence and its secondary structure expressed with dot-paren notation would be

```
5'- GGCCAUCGUUGCC -3'
    ((.((....))))
```
:::

(box:dot_paren_complex)=
:::{hint} A complex example: the input used to generate [](#fig:secondary_structure)

[](#fig:secondary_structure) has been generated using the following prompt, where the first line is the sequence and the second line is the dot-paren representation of the associated secondary structure:

```
CGCUUCAUAAUGCACAUCCUCAAGCUGAUAGUGUGCUUGGGAAUGUCUGCACCAAGAGCCUUAAACUCUUGUGAUUAUGAAGUG
...((((((.(((...((((.((((.....[[[.)))))))).]]]..))).((((((.......))))))....))))))...
```

Note the square brackets: this is a (non-standard but common) way of representing pseudoknots.
:::

There are other formats to store nucleic acid secondary structures that overcome some of the problems with the dot-paren notation (pseudoknot representation, lack of multi-strand support, poor readability, *etc.*), but none has become a standard yet. You can look at some examples [here](https://viennarna.readthedocs.io/en/latest/io/rna_structures.html) or [here](https://it.mathworks.com/help/bioinfo/ug/predicting-and-visualizing-the-secondary-structure-of-rna-sequences.html).

[^no_triplex]: This is already an approximation, since some non-canonical (*e.g.* Hoogsteen) base pairings can connect a nucleotide to another which is already involved in a base pair.

(sec:NN_models)=
## Nearest-neighbour models

Nearest-neighbour (NN) models for nucleic acids are computational frameworks used to predict the thermodynamic stability and secondary structure of DNA and RNA molecules. These models operate on the principle that the stability of a base pair is influenced primarily by its immediate neighbors rather than by more distant sequences. As we [already know](#sec:hybridisation), base stacking and hydrogen bonding are the main drivers for the formation of helical structures. Therefore, assuming that the stability of nucleic acid structures is determined primarily by the local sequence context (the identity and orientation of adjacent base pairs), the complex interactions within nucleic acids can be simplified by considering only the contributions of dinucleotide pairs, hence the "nearest-neighbour" name.

The most used NN models are:

* The SantaLucia model for DNA ([](doi:10.1146/annurev.biophys.32.110601.141800))
* The Turner model for RNA ([](doi:10.1093/nar/gkl472))

In addition, [here](https://rna.urmc.rochester.edu/NNDB/) is a useful website where several NN models are briefly described, and their parameters can be downloaded.

In general, a NN model is defined by a list of contributions that make it possible to assign a free-energy cost to the secondary structure of a specific strand (or system of strands) in an additive way: the total free-energy cost of the structure is given by a sum of terms that refer to the sequence and type of each local secondary structure. The specific values that enter into the calculations are being constantly improved upon by means of careful experiments on many different sequences (similar in spirit to  those performed to obtain the sequence-dependent stacking strength, see [the box above](#box:stacking)). A recent example is [](doi:10.1093/nar/gkac261)). By contrast, the functional forms and the nature of the different free-energy terms are rather stable and did not change much in the last 20+ years.

The contributions are given in terms of $\Delta H^\circ$ and $\Delta S^\circ$ or $\Delta G^\circ$ and $\Delta H^\circ$, which are linked by the relation

$$
\Delta G^\circ = \Delta H^\circ - T \Delta S^\circ.
$$

where the $^\circ$ superscript signals that these values refer to the free-energy differences estimated at the "standard" strand concentration of 1 molar (*i.e.* one mole per liter), $C^\circ$. If the strand concentration $C$ is different from $C^\circ$[^standard_conc], the final free-energy difference contains the additional entropic term

$$
\Delta S_C = - R \log C,
$$

where $R \approx 2$ cal / mol K is the gas constant.

We now analyse the main free-energy contributions to the formation of secondary structures used in NN models.

[^standard_conc]: In DNA nanotechnology the strand concentration is often of the order of $10^{-9}$ M.

### Nearest-Neighbor Interactions

Double-strand formation is driven by the combined effects of base stacking and hydrogen bonding between adjacent base pairs. In NN models, these are accounted for by summing up the enthalpy and entropy contributions of each dinucleotide base step. Therefore, for a sequence of length $N$ there are $N - 1$ terms. As an example, consider the two fully-complementary strands

```
5'- TACCTG -3'
3'- ATGGAC -5'
```

In order to estimate the free-energy, we split the sequence into dinucleotide steps:

```
5'- TA AC CC CT TG -3'
3'- AT TG GG GA AC -5'
```

so that the total contributions due to this term are given by

$$
\begin{aligned}
\Delta H & = \Delta H_{\rm TA/AT} + \Delta H_{\rm AC/TG} + \Delta H_{\rm CC/GG} + \Delta H_{\rm CT/GA} + \Delta H_{\rm TG/AC}\\
\Delta S & = \Delta S_{\rm TA/AT} + \Delta S_{\rm AC/TG} + \Delta S_{\rm CC/GG} + \Delta S_{\rm CT/GA} + \Delta S_{\rm TG/AC}
\end{aligned}
$$

:::{warning}
Remember that DNA and RNA strands have a **polarity**: the contributions due to CT/GA and TC/AC base steps differ! Always arrange the strands so that the top strand is listed in the 5' $\to$ 3' direction and the bottom strand in the 3' $\to$ 5' direction before splitting the sequence.
:::

### Terminal penalty

Most NN models requires adding a correction term to the free energy that depends on the base pair at each end of the helix. In most cases the penalty is present only if the terminal base pair is AT (or TA) and it tends to be rather small. Therefore, it is important only for short sequences.

### Initiation and Symmetry Terms

The formation of a double helix requires an initial free-energy input to overcome entropy and start the pairing process. The term has to be added once per contiguous double-stranded region or helix.

In the presence of self-complementary duplexes (*i.e.* if the sequence of the two strands composing a duplex is palindromic), there is an additional purely entropic term that takes into account the fact that a strand can pair with a copy of itself. Here is an example:

```
5'- AGCGCT -3'
3'- TCGCGA -5'
```

### Loop Free Energies

Forming hairpins, bulges, external, internal or multibranched loops costs both entropy (since we constrain the strand backbone) and enthalpy (since the backbone has to be bent or twisted to some extent). These free-energy penalties are, in general, sequence-dependent, but the number of possible combinations grows exponentially with the loop size. Therefore, the free-energy cost of these motifs is usually approximated to depend only on the loop size (with some exceptions such as hairpins with loops of length three and four, see *e.g.* [](doi:10.1146/annurev.biophys.32.110601.141800)).

(sec:mismatches)=
### Penalties for Mismatches and Terminal Mismatches

Non-complementary paired bases have a destabilising effect on secondary structure. The specific energetic penalty that applies depend on the sequence of the mismatch, but also whether it occurs at the ends of helices (terminal mismatch), or within the helix (internal mismatches, sometimes referred to as "1x1 internal loops". Both reduce the overall stability of the structure, but to different extents. The following example contains both types of mismatches (look at the first and fifth base pairs):

```
5'- CTACACTG -3'
3'- TATGCGAC -5'
```

### Dangling Ends

Dangling ends refer to unpaired nucleotides at the 5' or 3' ends of a helix that can stabilise specific secondary structures such as multibranched and exterior loops through additional stacking interactions. For instance, in this example the top strand has both a 5' and a 3' dangling end:

```
5'- ATACCTGC -3'
3'-  ATGGAC  -5'
```

These terms tend to be sequence-dependent: in the example above, the contribution of the 5' dangling end would be different if the sequence was TT/A or AC/A instead of AT/A. Note that a helix end extended on both strands has a [terminal mismatch](#sec:mismatches) rather than two dangling ends.

### Coaxial Stacking Parameters

Coaxial stacking refers to the stacking interactions between adjacent helices in branched or complex secondary structures, such as those found in multibranched loops or multi-strand systems. These interactions can significantly stabilize the overall structure by allowing helices to stack on top of each other in a manner similar to base stacking within a single helix.

In most NN models, two types of coaxial stacking are handled: when two helices are directly adjacent and no intervening unpaired nucleotides are present ("flush coaxial stacking"), or when a single mistmatch occurs between the stacked helices ("mismatch-mediated coaxial stacking"). The following examples (taken from [here](https://rna.urmc.rochester.edu/NNDB/rna_2004/rna_2004_coaxial_stacking.html)) show the two cases, where the backbone is explicitly drawn using dashes:

:::{card} Flush coaxial stacking
```
5'- C-A-G-A -3'
3'- G-U C-U -5'
      | |
      G A
      | |
      5'3'
```
:::

:::{card} Mismatch-mediated coaxial stacking
```
5'- C-A-G-A -3'
3'- G-U A-U -5'
      | |
      G A
      | |
      5'3'
```
:::

### Salt-Dependent Terms

The stability of nucleic acid structures is influenced by the ionic environment, since cations like Na$^+$ and Mg$^{2+}$ shield the negative charges on the phosphate backbone and reduce electrostatic repulsion between strands. Note that, as far as I know, only the SantaLucia NN model for DNA, presented in [](doi:10.1146/annurev.biophys.32.110601.141800), takes into account this contribution through the following entropic term:

$$
\Delta S_{\rm salt} = 0.368 \frac{N_p}{2} \log C_S,
$$

where $N_p$ is the number of phosphates in the duplex, so that $N_p / 2 = N$ is, under usual conditions, the duplex length, $C_S$ is the molar concentration of monovalent cations[^magnesium], and the resulting contribution is in units of cal / mol K.

[^magnesium]: Magnesium and other multivalent cations are not supported.

(sec:two_state)=
## The two-state model

### Duplex formation

One of the most straightforward applications of any NN model is to model the thermodynamics of duplex formation in a system composed of just two (perfectly or partially) complementary strands, A and B. If the sequences are such that the possibility of stable intermediates can be neglected (that is, if the strands spend most of the time either free in solution or bound to each other), hybridisation can be described as a two-state process. The latter, in turn, is formally equivalent to the chemical equilibrium between two reactants and a product, *viz.*

$$
A + B \rightleftharpoons AB.
$$

The two-state model works particularly well for short strands (*i.e.* oligonucleotides) since, if they do not contain repeating patterns or similar "pathological" sequences, are more unlikely to exhibit metastable intermediates, *i.e.*, states with secondary structures whose stability can almost match that of the product.

Under the two-state assumption, equilibrium is described by the law of mass action equation

$$
K_{AB} = \frac{C_A C_B}{C_{AB}},
$$

where $C_X$ is the (molar) concentration of $X$ and $K_{AB}$ is the dissociation constant associated to the reaction, which is in turn connected to the free-energy difference between the two states, $\Delta G_{AB}^\circ$, through

$$
K_{AB} = C^\circ \exp\left( \beta \Delta G_{AB}^\circ \right),
$$

hence

$$
\frac{C_A C_B}{C_{AB}} = C^\circ \exp\left( \beta \Delta G_{AB}^\circ \right).
$$ (eq:equilibrium)

:::{hint} The dissociation constant
The dissociation constant, often denoted as $K_d$​, is a measure of the affinity between two molecules in a binding interaction, such as between between DNA strands in a duplex or a protein and its ligand. It is expressed in units of concentration (typically molarity, M), representing the concentration at which half of the strands are bound in duplexes, or half of the binding sites are occupied by the ligand. A lower $K_d$​ value indicates higher affinity, meaning the molecules are more likely to remain bound together at lower concentrations. Conversely, a higher $K_d$ suggests weaker binding, indicating that higher concentrations are needed for significant binding to occur.
:::

We define the concentrations of the two isolated strands (*i.e.* when $C_{AB} = 0$) as $C_{A,0}$ and $C_{B,0}$, so that the total strand concentration is $C_0 \equiv C_{A,0} + C_{B,0}$. Without loss of generality we assume that $C_{A,0} \leq C_{B,0}$.

Every time a duplex forms, the number of A and B strands decreases by one each, so that the total strand concentration can be written as $C_0 = C_A + C_B + 2 C_{AB}$. We first compute the *melting temperature* $T_m$, which is the temperature at which half of the duplexes are formed. Under this condition, $C_{AB} = C_A$, since the number of duplexes that can form is controlled by the number of strands in the minority species, so that Eq. [](#eq:equilibrium) becomes

$$
\frac{C_B}{C^\circ} = \exp\left( \frac{\Delta G_{AB}^\circ}{R T_m} \right),
$$

which, recalling that $\Delta G_{AB}^\circ = \Delta H_{AB}^\circ - T \Delta S_{AB}^\circ$, can be used to obtain the following expression for $T_m$:

$$
T_m = \frac{\Delta H_{AB}^\circ}{\Delta S_{AB}^\circ + R \log\left( \frac{C_B}{C^\circ} \right)}.
$$

However, since $C_B = C_0 - C_A - 2 C_{AB} = C_0 - 3 C_A = C_{A,0} + C_{B,0} - 3 C_A$, and, by definition of melting temperature, $C_A = C_{A,0} / 2$, so that $C_B = C_{B,0} - C_{A,0} / 2$, we find

$$
T_m = \frac{\Delta H_{AB}^\circ}{\Delta S_{AB}^\circ + R \log\left( \frac{2C_{B,0} - C_{A,0}}{2 C^\circ} \right)}.
$$ (eq:T_m_general)

In the common case of equimolarity, *i.e.* when $C_{B,0} = C_{A,0} = C_0 / 2$, Eq. [](#eq:T_m_general) simplifies to

$$
T_m = \frac{\Delta H_{AB}^\circ}{\Delta S_{AB}^\circ + R \log\left( \frac{C_0}{4 C^\circ} \right)}.
$$ (eq:T_m_equi)

:::{warning} The $C^\circ$ factor
In most of the cases, the $C^\circ$ factor in Eq. [](#eq:T_m_equi) is omitted, since $C^\circ = 1$ M. However, as a physicist, you should always distrust equations where the arguments of mathematical functions such as $\log$, $\exp$, $\cos$, *etc.* have physical dimensions. As a rule, you can disregard dimension issues only if you understand where they come from.
:::

The denominator of Eq. [](#eq:T_m_equi) can be rewritten as $\Delta S_{AB}^\circ + R \log\left( \frac{C_0}{4 C^\circ} \right) = \Delta S_{AB}^\circ + R \log\left( \frac{C_{0,A}}{2 C^\circ} \right) = \Delta S_{AB} - R \log 2$, where

$$
\Delta S_{AB} \equiv \Delta S_{AB}^\circ + R \log \left( \frac{C_{0,A}}{C^\circ} \right)
$$ (eq:renormalised_entropy)

is a renormalised entropy difference that takes into account the concentration at which the reaction takes place. Using Eq. [](#eq:renormalised_entropy) makes it possible to directly derive the free-energy difference $\Delta G_{AB} = \Delta H_{AB}^\circ - T \Delta S_{AB}$ between the $A + B$ and $AB$ states at any strand concentration. In the general case $C_{A,0} \leq C_{B,0}$, the additional entropic factor is $R\log \left( \frac{2C_{B,0} - C_{A,0}}{C^\circ} \right)$.

```{figure} figures/SantaLucia.png
:name: fig:SantaLucia
:align: center
:width: 600px

Experimental versus predicted melting temperatures of DNA oligonucleotides. (a) Data for 264 duplexes of length 4 to 16 bp in solution with 1 M of NaCl. (b) Data for 81 duplexes of length 6 to 24 bp in solution with variable NaCl concentration ranging from 0.01 to 0.5 M. Adapted from [](doi:10.1146/annurev.biophys.32.110601.141800).
```

[](#fig:SantaLucia) shows a comparison between the experimental and theoretical melting temperatures of hundreds of DNA oligonucleotides predicted with the [SantaLucia model](doi:10.1146/annurev.biophys.32.110601.141800). The average absolute deviation is smaller than $2.3$ K.

(sec:hairpin_formation)=
### Hairpin formation

Hairpin formation can also be modelled as a two-state process, where the equilibrium is between the open (random coil) and closed (stem-loop or hairpin) states:

$$
{\rm Coil} \rightleftharpoons {\rm Hairpin}.
$$

In this case the overall strand concentration does not play any role[^hairpin_conc], and the dissociation equilibrium is given by

$$
\frac{C_c}{C_h} = e^{\beta \Delta G_{ch}^\circ},
$$

where $C_c$ and $C_h$ are the concentrations of strands in the coil and hairpin conformations, respectively, and $\Delta G_{h}^\circ$ is the free-energy difference between the two states. The condition for the melting temperature is $C_c = C_h$, which yields

$$
T_m = \frac{\Delta H^\circ}{\Delta S^\circ}.
$$

[^hairpin_conc]: If we can neglect the interaction between different strands, *i.e.* if the overall concentration is low enough that we can assume ideal gas behaviour.

### Melting curves

We now consider a generic nucleic acid system where one or more strands can pair and/or fold into a product P. We define the melting curve as the yield of P, in terms of concentration or fraction of formed product, as a function of temperature or another thermodynamic parameter that is changed experimentally, such as pH or salt concentration.

I will now show how the same two-state formalism introduced ealier can be used to predict the melting curve of a system. For the sake of simplicity I will use the $A + B \rightleftharpoons AB$ system with $C_{A,0} = C_{B,0}$. Defining $X$ as the probability that a strand is *not* part of a duplex, the equilibrium concentrations can be written as $C_A = C_B = X C_{A,0}$ and $C_{AB} = (1 - X) C_{A,0}$. Substituting these relations in Eq. [](#eq:equilibrium) and simplifying common factors we find

$$
\frac{X^2}{1 - X} = \frac{C^\circ}{C_{A,0}} e^{\beta \Delta G_{AB}^\circ} = e^{\beta \Delta G_{AB}},
$$

where we have used the renormalised entropy to obtain the latter relation. Resolving for $X$ and taking the positive root we find

$$
X = \frac{-e^{\beta \Delta G_{AB}} + \sqrt{e^{2\beta \Delta G_{AB}} + 4e^{\beta \Delta G_{AB}}}}{2} = \frac{-1 + \sqrt{1 + 4 e^{-\beta \Delta G_{AB}}}}{2 e^{-\beta \Delta G_{AB}}}.
$$

From $X$ the duplex yield can be computed as $p_b = 1 - X$ (which is the probability that a strand is in a duplex), and the duplex concentration as $C_{AB} = p_b C_{A,0}$.

```{figure} figures/DNA_melting_comparison.png
:name: fig:DNA_melting_comparison
:align: center
:width: 250px


The yield of a 10-bp duplex as predicted by SantaLucia (red line), and two coarse-grained simulation models, oxDNA and oxDNA2. Adapted from [](doi:10.1063/1.4921957).
```

[](#fig:DNA_melting_comparison) shows the yield of a oligonucleotide of 10 bp as predicted by the SantaLucia NN model and by simulations of two coarse-grained models.

# Tertiary structure

The tertiary structure of nucleic acids is, in general, much simpler than that of proteins. This is due to the more limited variety of building blocks involved (4 *vs* 20), and to the charged (and, in general, hydrophilic) nature of the DNA and RNA backbones, which tend to destabilise the type of "super-secondary structures" that are so common in proteins. Moreover, the lack of tightly-packed structures blurs the difference between secondary and tertiary structures, since most of the times, especially in simple systems, the tertiary structure is straightforwardly implied by the secondary structure.

In general, the most dominant tertiary structure of RNA and DNA is the double helix, and is the main one we will consider going forward. Many different helical structures have been discovered in biological contexts or synthesised artificially under specific conditions. However, here I will present only the three main ("canonical") ones. Note that they I will provide only an average characterisation for each helix, as it is known that their properties can depend (even strongly in some cases) on the local sequence. You can use [this webserver]() (which is presented in [](doi:10.1016/j.jmb.2023.167978)) to visualise the 3D structures of a given sequence, as estimated with a coarse-grained, realistic model.

(sec:canonical_helices)=
## Canonical helices

When two complementary strands pair together, they wrap around each other and form a periodic helical structure. The type of the resulting helix depends, in general, on the strand type (DNA or RNA), external conditions (ionic strength, pH, water concentration, *etc.*), and sequence.

### B-DNA

```{figure} figures/bdna.png
:name: fig:bdna
:align: center
:width: 500px

A B-DNA double helix composed of 12 base pairs as seen from the side, (a) and (b), and from the top (c). In (a) the backbone is represented as a ribbon, the sugar as a pentagon and the bases are outlined as pentagons and hexagons, while (b) and (c) shows the atoms as van der Waals spheres. The structure has been generated with the [3DNA 2.0 webserver](http://web.x3dna.org).
```

The single most important conformation of DNA is B-DNA, which is the famous double helix whose structure was described for the first time by [Watson and Crick](doi:10.1038/171737a0). B-DNA is a right-handed helical structure that consists of about 10.5 base pairs per turn, with a diameter of approximately 2 nanometers. The structure features major and minor grooves, which are essential for protein-DNA interactions; the major groove is wide and deep, while the minor groove is narrow and shallow. B-DNA is most stable under physiological conditions, making it the prevalent form in biological systems. As described in the box above, the famous experimental photo that was instrumental for Watson and Crick's discovery was taken by a student of [Rosalind Franklin](https://en.wikipedia.org/wiki/Rosalind_Franklin), who gave its name not only to B-DNA, but also to A-DNA, which is another canonical form of DNA.

[](#fig:bdna) shows a perfect B-DNA helix with different representations and from two different points of view. In panel (a) the representation makes it easy to see that

1. The strands run in an anti-parallel fashion: look at the orientation of facing pentagons, or at the way bases are connected to the sugars.
2. Base pairs have essentially zero inclination with respect to the helical axis.

### A-DNA and A-RNA

```{figure} figures/adna.png
:name: fig:adna
:align: center
:width: 500px

An A-DNA double helix composed of 12 base pairs represented as in [](#fig:bdna). The structure has been generated with the [3DNA 2.0 webserver](http://web.x3dna.org).
```

DNA turns into its A form when dyhdrated, which can happen as a result of human intervention (*e.g.* when crystalling samples to study with X-rays), but also under natural conditions, when water is scarce or specific proteins bind to the DNA, in place of some of the water molecules[^A-DNA]. A-DNA is a right-handed helix that is more compact than B-DNA: it has about 11 base pairs per turn and a diameter of approximately 2.3 nanometers. A-DNA features deeper major grooves and shallower minor grooves compared to B-DNA. The helical structure is more tightly wound, with the base pairs tilted by abouy $20^\circ$ relative to the helix axis, and the distance between base pairs is shorter. Moreover, the sugars are in the C3'-endo conformation rather than in the C2'-endo conformation like in B-DNA. A-DNA is less prevalent in biological systems but can be found in certain DNA-RNA hybrids. Additionally, double-stranded RNA molecules or regions typically adopt the A-form helical structure.

[](#fig:adna) shows a perfect A-DNA helix made of the same number of base pairs as the one shown in [](#fig:bdna). Compared to B-DNA, it is evident that A-DNA is shorter and thicker, and the base pairs are inclined by $\approx 20^\circ$ with respect to the helical axis.

[^A-DNA]: In both cases the A form [is thought to protect the genetic material from extreme conditions](https://en.wikipedia.org/wiki/A-DNA#Biological_function).

### Z-DNA

```{figure} figures/zdna.png
:name: fig:zdna
:align: center
:width: 600px

A Z-DNA double helix composed of 12 GC base steps (*i.e.* 24 base pairs) represented as in [](#fig:bdna). Note the zig-zag pattern drawn by the backbone, which gives Z-DNA its name. The structure has been generated with the [3DNA 2.0 webserver](http://web.x3dna.org).
```

Another important DNA conformation is Z-DNA, which is a less common and structurally distinct form of DNA characterized by a left-handed double helical structure, in contrast to the right-handed helices of A-DNA and B-DNA. It has a zigzag backbone, which gives Z-DNA its name, and consists of about 12 base pairs per turn with a diameter of approximately 1.8 nanometers. The major and minor grooves are less distinct in Z-DNA, and the structure is more elongated and thinner. Z-DNA requires a specific sequence of alternating purines and pyrimidines (*e.g.* repeating GC steps) and can form under high salt conditions or negative supercoiling, *i.e.* twisting of the DNA in the left-handed direction. Though less prevalent than B-DNA, Z-DNA is thought to have important biological functions in gene regulation, genetic recombination, and in the context of certain protein-DNA interactions.

[](#fig:zdna) shows a perfect Z-DNA helix made of 24 base pairs (12 repeating GC steps) which shows the main features of Z-DNA: the left handedness of the helix, the zig-zag pattern that gives this conformation its name, and the small BP inclination.

A comparison between the average properties of the three main helical conformations is shown in [](#tbl:helices).

:::{table} The main (average) properties of the three canonical forms of DNA and RNA.
:label: tbl:helices
:align: center

| Form | Handedness | BPs per turn | Diameter ($\angstrom$) | Pitch ($\angstrom$) | BP inclination | Sugar pucker |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| B-DNA | Right | 10 | 20 | 34 | 0$^\circ$ | C2'-endo |
| A-DNA/A-RNA | Right | 11 | 26 | 28 | 20$^\circ$ | C3'-endo |
| Z-DNA | Left | 12 | 18 | 45 | -7$^\circ$ | alternating[^Z-DNA_pucker] |
:::

[^Z-DNA_pucker]: The alternating pyrimidines and purines take the C2'-endo and C3'-endo conformations, respectively.

## Other tertiary motifs

Other important tertiary motifs are those that can connect more than two strands together, or two sections of the same strand that are far apart from each other. For instances, in triplexes (*i.e.*  triple-stranded DNA or RNA) a third strand can bind in the major groove of a duplex through Hoogsteen base pairings, or in RNA in the minor groove by leveraging the presence of the additional hydroxyl group in the sugar.

Hoogsteen base pairings can also lead to the formation of quadruplexes, which can occur in a variety of patterns, including intramolecular (within a single strand), intermolecular (between different strands), and hybrid types. A strand (or multiple strands) with a high number of consecutive guanine bases folds back on itself (or aligns with other strands) to bring the G bases into proximity. Four guanine bases form a planar structure known as a G-tetrad through Hoogsteen hydrogen bonding, and multiple G-tetrads stack on top of each other, stabilized by the $\pi-\pi$interactions between the aromatic rings of the guanine bases.

Finally, as mentioned earlier, pseudoknots and coaxial stacking interactions can also seen as mechanisms the can lead to the formation of particular tertiary structures, since they can stabilise multi-strand (or multi-loop) structures.

---
title: Nucleic acids
---

```{tip}
The main references for this part are @lehninger2005lehninger, and @schlick2010molecular.
```

# The basic structure of DNA and RNA

Deoxyribonucleic acid (DNA) and ribonucleic acid (RNA) are the fundamental molecules that carry and execute the genetic instructions essential for all forms of life. Like proteins, they are [macromolecules](./proteins#macromolecules). However, the repeating monomer is not an amino acid, but a nucleotide, which is a molecule consisting of a a cyclic sugar (a [pentose](https://en.wikipedia.org/wiki/Pentose)), a phosphate group, and a nitrogenous base. A nucleotide without the phosphate group is called a nucleoside.

```{figure} figures/DNA_schematics.png
:name: fig:DNA_schematics
:align: center
:width: 500px

(a) The different parts that make up a nucleotide, which is an RNA or DNA nucleotide depending on whether the bottom-right group of the sugar is a hydroxyl or a hydrogen, respectively. Here the nitrogenous base is a guanine, and the picture also shows the labeling of the sugar and base atoms. (b) A sketch showing two bonded DNA dinucleotides. The presence of a $-CH_3$ group connected to the pyrimidine ring highlighted in red makes the yellow-shaded base a thymine instead of a uracil. Adapted from [here](https://commons.wikimedia.org/wiki/File:0322_DNA_Nucleotides_Numbered.jpg).
```

:::{warning} TODO
Substitute panel (a) with Fig 5.2 from Schlick's book
:::

In the classic DNA double helix (or *duplex*) described by [Watson and Crick](doi:10.1038/171737a0), a flexible ladder-like structure is formed by two chains (strands) that wrap around a virtual central axis. The two chains' backbones consist of alternating sugar and phosphate units, while the rungs that connect them consist of nitrogenous bases held together by hydrogen bonds and stabilised by stacking interactions. A schematic of a single nucleotide and of the Watson-Crick mechanism are shown in [](#fig:DNA_schematics).

:::{seealso} The discovery of the double helix

```{figure} figures/photo_51.jpg
:name: fig:photo_51
:align: center

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

## Hybridisation and denaturation

Two nucleotides coming together and binding to each other for a base pair (BP). Since hydrogen bonds are [highly directional](#sec:hydrogen-bonds), HB formation requires that the two nucleotides approach each other with the right mutual orientation. However, note that bases on their own (*i.e.* not part of a strand) can pair with any other base, as well as with themselves, often with more than one arrangement ([](doi:10.1016/S0079-6603(08)60565-6)). It is only in an extended paired region that the correct orientation is obtained only if the two nucleotides are compatible (*i.e.* if they can form a canonical or a non-canonical base pair, as mentioned above), and the strands run anti-parallel to each other. In this geometry, which is adopted by the most common DNA and RNA helical conformations,

* the bases on one strand can align perfectly with their complementary bases on the opposite strand, allowing stable hydrogen bonding;
* the sugar-phosphate backbones run in opposite directions, allowing the bases to stack neatly on top of each other, contributing to the overall stability of the double helix through van der Waals forces and hydrophobic interactions;
* negatively charged phosphate groups on the sugar-phosphate backbone are positioned in a way that minimizes repulsion.

Interestingly the spontaneous process through which single strands pair and form a duplex, which is called *hybridisation*, is now taken for granted, but was initially met with skepticism[^without_enzyme], since researchers at the time believed that duplex formation required enzymes to overcome electrostatic and entropic penalties ([](doi:10.1017/S0033583509004776)).

In addition to hydrogen bonding, duplexes are stabilised by interactions acting between the aromatic rings of consecutive nitrogenous bases. These so-called base stacking interactions have a hydrophobic nature that is complemented by van der Waals and dipole-dipole termsthat contribute significantly to the stability and structural integrity of the DNA molecule. Rather counterintuitively, there is now ample evidence that base stacking is at least as important for duplex stability as base pairing, if not more (see *e.g.* [](doi:10.1093/nar/gkj454) and [](doi:10.1038/s41565-023-01485-1)).

:::{hint} Measuring the stacking strength

The classic study of [](doi:10.1016/j.jmb.2004.07.075) used a clever experimental setup based on the observation that a duplex where one of the composing strands is broken in two (*i.e.* it is *nicked*) has a mobility that is lower compared to an intact one. The reason is that a nicked strand can take two conformations: one straight, where all the stacking interactions are present that resembles that of an intact duplex, and one kinked, where all base pairs are formed but the stacking interaction at the nicked site is broken.

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

Assuming that each base can either be unbound or involved in a single base pair[^no_triplex], the secondary structure of nucleic acids is defined by the pairing information of the nucleotides. By assigning a unique index to each nucleotide, this information can be represented as a list of index pairs, where each entry represents a single base pair. While, as I mentioned already, there are key differences between DNA and RNA, the surge of DNA nanotechnology have blurred the boundaries between the two, and many of the secondary structures that were unique to RNA can also be found in synthetic DNA systems. As a result, in this part I will refer to generic "nucleic acids", if not explicitly stated otherwise.

[^no_triplex]: This is already an approximation, since some non-canonical (*e.g.* Hoogsteen) base pairings can connect a nucleotide to another which is already involved in a base pair.

## Nearest-neighbour models


(sec:canonical_helices)=
# Canonical helices

When two complementary strands pair together, they wrap around each other and form a periodic helical structure. The type of the resulting helix depends, in general, on the strand type (DNA or RNA), external conditions (ionic strength, pH, water concentration, *etc.*), and sequence.

## B-DNA

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

## A-DNA and A-RNA

```{figure} figures/adna.png
:name: fig:adna
:align: center
:width: 500px

An A-DNA double helix composed of 12 base pairs represented as in [](#fig:bdna). The structure has been generated with the [3DNA 2.0 webserver](http://web.x3dna.org).
```

DNA turns into its A form when dyhdrated, which can happen as a result of human intervention (*e.g.* when crystalling samples to study with X-rays), but also under natural conditions, when water is scarce or specific proteins bind to the DNA, in place of some of the water molecules[^A-DNA]. A-DNA is a right-handed helix that is more compact than B-DNA: it has about 11 base pairs per turn and a diameter of approximately 2.3 nanometers. A-DNA features deeper major grooves and shallower minor grooves compared to B-DNA. The helical structure is more tightly wound, with the base pairs tilted by abouy $20^\circ$ relative to the helix axis, and the distance between base pairs is shorter. Moreover, the sugars are in the C3'-endo conformation rather than in the C2'-endo conformation like in B-DNA. A-DNA is less prevalent in biological systems but can be found in certain DNA-RNA hybrids. Additionally, double-stranded RNA molecules or regions typically adopt the A-form helical structure.

[](#fig:adna) shows a perfect A-DNA helix made of the same number of base pairs as the one shown in [](#fig:bdna). Compared to B-DNA, it is evident that A-DNA is shorter and thicker, and the base pairs are inclined by $\approx 20^\circ$ with respect to the helical axis.

## Z-DNA

```{figure} figures/zdna.png
:name: fig:zdna
:align: center
:width: 600px

A Z-DNA double helix composed of 12 GC base steps (*i.e.* 24 base pairs) represented as in [](#fig:bdna). Note the zig-zag pattern drawn by the backbone, which gives Z-DNA its name. The structure has been generated with the [3DNA 2.0 webserver](http://web.x3dna.org).
```

Another important DNA conformation is Z-DNA, which is a less common and structurally distinct form of DNA characterized by a left-handed double helical structure, in contrast to the right-handed helices of A-DNA and B-DNA. It has a zigzag backbone, which gives Z-DNA its name, and consists of about 12 base pairs per turn with a diameter of approximately 1.8 nanometers. The major and minor grooves are less distinct in Z-DNA, and the structure is more elongated and thinner. Z-DNA requires a specific sequence of alternating purines and pyrimidines (*e.g.* repeating GC steps) and can form under high salt conditions or negative supercoiling, *i.e.* twisting of the DNA in the left-handed direction. Though less prevalent than B-DNA, Z-DNA is thought to have important biological functions in gene regulation, genetic recombination, and in the context of certain protein-DNA interactions.

[](#fig:zdna) shows a perfect Z-DNA helix made of 24 base pairs (12 repeating GC steps) which shows the main features of Z-DNA: the left handedness of the helix, the zig-zag pattern that gives this conformation its name, and the small BP inclination.

[^A-DNA]: In both cases the A form [is thought to protect the genetic material from extreme conditions](https://en.wikipedia.org/wiki/A-DNA#Biological_function).

:::{table} The main (average) properties of the three canonical forms of DNA and RNA.
:label: tbl:helices
:align: center

| Form | Handedness | BPs per turn | Diameter ($\AA$) | Pitch ($\AA$) | BP inclination | Sugar pucker
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| B-DNA | Right | 10 | 20 | 34 | 0$^\circ$ | C2'-endo |
| A-DNA/A-RNA | Right | 11 | 26 | 28 | 20$^\circ$ | C3'-endo
| Z-DNA | Left | 12 | 18 | 45 | -7$^\circ$ | alternating[^Z-DNA_pucker]
:::

A comparison between the average properties of the three main helical conformations is shown in [](#tbl:helices).

[^Z-DNA_pucker]: The alternating pyrimidines and purines take the C2'-endo and C3'-endo conformations, respectively.

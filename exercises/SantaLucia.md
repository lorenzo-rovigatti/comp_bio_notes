---
title: The SantaLucia model
---

:::{attention}
Remember to read the [homeworks](#sec:homeworks) section to understand what you need to submit for this assignment.
:::

# The main assignment

The SantaLucia (SL) model is the most famous (and widely-used) [nearest-neighbour model](#sec:NN_models) for DNA. The mandatory part of this assignment is to write two codes, one that computes the melting temperature of a double strand given its sequence, and another that computes the melting temperature of the secondary structure of a single strand.

The following instructions apply to the first code:

1. The code takes as input the two strands that make up the double strand and check that they are the same length and fully complementary. **Nota Bene:** this check is not required if you decide to also include the terms that take care of mismatches.
2. You can assume that the total strand concentration is $C^\circ = 1$ M for simplicity.

The following instructions apply to the second code:

1. The code takes as input the strand sequence and a secondary structure, specified with the [dot-paren notation](#sec:dot-paren).
2. The code should check that the secondary structure is compatible with the given sequence (*e.g.* check that they are of the same length, that the base pairs specified in the secondary structure are valid Watson-Crick pairs, and that each opening parenthesis has a closing partner).

The two codes should

1. Evaluate the enthalpic and entropic contributions to the given secondary structure $\Delta H^\circ$ and $\Delta S^\circ$, according to the SantaLucia model.
2. Apply the two-state model to plot the melting temperature of the secondary structure. **Nota Bene:** for complicated examples the main assumption behind the two-state model[^two_state] will not hold, but for the sake of this assignment we will pretend that this is never the case.

:::{important} Choose and document the supported features
A fully-fledged NN model has many parameters that take into account all possible secondary structures (internal mismatches, internal loops, bulges, hairpins, coaxial stacking, *etc.*). You do not have to necessarily to write a code that supports all these motifs, but:

1. the more the better (although supporting only the basic ones, shown in [](#tbl:SL_watson_crick), would be already good enough);
2. if you choose to support only a subset of structures, the documentation should be clear about it, and the code should explicitly fail when an unsupported secondary structure is supplied.
:::

[^two_state]: That there are no stable or metastable intermediate states.

# Possible extensions

* For the double-strand code, add the possibility of choosing the strand concentrations.
* Add the possibility of visualising the melting curve of the given secondary structure (*i.e.* the yield of the secondary structure as a function of temperature).
* Add the possibility of choosing the NN model: you can, for instance, add support for one of the RNA models listed [here](https://rna.urmc.rochester.edu/NNDB/).
* Add all the NN terms.

# Additional details

The table below contain the thermodynamic parameters of the model for the different base stacks, together with the initiation, terminal AT, and symmetry correction contributions. You will find the others, together with the explanation of how to use them, in [](doi:10.1146/annurev.biophys.32.110601.141800).

:::{table} SL parameters for Watson-Crick pairs at 1 M
:label: tbl:SL_watson_crick
:align: center

|Term|$\Delta H^\circ$ (kcal / mol)|$\Delta S^\circ$ (cal / mol K)|
|:---|:---:|:---:|
AA/TT | −7.6 | −21.3 |
AT/TA | −7.2 | −20.4 |
TA/AT | −7.2 | −21.3 |
CA/GT | −8.5 | −22.7 |
GT/CA | −8.4 | −22.4 |
CT/GA | −7.8 | −21.0 |
GA/CT | −8.2 | −22.2 |
CG/GC | −10.6 | −27.2 |
GC/CG | −9.8 | −24.4 |
GG/CC | −8.0 | −19.9 |
Initiation | +0.2 | −5.7 |
Terminal AT penalty | +2.2 | +6.9 |
Symmetry correction | 0.0 | −1.4 |
:::

For the hairpin code, you can model the free-energy cost of a loop of length $L$ as

$$
F_H(L) = 4.8 + 0.19 L - 0.0031 L^2,
$$

where the resulting number is expressed in kcal / mol.

Just to give you some numbers you can use to test your code, consider the following secondary structure (which is a duplex made of 6 base pairs):

```
5'-CGTTGA-3'
3'-GCAACT-5'
```

Its SL free energy contributions are $\Delta H = -409$ kcal/mol, $\Delta S = -114.6$ cal/mol K, so that $\Delta G = -5.35$ kcal/mol at 37$^\circ$. Its melting temperature at $C = 1$ M is $T_m = 79.5$ C.

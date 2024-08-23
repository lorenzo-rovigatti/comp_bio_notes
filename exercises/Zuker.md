---
title: The Zuker's algorithm
---

:::{attention}
Remember to read the [homeworks](#sec:homeworks) section to understand what you need to submit for this assignment.
::: 

# The main assignment

The [Zuker's algorithm](#sec:zuker) is the state-of-the-art deterministic method to obtain the MFE secondary structure of an RNA strand. However, it is only as good as the parameters used in its implementation. The aim of this exercise is to take [my implementation](../notebooks/RNA_folding.ipynb) (or your own, if you feel brave enough) and extend it so that it uses at least some of the real [nearest-neighbour model](#sec:NN_models) parameters. Since we are talking about RNA, the experimental values should be taken from [](doi:10.1093/nar/gkl472), which can also be found [here](https://rna.urmc.rochester.edu/NNDB/), in the "RNA (Turner 2004)" section.

In order to be positively evaluated, your code should at least include the proper free-energy contribution of a "perfect" dinucleotide stack, which is defined as two base pairs formed by nucleotides that are complementary. For instance, "AA/UU" or "AG/UC" are perfect stacks, while "AU/GC" is not.

Remember that NN models given the free energy contribution in terms of $\Delta H$ and $\Delta S$, or $\Delta G$ and $\Delta H$ (whence $\Delta S$ can be computed). This splitting makes it possible to calculate the $\Delta G = \Delta H - T \Delta S$ at any given temperature. This means that, once you have added the "real" stacking contribution, the MFE structure will somehow depend on temperature. Once you think your code works, prepare a plot of $N_{\rm bp}$ vs $T$ for an RNA strand, where $N_{\rm bp}$ is the number of base pairs in the MFE structure. You can use the following [tRNA sequence](https://www.rcsb.org/structure/1I9V): `GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA`.

# Possible extensions

* Add the possibility of obtaining the MFE structure of a subsequence of the target RNA strand, $S_{ij}$. Since with the Zuker's algorithm the MFE structures of *all* subsequences are computed during the fill-in phase, you will have to change the traceback function so that it acts only on the subsequence of interest.
* Add more parameters from Turner's model. For instance, you could use the correct values for hairpins, interior and bulge loops instead of the functions I fitted to the model. You can find all the parameters in "Python" format [here](https://github.com/Lattice-Automation/seqfold/blob/main/seqfold/rna.py), but be careful: beyond 30, the energy contribution of the loops is not tabulated. However, the $\Delta G$ of loops of size $N$ can be extrapolated *via* this relation:
$$
\Delta G_N = \Delta G_{30} + 2.44 R T \log(N / 30),
$$
where $R$ is the gas constant, $\Delta G_{30}$ is the free energy penalty of a loop of length $30$.

(sec:zuker-additional-details)=
# Additional details

For your convenience, here is a Python snippet that shows how to compute the $\Delta G$ of a "perfect stack":

```python
T = 273.15 + 37 # fixed temperature of 37° C

# Delta H and Delta S contributions to base steps (from Turner 2004)
dH_dS = {
    "AA/UU": (-6.8, -19.0),
    "AC/UG": (-11.4, -29.7),
    "AG/UC": (-10.5, -27.1),
    "AU/UA": (-9.4, -26.8),
    "CA/GU": (-10.4, -26.8),
    "CC/GG": (-13.4, -32.6),
    "CG/GC": (-10.6, -26.4),
    "CU/GA": (-10.5, -27.1),
    "GA/CU": (-12.4, -32.2),
    "GC/CG": (-14.9, -37.1),
    "GG/CC": (-13.4, -32.6),
    "GU/CA": (-11.4, -29.7),
    "UA/AU": (-7.7, -20.6),
    "UC/AG": (-12.4, -32.2),
    "UG/AC": (-10.4, -26.8),
    "UU/AA": (-6.8, -19.0),
    # wobble pairs
    "AG/UU": (-3.2, -8.4),
    "AU/UG": (-8.8, -23.9),
    "CG/GU": (-5.6, -13.5),
    "CU/GG": (-12.1, -32.2),
    "GA/UU": (-12.8, -37.1),
    "GC/UG": (-12.6, -32.6),
    "GG/CU": (-8.3, -21.9),
    "GG/UC": (-12.1, -32.2),
    "GG/UU": (-13.5, -41.9),
    "GU/CG": (-12.6, -32.6),
    "GU/UA": (-8.8, -23.9),
    "GU/UG": (-14.6, -51.3),
    "UA/GU": (-7.0, -19.3),
    "UC/GG": (-8.3, -21.9),
    "UG/AU": (-7.0, -19.3),
    "UG/GC": (-5.6, -13.5),
    "UG/GU": (-9.3, -31.0),
    "UU/AG": (-12.8, -37.1),
    "UU/GA": (-3.2, -8.4),
    "UU/GG": (-13.5, -41.9),
}

def dG(base_step):
    return dH_dS[base_step][0] - T * (dH_dS[base_step][1] / 1000.0)
    
print("The delta G associated to the 'UU/GA' stack is", dG("UU/GA"))
```

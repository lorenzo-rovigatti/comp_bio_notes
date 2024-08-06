---
title: Folding
---

# Proteins

```{tip}
The main references for this part are @finkelstein2002protein and @bialek2012biophysics. Note that here I refer mostly to globular proteins.
```

Proteins can be denatured by changing the solution conditions. Very high or low temperature and pH, or high concentration of a denaturant (like urea or guanine dihydrochloride) can "unfold" a protein, which loses its solid-like, native structure, as well as its ability to perform its biological function.

Investigations of the unfolding of small globular proteins showed that, as the denaturing agent (*e.g.* temperature or denaturant concentration) increases, many of the properties of the protein go through an "S-shaped" change, which is characteristic of cooperative transitions.

Furthermore, calorimetric studies show that the denaturation transition is an "all-or-none" transition, *i.e.* that the "melting unit" of the transition is the protein as a whole, and not some subparts. This of course applies to single-domain small proteins, or to the single domains of larger proteins. Here "all-or-none" means that the protein exist only in one of two states (native or denaturated), with all the other "intermediate" states being essentially unpopulated at equilibrium. Therefore, this transition is the microscopic equivalent of a first-order phase transition (*e.g.* boiling or melting) occurring in a macroscopic system. Of course, since proteins are finite systems and therefore very far from the thermodynamic limit, this is not a true phase transition, as the energy jump is continuous, and the transition width is finite.

:::{note} The van't Hoff criterion
Let $N$ and $D$ be two states in chemical equilibrium, *i.e.* $N \rightleftharpoons D$, at a certain temperature $T$. Recalling the two-state equilibrium concepts derived [before](#sec:two_state), and defining $p(T)$ as the fraction of, say, state $D$, we have

$$
K_{\rm eff}(T) = \frac{p(T)}{1 - p(T)} = e^{- \beta \Delta G_{ND}},
$$

where $K(T)$ is the effective equilibrium constant at temperature $T$ and $\Delta G_{ND} = \Delta H_{ND} - T \Delta S_{ND}$ is the free-energy differences between the two states. By simple algebra we have that

$$
\Delta H_{ND} = k_B T^2 \frac{d \log K_{\rm eff}(T)}{dT}.
$$

Note that $K_{\rm eff}(T)$ is connected to $p(T)$, which can be measured experimentally (for instance with Circular Dicroism). The value of $\Delta H_{ND}$ obtained can be interpreted as the heat consumed by the "melting unit" associated to the transition.

The enthalpy change can also be estimated by integrating the heat capacity, measured with calorimetry, over the range of temperature associated to the transition. This is the heat consumed by all proteins during the transition, $\Delta H_{\rm tot} = N_p \Delta H_{\rm cal}$, where $N_p = m / M$ is the number of proteins, with $m$ being the total mass and $M$ the protein's molecular mass.

If $H_{\rm cal} = \Delta H_{ND}$, it means that the "melting unit" is the whole protein, and therefore the transition is of the "all-or-none" type, as it is the case of small globular proteins.

Note that there are [some subtleties](10.1110/ps.8.5.1064) that should be taken into account when applying this method, but the general concept is still useful, provided that the results are confirmed by other experimental methods.
:::

```{figure} figures/molten_globule.png
:name: fig:molten_globule
:align: center

(a) Phase diagram of the conformational states of lysozyme at pH 1.7 as a function of denaturant (guanine dihydrochloride) concentration and temperature. The red lines correspond to the mid transitio, the dashed lines outline the transition zone. (b) Schematic model of the native and molten globule protein states. Here the protein consists of only two helices connected by a look, and the side chains are shown as shaded regions. Adapted from @finkelstein2002protein.
```

But how does the denatured state look like? As discovered by using a plethora of experimental methods, which often seemed contradicting each other, the answer depends on the denaturing conditions. [](#fig:molten_globule)(a) shows the phase diagram of lysozyme (a globular small protein) for different temperatures and concentrations of a denaturing agent. Increasing the latter leads to a transition to a disordered coil state where essentially no secondary structure is present. By contrast, the effect of temperature is qualitatively different, as the protein partially melts in a state (the *molten globule*, MG) that retains most of its secondary structure, but it is not solid like and it cannot perform any biological function. A schematic of the difference between the native and molten globule states is shown in [](#fig:molten_globule)(b).

The cooperative transition that leads to the molten globule is due to the high packing of the side chains, which is energetically favourable but entropically penalised, since it impedes the small-scale backbone fluctuations that trap many of the internal degrees of freedom of the protein. The liberation of these small-scale fluctations, and the resulting entropy gain, requires a slight degree of swelling, which is why the MG is not much larger than the native state. However, such a swelling is large enough to greatly decrease the van der Waals attractions, which are strongly dependent on the distance, and to let solvent (water) molecules in the core.

In general, not all proteins behave like this, as some (usually small) proteins that unfold directly into a coil, others that form the molten globule under the effect of specific denaturants, and coils under the effect of others, *etc.* However, the unfolding transition has nearly always an "all-or-none", cooperative nature.

:::{note} Coil-globule transition in polymers
The cooperative denaturing of a protein is *very* different from the coil-globule transition observed in polymers, even when no molten globule is involved and the protein turns directly into a coil. The reason is that the coil-globule transition is not a microscopic analog of a first-order phase transition, and cannot be described by two-state models.
:::

Generally speaking, the denaturation of a protein is a reversible process, provided that some experimental conditions are met (the protein is not too large, its chemical structure has not been altered by the denaturation process or by post-translational modifications, *etc.*). The fact that protein folding and unfolding is a reversible process implies that all the information required to build a protein is stored in its sequence, and that its native structure is thermodynamically stable. Therefore, the folding process can be, in principle, understood by using thermodynamics and statistical mechanics. The sequence $\to$ structure hypothesis is sometimes called [Anfisen's dogma](https://en.wikipedia.org/wiki/Anfinsen%27s_dogma), named after the Nobel Prize laureate Christian B. Anfisen, whose studies on refolding of ribonuclease were key to establish the link between "the amino acid sequence and the biologically active conformation" ([](doi:10.1073/pnas.47.9.1309), [](doi:10.1126/science.181.4096.223)).

:::{warning}
There are many "exceptions" to this dogma, mostly because large proteins need help to fold, because some native structures require post-translational modifications, and because *in vitro* renaturation is difficult and can be hindered by inter-molecular or (for large proteins) even intra-molecular aggregation. Nevertheless, thanks to the validity of this "thermodynamic hypothesis", many common features of denaturation and folding can be understood by means of rather simple arguments which, as I will show, can be used to build simple models.
:::

```{figure} figures/energy_gap.png
:name: fig:energy_gap
:align: center

(a) Typical curves for the energy, entropy and free-energy of a protein as functions of the (intramolecular) density. The D and N letters mark the positions of the denatured and native states, respectively. The hash pound marks the position of the maximum of the free-energy barrier. (b) A sketch of the energy spectra (bottom panels) and associated $S(E)$ curves (top panels) for (A) a protein and (B) a random heteropolymer. Adapted from @finkelstein2002protein.
```

[](#fig:energy_gap)(a) shows the origin of the "all-or-none" transition observed in protein denaturation in thermodynamic terms. Here I will use the description given by @finkelstein2002protein (whence the figure comes), almost word by word. The panel shows the energy $E$ and the entropy $S$ of a protein as a function of its internal density (in some arbitrary units). The energy is at its minimum at close packing, while the entropy S increases with decreasing density. Initially, this increase is slow, since the side chains can only vibrate, but their rotational isomerization is impossible. When the density decreases enough that the latter becomes possible the increase is much steeper. Finally, when free isomerization has been reached the entropy growth becomes slow again.

The nonuniform growth of the entropy is due to the following. The isomerization needs some vacant volume near each side group, and this volume must be at least as large as a $CH_3$ group. And, since the side chain is attached to the rigid backbone, this vacant volume cannot appear near one side chain only but has to appear near many side chains at once. Therefore, there exists a threshold (barrier) value of the globule's density after which the entropy starts to grow rapidly with increasing volume (decreasing density) of the globule. This density is rather high, about $80\%$ of the density of the native globule, since the volume of the $CH_3$ group amounts to $\approx 20\%$ of the volume of an amino acid residue. Therefore, the dependence of the resulting free energy $F = E - TS$ on the globule's density $\rho$ has a maximum (the so-called free-energy barrier) separating the native, tightly packed globule (N) from the denatured state (D) of lower density. It is because of the presence of this free energy barrier that protein denaturation occurs as an "all-or-none" transition akin to melting or nucleation, independent of the final state of the denatured protein.

It is useful, when talking about proteins, to also introduce the concept of random heteropolymers, which are polymers with sequences that have not been selected by nature, and therefore have many (degenerate or quasi-degenerate) ground states rather than a single native conformation. [](#fig:energy_gap)(b) shows the qualitative difference between a protein and a random heteropolymer. The bottom panels of the figure show the typical energy spectra of the two chains, where each line corresponds to a single conformation. At high energy both spectra are essentially continuous. However, the low-energy, discrete parts look very different: the energy difference between the ground state and the continuous part of the spectrum, where the majority of the "misfolded" states lie, is much larger in a protein than in a random heteropolymer. It is this energy gap that creates the free-energy barrier of [](#fig:energy_gap)(a).

We can also define the entropy $S(E)$, which is proportional to the logarithm of the number of structures having energy $E$, and the temperature corresponding to the energy $E$, implicitly defined by the relation $d S(E) / dE = 1 / T(E)$. As sketched in [](#fig:energy_gap)(b), the peculiar shape of the protein energy spectrum gives raise to a concave part of $S(E)$, which means that there is a temperature at which the lowest-energy, native conformation coexists with many high-energy (denatured) ones. In this simple picture, this is the melting temperature $T_M$. This temperature should be compared with the so-called "glass temperature" $T_G$, called $T_C$ in the figure, which is the temperature associated to the lowest energy of the continuous band of the spectrum. Above $T_G$ the chain behaves, from the kinetic point of view, like a viscous liquid, while below the glass temperature the behaviour is glassy-like (*i.e.* non-Arrenhius dependence on temperature, aging, non ergodicity, *etc.*).

If $T_M < T_G$, the slow dynamics asssociated to the glass state makes it kinetically impossible to reach the ground state, and the chain remains trapped in higher-energy misfolded states. By contrast, if $T_M > T_G$, folding is possible, and results obtained with minimalistic models show that the fraction $T_M / T_G$ is a good proxy for "foldability": the larger this fraction, the faster the protein folds.

```{figure} figures/frustration.png
:name: fig:frustration
:align: center
:width: 600

(a) An example of frustrated interactions: if interaction is such that only like colours want to be close to each other, there is no way of arranging the three spheres so that all favourable contacts are made. As a result, the ground state is degenerate. (b) A schematic energy landscape of a 60 amino-acid helical protein. $A$ and $Q$ are the fractions of correct dihedral angles in the backbone and correct native-like contacts, respectively, $\Delta E$ is the "ruggedness" of the landscape, and $\delta E_s$ is the energy gap. Taken from [](doi:10.1073/pnas.92.8.3626).
```

All these concepts are not merely qualitative, but have been grounded in theory by using sophisticated statistical mechanics approaches that have generated the "funnel landscape" folding picture, which is schematically shown in [](#fig:frustration)(b).

The basic idea is to leverage the concept of frustration, which happens when the interactions between the different parts of a system are such that they cannot be satisfied, from the energetic point of view, all at the same time. If this is the case, then there is no single ground state, but rather many low-energy stable states that, in a many-body system, are uncorrelated from each other and separated by (possibly large) energy barriers. An example is provided in [](#fig:frustration)(a). The resulting "energy surface" (or landscape) is said to be rugged (or rough), and generates a dynamics where the system sits in a valley for a long time before being able to jump over the barrier and move to a different stable conformation. This is a hallmark of glassiness.

To draw a parallel with proteins, we can consider a random heteropolymers, where the interactions among the different amino acids will be frustrated, blocking the system from finding a single well isolated folded structure of minimum energy. A candidate principle for selecting functional sequences is thus the minimization of this frustration. Of course even in the absence of frustration, there are still energetic barriers on the path towards the native conformation due to local structural rearrangements which still give raise to a rugged landscape, slowing down the kinetics towards the folded state. This scenario has come to be called a folding "funnel", emphasizing that there is a single dominant valley in the energy landscape, into which all initial configurations of the system will be drawn.

A classic funnel diagram, such as the one shown in [](#fig:frustration)(b), represents an energy-entropy landscape, with the width being a measure of the entropy, whereas the depth represents both the energy and two correlated order parameters $Q$ and $A$, which are the fractions of native contacts and correct dihedral angles in the protein backbone. Although in reality the landscape is multidimensional, here the projection attempts to retain its main features, such as the barrier heights, which are a measure of the "ruggedness" or "roughness" of the landscape. The typical height of these barriers, together with the energy gap and the number of available conformations, are the three main parameters of the Random Energy Model, which is one of the main theoretical tools in this context.

## Some simple lattice models

### The helix-to-coil transition

As we discussed, polypeptides can, in general, be reversibly denatured (with temperature or pH, for instance). Here we theoretically investigate this phenomenon in the context of secondary structure with the [Zimm-Bragg model](doi:10.1063/1.1730390), which is a toy model for a chain that can exhibit a continuous coil-to-helix transition. Note that it is possible to synthesise polypeptide chains that undergo such a transition, see *e.g.* [](doi:10.1021/ja01583a070).

Consider a polymer chain composed of $N$ residues. Each residue can be in one of two states, "coil" (c) or "helix" (h), so that a microscopic state of the chain can be expressed as a one-dimensional sequence, like `ccchhhhhcchhhchcc`. In equilibrium, the probability (or statistical weight) to observe a given chain microstate $i$ is just the Boltzmann factor, $e^{-\beta \Delta F_i}$, where $F_i$ is the free-energy cost of the microstate, so that the total (configurational) partition function is

$$
Q = \sum_{i \in {\rm config}} e^{-\beta \Delta F_i}.
$$

Since (free) energies are always defined up to a constant, we define the free-energy of a coil of any length to be zero by definition, so that the statistical weight of a coil residue is always 1. By contrast, forming a sufficiently long helical segment should be advantageous, otherwise no transition would occur. However, we know that in [the helices](#sec:protein_helices) formed by polypeptides, each $i$-th amino acid is hydrogen-bonded to the $(i + k)$-th, where $k = 4$ for a $\alpha$-helix. As a consequence, since one full helical turn needs to be immobilised before seeing any benefit from the formation of hydrogen bonds, the formation of a helix incurs a large entropic penalty that needs to be paid every time a coil segment turns into a helical one, $f_{\rm init}$. However, once a helix is formed and the initial price has been paid[^nucleus], adding an H contributes a fixed amount $f_h$ to the helix free energy. As a result, the free-energy contribution of a helical stretch of size $n$ is $\Delta F = f_{\rm init} + n f_h$, so that is statistical weight is

$$
\exp(-\beta \Delta F) = \exp(-\beta f_{\rm init}) \left[ \exp(-\beta f_h) \right]^n \equiv \sigma s^n,
$$

where we have defined the cooperativity (or initiation) parameter, $\sigma \equiv \exp(-\beta f_{\rm init})$, and the helix-elongation (or hydrogen-bonding) parameter, $s \equiv \exp(-\beta f_h)$. Note that, since $f_{\rm init}$ is purely entropic, $\beta f_{\rm init} = -S_{\rm init} / k_B$ and therefore $\sigma$ does not depend on temperature.

:::{warning}
Strictly speaking, it is not correct to add the $n f_h$ term is $n < 4$, since at least one turn ($\approx 4$ amino acids) should be present in order to stabilise an $\alpha$-helix. There are some more complicated models that take this into account (see *e.g.* the [Lifson-Roig model](https://en.wikipedia.org/wiki/Lifson%E2%80%93Roig_model) or the [Generalized Model of Polypeptide Chain](doi:10.1140/epje/i2013-13046-7)), but here we will keep it simple and study the simplest model that exhibits the main features of the transition.
:::

Following @finkelstein2016protein, we define two accessory quantities: $C_i$ is the partition function of a chain of size $i$, where the last monomer is in the coil state, and $H_i$ is the partition function of a chain of size $i$, where the last monomer is in the helix state. With this definition, the partition function of the chain is $Q_i = C_i + H_i$. We can explicitly compute the first terms:

* For $i = 1$, the sequence is either `c` or `h`, so that $C_1 = 1$ and $H_1 = \sigma s$.
* For $i = 2$, we have two sequences ending with c, which are `hc` and `cc`, and two sequences ending with h, which are `hh` and `ch`. Each partition function is given by a sum of the statistical weights of the allowed microstates, so that we have $C_2 = 1 + \sigma s$ and $H_2 = \sigma s^2 + \sigma s$.
* For $i = 3$, there are 8 available sequences, which are built by taking each $i = 2$ sequence and adding either `c` or `h` at its end. Each of the resulting sequence will have a statistical weight that is that of the base $i = 2$ sequence, multiplied by 1 if the sequence ends with `c`, and by $s$ or $\sigma s$ if the sequence ends with `h` and the preceeding character is `h` or `c`, respectively. Summing up the contribution we obtain $C_3 = 1 + 2 \sigma s$ and $H_3 = \sigma s + \sigma s^2 + \sigma^2 s^2 + \sigma s^3$.

In general the operations carried out to compute the $i = 3$ partial partition functions can be applied to any other value of $i$, meaning that we can write down recursive relationships connecting $H_{i}$ and $C_{i}$ to $H_{i-1}$ and $C_{i-1}$. These relationships take a very simple form if expressed in matricial form:

$$
(C_i, H_i) =
(C_{i-1}, H_{i-1})
\begin{pmatrix}
1 & \sigma s\\
1 & s
\end{pmatrix}.
$$ (eq:zimm_bragg_recursive)

It is convenient to define the matrix

$$
\matr{M} \equiv \begin{pmatrix}
1 & \sigma s\\
1 & s
\end{pmatrix}
$$

so that the the total partition function of a chain of $N$ monomers can be written as[^zimm_bragg_difference]

$$
Q_N = (1, 0) \matr{M}^N
\begin{pmatrix}
1\\1
\end{pmatrix},
$$ (eq:zimm_bragg_Q)

where the first vector-matrix multiplication on the right-hand side gives $(C_N, H_N)$, which multiplied by the $(1, 1)$ column vector yields $C_N + H_N$. The partition function can be written in terms of the eigenvalues of $\matr M$ if we note that 

$$
\matr M^N = \matr M \matr M \ldots \matr M = \matr A \matr\Lambda \matr A^{-1} \matr A \matr\Lambda \matr A^{-1} \ldots \matr A \matr\Lambda \matr A^{-1} = \matr A \matr\Lambda^N \matr A^{-1},
$$

where $\matr \Lambda = \matr A^{-1} \matr M \matr A = \begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix}$, and $\lambda_1$ and $\lambda_2$ (with $\lambda_1 > \lambda_2$) are the eigenvalues of $\matr M$, which can be found by solving the secular equation, *viz.*:

$$
\det \begin{pmatrix}
1 - \lambda & \sigma s\\
1 & s - \lambda
\end{pmatrix} =
(1 - \lambda)(s - \lambda) - \sigma s = 0,
$$

whence

$$
\begin{align}
\lambda_1 &= 1 + \frac{s - 1}{2} + \sqrt{\left( \frac{s - 1}{2} \right)^2 + \sigma s}\\
\lambda_2 &= 1 + \frac{s - 1}{2} - \sqrt{\left( \frac{s - 1}{2} \right)^2 + \sigma s}.
\end{align}
$$

Note that the two eigenvalues are linked by the relation $\lambda_1 + \lambda_2 = 1 + s$. 

Since they can be of any length, the two eigenvectors $\vec{A}_1$ and $\vec{A}_2$ of a 2x2 matrix each have a single independent element, $a_k$ (where $k = 1, 2$), which is given by

$$
\begin{pmatrix}
1 - \lambda_k & \sigma s\\
1 & s - \lambda_k
\end{pmatrix}
\begin{pmatrix}
a_k\\
1
\end{pmatrix} = 0.
$$

We obtain $a_1 = \lambda_1 - s = 1 - \lambda_2$ and $a_2 = \lambda_2 - s = 1 - \lambda_1$, and therefore $\vec{A}_1 = \begin{pmatrix} 1 - \lambda_2\\1\end{pmatrix}$ and $\vec{A}_2 = \begin{pmatrix} 1 - \lambda_1\\1\end{pmatrix}$. The two eigenvectors are used to construct the diagonalising matrix and its inverse[^finkelstein_typo]

$$
\matr A = \begin{pmatrix}
1 - \lambda_2 & 1 - \lambda_1\\
1 & 1
\end{pmatrix},
\matr A^{-1} = \frac{1}{\lambda_1 - \lambda_2} \begin{pmatrix}
1 & \lambda_1 - 1\\
-1 & 1 - \lambda_2
\end{pmatrix}.
$$

We can now write Eq. [](#eq:zimm_bragg_Q) explicitly in terms of $\lambda_1$ and $\lambda_2$:

$$
\begin{align}
Q_N & = (1, 0) \matr{A}
\begin{pmatrix}
\lambda_1^N & 0\\
0 & \lambda_2^N
\end{pmatrix}
\matr{A}^{-1}
\begin{pmatrix}
1\\1
\end{pmatrix} =\\
&= \frac{1}{\lambda_1 - \lambda_2} \left[ (1 - \lambda_2)\lambda_1^{N + 1} - (1 - \lambda_1) \lambda_2^{N + 1} \right].
\end{align}
$$ (eq:zimm_bragg_Q_lambda)

From the partition function we can compute any observable we need. The most interesting quantity is the average fraction of residues that are in the helical state h, which can be compared to experiments. First of all, we note from the $C_i$ and $H_i$ we explicitly calculated and from Eq. [](#eq:zimm_bragg_recursive) that the partition function is a polynomial in $s$, and therefore can be written as 

$$
Q_N = \sum_k g_k s^k,
$$

where $g_k$ are coefficients that do not depend on $s$, but only on $\sigma$ and on the multiplicity of each state with $k$ helical residues. Therefore, by definition the probability that the chain has exactly $k$ helical residues is $P_N(k) = \frac{g_k s^k}{Q_N}$, which means that the average number of helical residues is

$$
\langle n_h \rangle = \frac{\sum_k k g_k s^k}{Q_N} = \frac{1}{Q_N} \sum_k s g_k \frac{d s^k}{ds} = \frac{s}{Q_N} \frac{d}{ds} \sum_k g_k s^k = \frac{s}{Q_N} \frac{\partial Q_N}{\partial s},
$$

since $\frac{ds^k}{ds} = k s^{k-1}$ and only $s^k$ depends on $s$. Defining the fraction of helical residues, $\theta_N \equiv \langle n_h \rangle / N$, and exploiting the known property of logarithms we find

$$
\theta_N = \frac{s}{N} \frac{\partial \log Q_N}{\partial s} = \frac{1}{N} \frac{\partial \log Q_N}{\partial \log s}.
$$ (eq:theta_N)

Substituting Eq. [](#eq:zimm_bragg_Q_lambda) it is possible to write $\theta_N$ in the following closed, albeit rather cumbersome, form :

$$
\theta_N = \frac{\lambda_1 - 1}{\lambda_1 - \lambda_2} \frac{1 + \left( \frac{\lambda_2}{\lambda_1} \right)^{N+1} - \frac{2}{N}\frac{\lambda_2}{\lambda_1 - \lambda_2} \left[ 1 - \left( \frac{\lambda_2}{\lambda_1} \right)^{N} \right]}{1 + \frac{\lambda_1 - 1}{1 - \lambda_2}\left( \frac{\lambda_2}{\lambda_1} \right)^{N+1}}
$$ (eq:zimm_bragg_theta)

For sufficiently long chains, $(\lambda_2 / \lambda_1)^N \ll 1$ and the expression can be greatly simplified:

$$
\theta_N \approx \frac{\lambda_1 - 1}{\lambda_1 - \lambda_2} \left( 1 - \frac{2}{N} \frac{\lambda_2}{\lambda_1 - \lambda_2} \right) \xrightarrow[N \to \infty]{} \frac{\lambda_1 - 1}{\lambda_1 - \lambda_2} = \frac{1}{2} + \frac{s - 1}{\sqrt{\left(\frac{s-1}{2}\right)^2 + \sigma s}}.
$$

The model as described and solved takes into account the cooperativity of the transition through the parameter $\sigma$. It is instructive to look at the non-cooperative solution, which can be obtained by setting $\sigma = $. In this case $\lambda_1 = 1 + s$ and $\lambda_2 = 0$, so that 

$$
\theta_N = \frac{s}{s + 1},
$$ (eq:zimm_bragg_non_coop)

which is independent of $N$ and describes the simple chemical equilibrium of a $c \rightleftharpoons h$ reaction. Indeed, defining the equilibrium constant $K_{hc} \equiv [h] / [c]$, where $[x]$ is the concentration of $x$, we find that the fraction of residues in the helical state is 

$$
\theta = \frac{[h]}{[h] + [c]} = \frac{K_{hc}}{K_{hc} + 1},
$$

which, if compared with Eq. [](#eq:zimm_bragg_non_coop), shows that $s = K_{hc}$.

```{figure} figures/zimm_bragg_fits.png
:name: fig:zimm_bragg_fits
:align: center
:width: 500

Comparison between experimental data (points) and theoretical fits to Eq. [](#eq:zimm_bragg_theta) (lines). The number of monomers is $N = 1500$ and $n = 26$ for the top and bottom curves, respectively. $s$ is fixed by the relation $\log s = 0.00614 (T - T_c)$, where $T_c$ is the temperature at which $\theta = 0.5$. Solid and dashed lines correspond to $\sigma = 2 \times 10^{-4}$ and $\sigma = 10^{-4}$, respectively. Adapted from [](doi:10.1063/1.1730390).
```

A comparison between the Zimm-Bragg theory and experiments on poly-$\gamma$-benzyl-L-glutamate of different lengths are reported in [](#fig:zimm_bragg_fits). It is clear that there is a strong $N$-dependence that shows the cooperative nature of the transition. Moreover, in fitting the curves it is assumed that $\sigma$ does not depend on $T$, since it has a purely entropic origin.

:::{note} The Zimm-Bragg *vs.* 1D Ising model
The Zimm-Bragg model is often said to be equivalent to the 1D Ising model, where the coupling between the spins is linked to the initiation free energy and the external field is linked to the hydrogen-bonding contribution. However, this is not correct, as it has been formally demonstrated that the Zimm-Bragg model is equivalent to a [Potts-like model](https://en.wikipedia.org/wiki/Potts_model), which is a generalisation of the Ising model that can be, in turn, mapped into a 1D Ising model with a temperature-dependent coupling (see [](doi:10.1103/PhysRevE.81.021921) for additional details).

However, the regular 1D Ising model can still be used to study the coil-helix transition. Check [the notes](https://www.roma1.infn.it/~sciortif/didattica/SOFTSTRUTTURA/SOFTSTRUTTURA/stretching-helixcoil.pdf) of the [Soft and Biological Matter](https://www.roma1.infn.it/~sciortif/didattica/SOFTSTRUTTURA/softstruttura.html) by Prof. Sciortino, where the model is solved with a similar method as the one used here (the [transfer matrix method](https://en.wikipedia.org/wiki/Transfer-matrix_method)).
:::

:::{seealso} Python implementation
Head over [here](./notebooks/zimm_bragg.ipynb) for a Python notebook where the average helicity is computed exactly with Eq. [](#eq:zimm_bragg_theta), and numerically by performing a simulation. The simulation makes use of the Monte Carlo Metropolis method, which is also introduced.
:::

[^nucleus]: In another context, that of first-order phase transition, we would say "once a nucleus forms".
[^zimm_bragg_difference]: Note that in the original Zimm-Bragg model the first three segments are always in the coil state, so that the exponent is $N - 3$ rather than $N$.
[^finkelstein_typo]: There is a typo in the Appendix B of @finkelstein2016protein: the $(2, 1)$ element of $A^{-1}$ should be $-1$ and not $1$.

(sec:HP_model)=
## The HP model

Globular proteins exhibit a unique characteristic compared to other isolated polymer molecules: their native structure is unique, meaning that the path formed by the backbone is essentially fixed, bar some small fluctuations. The question arises: what forces, determined by the amino acid sequence, contribute to this remarkable specificity? Of course, as discussed before, a significant constraint is provided by the fact that native structures are compact. However, even for compact chain molecules, many conformations are possible, with the number of maximally compact conformations increasing [exponentially with chain length](doi:10.1021/ma00202a031). Among the various physically accessible compact conformations, there are in principle many different interactions that can impact their thermodynamic stability and thus select the native structure. However, it is now well accepted that the interaction type that most contribute to this conformational reduction is hydrophobicity: given a specific amino acid sequence, the number of compact conformations that it can take is greatly reduced by the constraint that residues with hydrophobic side chains should reside inside the globule, while polar and hydrophilic residues should be on the protein's surface (see *e.g.* @finkelstein2002protein and [](doi:10.1021/bi00483a001)).

A simple model that has been used to demonstrate the importance of hydrophobic interactions in selecting the native structure is the [Hydrophobic-Polar (HP) protein lattice model](doi:10.1021/ma00200a030). In this model, proteins are abstracted as sequences composed solely of hydrophobic (H) and polar (P) amino acids. In the HP model, the protein sequence is mapped onto a grid, or lattice, which can be two-dimensional (2D) or three-dimensional (3D). Common lattice types used include square lattices for 2D models and cubic lattices for 3D models, but other choices are possible. 

Here I will present the 2D square lattice version of the model[^2D_HP]. Each of the $N$ amino acids of the protein is assigned a type: "H" means hydrophobic, and "P" means polar, so that the sequence of the protein is then the ordered list of amino acids (*e.g.* "HHPHPPPHHP"). The $i$-th amino acid in the sequence is identified by its index $i$ and occupies a point on this grid, and two amino acids cannot occupy the same grid point. Two amino acids $i$ and $j$ for which $|j - i| = 1$ are "connected neighbours" and share a backbone edge representing the peptide bond, which on the square lattice can be either horizontal or vertical. By contrast, two amino acids $i$ and $j$ for which $|i - j| > 1$ that are adjacent in space but not along the sequence are said to be "topological neighbours". Every pair of HH topological neighbours form a topological contact that contributes an amount of free energy $\epsilon = \equiv \beta \Delta F_{HH} < 0$, while the other contacts give $\epsilon_{HP}$ and $\epsilon_{PP}$. In the simplest version of the model $\epsilon_{HP} = \epsilon_{PP} = 0$. In this case, the total (dimensionless) free energy of a conformation is

$$
\beta E_i = m \epsilon,
$$

where $m$ is the number of HH topological contacts. Following [](doi:10.1021/ma00200a030), the lowest-energy (highest-$m$) conformations of a chain with a specific sequence and a length $N$ are called "native".

The concept of topological contacts is also useful to describe the degree of compactness of a conformation, defined as

$$
\rho \equiv \frac{t}{t^{\rm max}},
$$

where $t$ is the number of topological contacts, irrespective of the type of neighbouring residues, and $t^{\rm max}$ is the largest possible number of topological contacts for a chain of size $N$, which is

$$
t^{\rm max} = N + 1 - \frac{P_{\rm min}}{2},
$$

where $P_{\rm min}$ is the perimeter of the smallest box that can completely surround the protein. By this definition, (maximally) compact conformations have $\rho = 1.0$.

```{figure} figures/HP_model_example.png
:name: fig:HP_model_example
:align: center
:width: 600

The sequence "HPPHHPHHPHPHHH" has a single conformation with 7 topological contacts (shown on the left). All other conformations, like the two on the right, have fewer hydrophobic contacts. Hydrophobic and polar residues are shown as black and white circles, respectively. Taken from [](doi:10.1021/bi00483a001).
```

As an example, [](#fig:HP_model_example) shows the single lowest-energy conformation of a given sequence, together with two other compact conformations of higher free energy. It is rather evident that the native conformation has a core that is richer in hydrophobic residues compared to the other two.

In order to be more quantitative, some useful observables can be introduced. In general, the thermodynamic properties of a chain of length $N$ and given sequence can be computed from the partition function:

$$
Q = \sum_{i = 1}^{\Omega_N} e^{-E_i}
$$

where $\Omega_N$ is the number of distinct conformations of an $N$-chain, so that the sum runs over all the conformations. If the length of the chains is not too large[^HP_length], it is possible to numerically perform a complete (exhaustive) enumeration of all possible chain conformations, from which any quantity of interest can be computed as an ensemble average, *viz.*

$$
\langle A \rangle = \frac{1}{Q} \sum_{i}^{\Omega_N} A e^{E_i}.
$$

The partition function can also be written as a sum over conformations of fixed number of HH topological contacts $m$:

$$
Q = \sum_{m = 0}^{m^{\rm max}} g(m) e^{-m \epsilon},
$$

where $m^{\rm max}$ is the maximum number of HH topological contacts and $g(m)$ is the number of distinct chain conformations having $m$. Note that by this definition, the native conformations are those for which $m = m^{\rm max}$. We also define $G(m, u)$ as the number of conformations that have $m$ HH topological contacts and $u = t - m$ non-HH topological contacts, which is connected to $g(m)$ by

$$
g(m) = \sum_{u = 0}^{t^{\rm max} - m} G(m, u).
$$

Since partition functions are defined up to a multiplicative constant, we multiply $Q$ by $e^{m^{\rm max} \epsilon}$ to obtain[^HP_notation_abuse]

$$
Q = \sum_{m = 0}^{m^{\rm max}} g(m) e^{(m^{\rm max} - m) \epsilon}.
$$

With this definition it is straightforward to take averages over native conformations only. Indeed, if we take the $\epsilon \to -\infty$ limit, $e^{(m^{\rm max} - m) \epsilon}$ is non-zero only for $m = m^{\rm max}$, and the partition function becomes

$$
Q_{\infty} = g(m^{\rm max}).
$$

We now define some useful average values to characterise the different states. The average compactness over all conformations is

$$
\langle \rho \rangle = \frac{1}{Q} \sum_{m = 0}^{m^{\rm max}} \sum_{u = 0}^{t^{\rm max} - m} \frac{u + m}{t^{\rm max}} G(m, u) e^{(m^{\rm max} - m) \epsilon},
$$

and the average compactness over the native conformations is

$$
\langle \rho \rangle_{ns} = \frac{1}{Q_{\infty}} \sum_{u = 0}^{t^{\rm max} - m^{\rm max}} \frac{u + m^{\rm max}}{t^{\rm max}} G(m^{\rm max}, u) 
$$

where we used the fact that

$$
\lim_{\epsilon \to -\infty} \sum_{m=0}^{m^{\rm max}} G(m, u) e^{(m^{\rm max} - m) \epsilon} = G(m^{\rm max}, u)
$$

The average energy is proportional to the average number of HH contacts, which is

$$
\langle m \rangle = \frac{1}{Q} \sum_{m = 0}^{m^{\rm max}} m g(m) e^{(m^{\rm max} - m) \epsilon}
$$

Finally, we know already that most of the hydrophobic residues are buried in the core of globular proteins. In the HP model we define the core as the set of residues that are completely surrounded by other residues. Let the number of core residues be $n_i$, and the number of hydrophobic core residues be $n_{hi}$, then the degree of hydrophobicity of a given conformation can be estimated as

$$
x \equiv \frac{n_{hi}}{n_i},
$$

and its averages over all and native conformations can be computed as

$$
\begin{align}
\langle x \rangle & = \frac{1}{Q} \sum_{i = 0}^{\Omega_N} x e^{(m^{\rm max} - m_i) \epsilon}\\
\langle x \rangle_{ns} & = \frac{1}{g(m^{\rm max})} \sum_{i = 1}^{g(m^{\rm max})} x.
\end{align}
$$

```{figure} figures/HP_model_dill_1.png
:name: fig:HP_model_dill_1
:align: center

Ensemble averages of folding (top) and non-folding (bottom) sequences. Here $Z$ is the value of the partition function, which I call $Q$. All quantities are plotted as a function of $-\epsilon$. Adapted from [](doi:10.1021/ma00200a030).
```

The ensemble averages of these quantities, as well as of the partition function, which is a measure of the number of available conformations, as a function of $-\epsilon$ for two sequences are shown in [](#fig:HP_model_dill_1). Looking at the top panels it is clear that if the HH interaction is strong the preferred conformations of the `HPHPPHPPHH` sequence have low energy, are few ($Q \to 1$), compact, $\langle \rho \rangle = 1$, and have a purely hydrophobic core, $\langle x \rangle = 1$: this is a "folding" sequence.

By contrast, the bottom panels show a sequence that behaves very differently: the low-energy (native) conformations are not compact and large in number. Looking at the sequence itself, `PPPPPHHHHH`, it is clear that the hydrophilic "tail" does not make it possible to form a conformation with a hydrophobic core, and the number of HH topological contacts itself cannot be too high. These are all clear signs of a "non-folding" sequence.

Note that the quantities reported in [](#fig:HP_model_dill_1) are somewhat correlated: folding sequences tend tend to have few (sometimes one) low-energy compact conformations with hydrophobic cores, whereas non-folding sequences tend not to have any of these attributes.

The HP model (or its extensions) can be used to understand some general or even specific aspects of protein folding because of its simplicity. The 2-letter alphabet and the greatly reduced residue conformations due to the underlying lattice makes it possible to fully explore (at the same time) the conformational space (*i.e.* the set of all possible internal conformations of a molecule) and the sequence space (*i.e.* the set of all possible sequences of H and P residues). Both spaces have sizes that grow exponentially with $N$, and therefore the maximum length of chains that can be investigated by exhaustive enumeration is somewhat limited, but there are computational techniques that can be used to circumvent this issue.

### Folding HP sequences

If we now use a stricter definition of a "folding" sequence as a sequence that has a single native conformation, we can draw a comparison between lattice folding sequences and real proteins. First of all, how many folding sequences there are?

```{figure} figures/HP_model_dill_2.png
:name: fig:HP_model_dill_2
:align: center

The distribution of the number of sequences having $g(s)$ native conformations for chain length $N = 10, 12$ and $24$. In the figure $s = m^{\rm max}$. Adapted from [](doi:10.1021/ma00200a030).
```

[](#fig:HP_model_dill_2) shows the distributions of the number of native conformations $g(m^{\rm max})$, *i.e.* how many sequences have the given degeneracy, for different chain lengths. As $N$ increases, the distribution becomes more and more peaked on $1$: far more sequences have a single lowest-energy conformation rather than, say, 10 or 100. These results suggest that hydrophobicity alone is enough to greatly reduce the number of compact configurations into a small set of folding candidate structures ([](doi:10.1021/bi00483a001)), and confirms that the hydrophobic force is a major driving force for protein folding (see *e.g.* @finkelstein2002protein).


```{figure} figures/HP_model_comparison.png
:name: fig:HP_model_comparison
:align: center

(a) Compactness and (b) hydrophobicity as a function of the chain length for two HP lattice models and for real proteins. Adapted from [](doi:10.1002/prot.24067).
```

But to what extent lattice and real proteins are similar? I will present some results reported in [](doi:10.1002/prot.24067). Therein, an exhaustive enumeration of all conformations of chains of lengths up to 25 has been carried out with two HP models:

* In HP100, which is the one described above, only HH contacts contribute to a conformation's free energy.
* In HP211 topological contacts between H and P and P and P contribute $-1$ (*i.e.* $\epsilon_{HP} = \epsilon_{PP} = -1$) to the free energy, while HH contacts contribute $\epsilon = -2$. As a result, compared to the HP100 model, all residues feel an attraction that makes compact conformations more favourable.

Some interesting results pertaining to folding sequences are reported in [](#fig:HP_model_comparison). [](#fig:HP_model_comparison)(a) shows, as a proxy for compactness, the normalised radius of gyration $R_g/{\rm Min}(R_g)$, where ${\rm Min}(R_g)$ is the smallest radius of gyration among the structures of a given length, for the lattice and real proteins. For real proteins the authors have considered a non-redundant set of 2401 single-domain proteins from the PDB with lengths between 50 and 200 residues long. By comparing the two lattice models, it is evident that the additional attraction between non-hydrophobic residues enhances the compactness in HP211 more, but the variation in compactness, given by the error bars, is very small in both cases, and the compactness itself does not depend on the chain length. The same behaviour is observed in real proteins.

[](#fig:HP_model_comparison)(b) shows the percentage of residues that are hydrophobic in the lattice structures and for the same real proteins used for the radius of gyration. Both lattice and real proteins have an average hydrophobicity of $\approx 50\%$, independent of $m$, and the error bars and total range (dashed lines) are comparable in the two cases (although both are a bit larger for the HP model). Folding sequences in the HP model are "protein-like".

### Designable HP conformations

Now we focus on compact conformations. How many sequences have a given compact conformation as their native (*i.e.* lowest-energy) state? Does the answer depend on some properties of the selected conformation itself? In the context of the HP model, some of the first answers were provided [here](10.1126/science.273.5275.666) with a HP model with $\epsilon = -2.3$, $\epsilon_{HP} = -1$ and $\epsilon_{PP} = 0$, chosen in order to fulfill the following physical constraints, rooted in the analysis of real proteins:

1. Compact conformations have lower energies compared to non-compact ones.
2. Hydrophobic residues should have the tendency to be buried in the core ($\epsilon_{PP} > \epsilon_{HP} > \epsilon$).
3. Different types of monomers should tend to segregate ($2 \epsilon_{HP} > \epsilon_{PP} + \epsilon_{HH}$).

The output of the computation was to enumerate all the native structures of each possible sequence on a 3D 3x3x3 cube and on 2D 4x4, 5x5, 6x5 and 6x6 boxes. The result of this complete enumeration is the list of all possible sequences that can "design" a given structure, or, in other words, that have that structure as their unique native state. The size of this list, $N_s$ is a measure of the designability of a given structure.

```{figure} figures/HP_model_folding_1.png
:name: fig:HP_model_folding_1
:align: center

(a) The number of structures with the given $N_s$ for the (top) 3x3x3 cube and (bottom) 6x6 square. (b) The most designable (A) 3D and (B) 2D structures. Hydrophobic and polar residues are coloured in black and grey, respectively. (c) The probability of finding a polar residue as a function of the sequence index for the same two structures. Adapted from [](10.1126/science.273.5275.666).
```

Compact structures differ markedly in terms of their designability: there are structures that can be designed by a large number of sequences, and there are "poor" structures that can be designed by only a few or even no sequences. In fact, $\approx 10\%$ of the conformations for which no sequence has that structure as its ground state. The majority of the sequences ($\approx 95\%$ in 3D) have degenerate ground states, *i.e.* more than one compact conformation of lowest energy, and the number of structures with a given $N_s$ value decreases continuously and monotonically as $N_S$ increases. The data for the 3D and largest 2D cases is shown in [](#fig:HP_model_folding_1)(a), where the long tails of the distributions, with some structures being the ground states of thousands of sequences, highlight the presence of "highly-designable" compact conformations.

Structures with large $N_S$ exhibit specific motifs (*i.e.* secondary structures) that small $N_S$ compact structures lack. For instance, the 3D compact structures with the 10 largest $N_S$ values contain parallel running lines packed regularly and eight or nine strands (three amino acids in a row), sensibly more than the average compact structure. [](#fig:HP_model_folding_1)(b) shows the most designable structures in 3D and in 2D, where it is easy to spot the regularity of the "secondary structures" and, for the 3D structure, also the "strands".

With the simple lattice model is also possible to directly assess the effect of mutations, which in real proteins is particularly important in the context of homologous sequences (sequences related by a common ancestor). Indeed, focussing on highly designable structures and referring to the $N_s$ different sequences that fold into them as "homologous", it is possible to observe phenomena that are qualitatively similar to those observed in real proteins. For example, sequences that differ by more than half of their residues can design the same structure. Looking at the effect of mutations, the rightmost panel of [](#fig:HP_model_folding_1)(c) shows that some in very designable conformations residues are highly mutable, whereas others are highly conserved, with the conserved sites being those sites with the smallest or largest number of sides exposed to water.

```{figure} figures/HP_model_folding_2.png
:name: fig:HP_model_folding_2
:align: center

The average energy gap between the lowest-energy and first excited states of the (a) 3D 3x3 cube and (b) 2D 6x6 HP models, as well as of the (c) 2D 6x6 MJ model. Adapted from [](doi:10.1016/S1093-3263(00)00137-6)).
```

Another important feature of proteins that is also reproduced by the model is the large energy gap between the lowest-energy and first excited conformations, which stabilises the native structures from the thermodynamic point of view. In the lattice model, the average energy gap $\bar \delta_S$ is defined as the minimum energy required to change a native conformation to a different compact structure, averaged over all the $N_S$ sequences that design that native conformation. [](#fig:HP_model_folding_2) shows that there is a sudden jump (in 3D) or change of slope (in 2D) at a specific value of $N_S$. Therefore, the conformations that are highly designale also have large energy gaps, making them highly stable from a thermodynamic point of view.

:::{note}
In [](doi:10.1016/S1093-3263(00)00137-6) the energy gap was defined as the energy difference between the ground and first excited states of a sequence. However, a much better indicator for "good foldability" and overall thermodynamic stability is the so-called stability gap, which is defined as the difference between the energy of the native state and the **average energy** of compact (misfolded) states (see *e.g.* [](doi:10.1146/annurev.physchem.48.1.545)).
:::

As shown in the rightmost panel of [](#fig:HP_model_folding_2), the same behaviour is also observed when a 20-letter (rather than 2-letter) alphabet is used, with the interactions between the residues being modelled using the [MJ interaction matrix](doi:10.1021/ma00145a039), which was derived by analysing the contact frequencies of amino acids from protein crystal structures.

Although the ingredients in the simple lattice model are minimal, lacking many of the factors that contribute to the structure of proteins, the results that have been obtained with it contributed to the understanding that there is a link between designability and thermodynamic stability: "From an evolutionary point of view, highly designable structures are more likely to have been chosen through random selection of sequences in the primordial age, and they are stable against mutations" ([](10.1126/science.273.5275.666)).

[^2D_HP]: For the chain lengths usually considered, the ratio of the numbers of residues exposed to the solvent or buried in the protein core is closer to that in real proteins compared to 3D lattices (see [](doi:10.1002/prot.24067) and references therein).
[^HP_length]: Which values are "too large" depend on the hardware and algorithms used to write the code.
[^HP_notation_abuse]: Here I'm abusing the notation by using the same symbol $Q$ to represent two different (but very related) quantities.

:::{seealso} Python implementation
Head over [here](./notebooks/HP_model.ipynb) for a Jupyter notebook containing the code to exhaustively explore the conformation and sequence space of the 2D HP model.
:::

## Sequence alignment

:::{seealso} Python implementation
Head over [here](./notebooks/sequence_alignment.ipynb) for Jupyter notebook containing code implementing the various alignment algorithms discussed in this section.
:::

## AlphaFold

# Nucleic acids

## Secondary structure

## NUPACK and ViennaRNA


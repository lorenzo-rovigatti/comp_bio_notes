---
title: Structure and folding prediction
exports:
   - format: pdf
---

# Proteins

```{tip}
The main references for this part are @finkelstein2002protein and @bialek2012biophysics. Note that here I refer mostly to globular proteins.
```

Let's consider a $100$-residue peptide chain with a unique folded state. If each amino acid had only two available conformations, the number of configurations available to the chain would be $2^{100} \sim 10^{30}$. If the time required to switch between any two configurations was $10^{-12}$ s (a picosecond), and we assume that no configuration is visited twice, it would take approximately $10^{18}$ seconds to explore all the "phase space" (and therefore, on average, to fold): this is close to the age of the universe! In reality, proteins fold on time scales ranging from microseconds to hours. This is the gist of the famous ["Levinthal's paradox"](https://en.wikipedia.org/wiki/Levinthal%27s_paradox), which implies that the search for the folded structure does not happen randomly, but it is guided by the energy surface of residue-residue interactions.

Proteins can be denatured by changing the solution conditions. Very high or low temperature and pH, or high concentration of a denaturant (like urea or guanine dihydrochloride) can "unfold" a protein, which loses its solid-like, native structure, as well as its ability to perform its biological function. Investigations of the unfolding of small globular proteins showed that, as the denaturing agent (*e.g.* temperature or denaturant concentration) increases, many of the properties of the protein go through an "S-shaped" change, which is characteristic of cooperative transitions.

Furthermore, calorimetric studies show that the denaturation transition is an "all-or-none" transition, *i.e.* that the "melting unit" of the transition is the protein as a whole, and not some subparts. This of course applies to single-domain small proteins, or to the single domains of larger proteins. Here "all-or-none" means that the protein exist only in one of two states (native or denaturated), with all the other "intermediate" states being essentially unpopulated at equilibrium. Therefore, this transition is the microscopic equivalent of a first-order phase transition (*e.g.* boiling or melting) occurring in a macroscopic system. Of course, since proteins are finite systems and therefore very far from the thermodynamic limit, this is not a true phase transition, as the energy jump is continuous, and the transition width is finite.

:::{note} The van't Hoff criterion
Let $N$ and $D$ be two states in chemical equilibrium, *i.e.* $N \rightleftharpoons D$, at a certain temperature $T$. Recalling the two-state equilibrium concepts derived [before](#sec:two_state), and defining $p(T)$ as the fraction of, say, state $N$, we have

$$
K_{\rm eff}(T) = \frac{p(T)}{1 - p(T)} = e^{\beta \Delta G_{ND}},
$$

where $K_{\rm eff}(T)$ is the effective dissociation constant at temperature $T$ and $\Delta G_{ND} = \Delta H_{ND} - T \Delta S_{ND}$ is the free-energy differences between the two states. By simple algebra we have that

$$
\Delta H_{ND} = -k_B T^2 \frac{d \log K_{\rm eff}(T)}{dT}.
$$

Note that $K_{\rm eff}(T)$ is connected to $p(T)$, which can be measured experimentally (for instance with Circular Dicroism). The value of $\Delta H_{ND}$ obtained can be interpreted as the heat consumed by the "melting unit" associated to the transition.

The enthalpy change can also be estimated by integrating the heat capacity, measured with calorimetry, over the range of temperature associated to the transition. This is the heat consumed by all proteins during the transition, $\Delta H_{\rm tot} = N_p \Delta H_{\rm cal}$, where $N_p = m / M$ is the number of proteins, with $m$ being the total mass and $M$ the protein's molecular mass.

If $\Delta H_{\rm cal} = \Delta H_{ND}$, it means that the "melting unit" is the whole protein, and therefore the transition is of the "all-or-none" type, as it is the case of small globular proteins.

Note that there are [some subtleties](10.1110/ps.8.5.1064) that should be taken into account when applying this method, but the general concept is still useful, provided that the results are confirmed by independent experimental measurements.
:::

```{figure} figures/molten_globule.png
:name: fig:molten_globule
:align: center

(a) Phase diagram of the conformational states of lysozyme at pH 1.7 as a function of denaturant (guanine dihydrochloride) concentration and temperature. The red lines correspond to the mid transition, the dashed lines outline the transition zones. (b) Schematic model of the native and molten globule protein states. Here the protein consists of only two helices connected by a loop, and the side chains are shown as shaded regions. Adapted from @finkelstein2002protein.
```

But how does the denatured state look like? As discovered by using a plethora of experimental methods, which often seemed contradicting each other, the answer depends on the denaturing conditions. [](#fig:molten_globule)(a) shows the phase diagram of lysozyme (a globular small protein) for different temperatures and concentrations of a denaturing agent. Increasing the latter leads to a transition to a disordered coil state where essentially no secondary structure is present. By contrast, the effect of temperature is qualitatively different, as the protein partially melts in a state (the *molten globule*, MG) that retains most of its secondary structure, but it is not solid like and it cannot perform any biological function. A schematic of the difference between the native and molten globule states is shown in [](#fig:molten_globule)(b).

The cooperative transition that leads to the molten globule is due to the high packing of the side chains, which is energetically favourable but entropically penalised, since it impedes the small-scale backbone fluctuations that trap many of the internal degrees of freedom of the protein. The liberation of these small-scale fluctations, and the resulting entropy gain, requires a slight degree of swelling, which is why the MG is not much larger than the native state. However, such a swelling is large enough to greatly decrease the van der Waals attractions, which are strongly dependent on the distance, and to let solvent (water) molecules in the core.

In general, not all proteins behave like this, as there exist some (usually small) proteins that unfold directly into a coil, others that form the molten globule under the effect of specific denaturants, and coils under the effect of others, *etc.* However, the unfolding transition has nearly always an "all-or-none", cooperative nature.

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

(a) An example of frustrated interactions: if the particle-particle interaction is such that only like colours want to be close to each other, there is no way of arranging the three spheres so that all favourable contacts are made. As a result, the ground state is degenerate. (b) A schematic energy landscape of a 60 amino-acid helical protein. $A$ and $Q$ are the fractions of correct dihedral angles in the backbone and correct native-like contacts, respectively, $\Delta E$ is the "ruggedness" of the landscape, and $\delta E_s$ is the energy gap. Taken from [](doi:10.1073/pnas.92.8.3626).
```

All these concepts are not merely qualitative, but have been grounded in theory by using sophisticated statistical mechanics approaches that have generated the "funnel landscape" folding picture, which is schematically shown in [](#fig:frustration)(b).

The basic idea is to leverage the concept of frustration, which happens when the interactions between the different parts of a system are such that they cannot be satisfied, from the energetic point of view, all at the same time. If this is the case, then there is no single ground state, but rather many low-energy stable states that, in a many-body system, are uncorrelated from each other and separated by (possibly large) energy barriers. An example is provided in [](#fig:frustration)(a). The resulting "energy surface" (or landscape) is said to be rugged (or rough), and generates a dynamics where the system sits in a valley for a long time before being able to jump over the barrier and move to a different stable conformation. This is a hallmark of glassiness.

To draw a parallel with proteins, we can consider a random heteropolymer, where the interactions among the different amino acids will be frustrated, blocking the system from finding a single well isolated folded structure of minimum energy. A candidate principle for selecting functional sequences is thus the minimization of this frustration. Of course even in the absence of frustration, there are still energetic barriers on the path towards the native conformation due to local structural rearrangements which still give raise to a rugged landscape, slowing down the kinetics towards the folded state. This scenario has come to be called a folding "funnel", emphasizing that there is a single dominant valley in the energy landscape, into which all initial configurations of the system will be drawn.

A classic funnel diagram, such as the one shown in [](#fig:frustration)(b), represents an energy-entropy landscape, with the width being a measure of the entropy, whereas the depth represents both the energy and two correlated order parameters $Q$ and $A$, which are the fractions of native contacts and correct dihedral angles in the protein backbone. Although in reality the landscape is multidimensional, here the projection attempts to retain its main features, such as the barrier heights, which are a measure of the "ruggedness" or "roughness" of the landscape. The typical height of these barriers, together with the energy gap and the number of available conformations, are the three main parameters of the Random Energy Model, which is one of the main theoretical tools in this context.

(sec:zimm-bragg)=
## The Zimm-Bragg model for the helix-to-coil transition

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

In general the operations carried out to compute the $i = 3$ partial partition functions can be applied to any other value of $i$, meaning that we can write down recursive relationships connecting $H_{i}$ and $C_{i}$ to $H_{i-1}$ and $C_{i-1}$. These relationships can be neatly expressed in matricial form:

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

Some interesting results pertaining to folding sequences are reported in [](#fig:HP_model_comparison). [](#fig:HP_model_comparison)(a) shows, as a proxy for compactness, the normalised radius of gyration $R_g/\min(R_g)$, where $\min(R_g)$ is the smallest radius of gyration among the structures of a given length, for the lattice and real proteins. For real proteins the authors have considered a non-redundant set of 2401 single-domain proteins from the PDB with lengths between 50 and 200 residues long. By comparing the two lattice models, it is evident that the additional attraction between non-hydrophobic residues enhances the compactness in HP211 more, but the variation in compactness, given by the error bars, is very small in both cases, and the compactness itself does not depend on the chain length. The same behaviour is observed in real proteins.

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

```{tip}
The main references for this part are @kellis_comp_bio and @miklos2016introduction.
```

Simple models are useful to understand the underlying physics of some particular phenomena. However, how can we understand something very specific, like what is the 3D structure of a particular sequence? The simplest way is to look for similarities: if we already have a list of sequence $\to$ structure connections, we can try to look whether the new sequence, for which the 3D structure is unknown, is similar, and to what degree, to another for which the 3D structure is already known. This operation is called "sequence alignment" (SA).

Sequence alignment is a fundamental technique in bioinformatics. The primary goal of sequence alignment is to identify regions of similarity that may indicate functional, structural, or evolutionary relationships between the sequences being compared. In the context of proteins, sequence alignment can be useful for reasons that go beyond the 3D structure prediction:

* It helps predict the function of unknown proteins based on their similarity to known proteins.
* It aids in understanding evolutionary relationships by identifying conserved sequences among different species, allowing the construction of phylogenetic trees.
* It can identify conserved domains that are crucial for the function of proteins.
* It helps in identifying potential drug targets by finding unique sequences in pathogens that differ from the host.
* It assists in annotating genomes by finding homologous sequences, thus providing insights into gene function and regulation.

Turning a biological problem into an algorithm that can be solved on a computer requires making some assumptions in order to obtain a model with relative simplicity and tractability. In practice, for our sequence alignment model we ignore realistic events that occur with low probability (*e.g.* duplications) and focus on the following three mechanisms through which sequences vary:

1. A (point) mutation, or substitution, occurs when some amino acid in a sequence changes to some other amino acid during the course of evolution.
2. A deletion occurs when an amino acid is deleted from a sequence during the course of evolution.
3. A insertion occurs when an amino acid is added to a sequence during the course of evolution.

There are many possible sequences of events that could change one genome into another, and we wish to establish an optimality criterion that allows us to pick the "best" series of events describing changes between sequences. We choose to invoke Occam's razor and select a maximum parsimony method as our optimality criterion[^parsimony]: we wish to minimize the number of events used to explain the differences between two sequences. In practice, it is found that point mutations are more likely to occur than insertions and deletions, and certain mutations are more likely than others. We now develop an algorithmic framework where these concepts can be incorporated explicitly by introducing parameters that take into account the "biological cost" of each of these changes.

(sec:homology)=
### Homology 

While here my focus is on protein structure, I also want to point out that in bioinformatics one of the key goals of sequence alignment is to identify homologous sequences (*e.g.*, genes) in a genome. Two sequences that are homologous are evolutionarily related, specifically by descent from a common ancestor. The two primary types of homologs are orthologous and paralogous. While they are outside the scope of these notes, I note here that other forms of homology exist (*e.g.*, xenologs, which are genes in different species that have arisen through horizontal gene transfer[^HGT] rather than through vertical inheritance from a common ancestor, and often perform similar functions in their respective organisms despite their distinct evolutionary origins). 

* Orthologs arise from speciation events, leading to two organisms with a copy of the same gene. For example, when a single species A speciates into two species B and C, there are genes in species B and C that descend from a common gene in species A, and these genes in B and C are orthologous (the genes continue to evolve independent of each other, but still perform the same relative function).
* Paralogs arise from duplication events within a species. For example, when a gene duplication occurs in some species A, the species has an original gene B and a gene copy B', and the genes B and B' are paralogus.

Generally, orthologous sequences between two species will be more closely related to each other than paralogous sequences. This occurs because orthologous typically (although not always) preserve function over time, whereas paralogous often change over time, for example by specializing a gene's (sub)function or by evolving a new function. As a result, determining orthologous sequences is generally more important than identifying paralogous sequences when gauging evolutionary relatedness.

[^HGT]: Horizontal gene transfer is the movement of genetic material between organisms in a manner other than traditional reproduction (vertical inheritance). This process allows genes to be transferred across different species, leading to genetic diversity and the rapid acquisition of new traits, and can occur through several mechanisms.

### Definitions

Following @kellis_comp_bio, we introduce a simplified version of the alignment problem and make it more complex by adding the required ingredients. First, some definitions:

Sequence
: In mathematics (and combinatorics) a sequence is a series of characters taken from an alphabet $\Sigma$. For what concerns us, DNA molecules are sequences over the alphabet $\lbrace A, C, G, T\rbrace$, RNA sequences over the alphabet $\lbrace A, C, G, U\rbrace$, and proteins are sequences over the alphabet $\lbrace A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V \rbrace$.

Substring
: A consecutive part of a sequence. A sequence of length $N$ has $N$ substrings of length 1, $N(N-1)/2$ substrings of length 2, $\ldots$, $2$ substrings of length $(N - 1)$. Given the sequence "SUPERCALIFRAGILISTICHESPIRALIDOSO", "FRAGILI" is a substring, but "SUPERDOSO" is not.

Subsequence
: A set of characters of a sequence. The characters do not have to be consecutive, but they should be ordered like in the original sequence. Given the sequence "SUPERCALIFRAGILISTICHESPIRALIDOSO", "FRAGILI" and "SUPERDOSO" are two subsequences, but "SOSUPER" is not. Formally, given a sequence $S = (s_1, s_2, \ldots, s_N)$, $K = (z_1, z_2, \ldots, s_n)$, with $n \leq N$, is a subsequence of $S$ if there exists a strictly increasing sequence $i_1 < i_2 < \ldots < i_n$ of indices of $S$ such that $S_{i_j} = K_j$ for each $j: 1 \leq j \leq n$.

Sequence alignment
: An alignment of two sequences $S$ and $T$, defined over the same alphabet $\Sigma$, is a $2 \times L$ table containing either a character from $\Sigma$, or the gap symbol "`-`". In the table there is no column in which both characters are "`-`", and the non-gap characters in the first line give back sequence $S$, while the non-gap characters in the second line give back sequence $T$.

### Problem formulation

We warm up by considering the following simple problem: given two sequences, $S$ and $T$, what is their longest common substring? We take the two following DNA fragments as examples: `ACGTCATCA` and `TAGTGTCA`. The problem can be solved by aligning the first characters of the two sequences, and then shifting one, say $T$, by an integer amount $i \in \mathcal{N}$. By computing the number of matching characters for every $i$ we can find the optimal alignment, which in this case is given by $i = -2$:

```
--ACGTCATCA
--xx||||---
TAGTGTCA---
```

The first and third rows of the diagram above describe the alignment itself, as introduced earlier, showing that gaps are characters that align to a space, while the central row shows the character matches: here `x` and `|` stand for mismatches and matches, respectively. In this example, the longest common substring is `GTCA`. The algorithm I just described has a complexity of $\mathcal{O}(n^2)$, where $n$ is the length of the shortest sequence.

A more complicated problem is to find the longest common *subsequence* (LCS) between two sequences. In the language we are introducing, this means allowing internal gaps in the alignment. The LCS between the $S$ and $T$ sequences defined earlier is

```
-ACGTCATCA
-|-||x-|||
TA-GTG-TCA
```

In this case the LCS is `AGTTCA`, which is longer than the longest commond substring.

The LCS problem as formulated here is a particular case of the full sequence-alignment problem. In order to generalise the problem, we introduce a cost function that makes it possible to recast it in terms of an optimisation problem. Given two sequences $S$ and $T$, we define $\delta(S, T)$ as the "biological cost" of turning $S$ into $T$. For the case just considered, I implicitly used as a cost function the [Levenshtein (or edit) distance](https://en.wikipedia.org/wiki/Levenshtein_distance), which is defined as the minimum number of single-character edits required to change one string into another. In other words, the problem has been formulated by implicitly considering that all possible operations (mutations, insertions and deletions) are equally likely. In a more generalised formulation, each edit operation is associated to a specific cost (penalty) that should reflect its biological occurrence, as discussed [later](#sec:substitution_matrices).

What about gaps? In general, the cost of creating a gap depends on many variables. There are varying degrees of approximations that can be taken in order to simplify the problem. At the zero-th order each gap can be assumed to have a fixed cost, as we implicitly did above. This is called the "linear gap penalty". An improvement can be done by considering that, biologically, the cost of creating a gap is more expensive than the cost of extending an already created gap. This is taken into account by the "affine gap penalty" method, whereby there is a large initial cost for opening a gap, and then a small incremental cost for each gap extension.
There are other (more complicated and more context-specific) methods that can be considered and that we will not analyse further here, such as the general gap penalty, which allows for any cost function, and the frame-aware gap penalty, which is applied to DNA sequences and tailors the cost function to take into account disruptions to the coding frame[^frame-aware].

### The Needleman-Wunsch algorithm

Once we allow for gaps, the enumeration of all possible alignments (as done for the longest common substring method) becomes unfeasible. Indeed, the number of non-boring alignments, *i.e.* alignments where gaps are always paired with characters, in a sequence of size $N$ containing $M$ gaps (where $N \approx M$) can be estimated as

$$
\binom{N + M}{M} = \frac{(N + M)!}{N!M!} \approx \frac{(2N)!}{(N!)^2} \approx \frac{2^{2N}}{\sqrt{\pi n}},
$$

where we used the second-order Stirling's approximation, $\log(N!) \approx N \log(N) - N + \frac{1}{2}\log (2 \pi N)$. Note that this number grows **very** fast: for $N = 30$ it is already larger than $10^{17}$. Considering that we also have to compute the score of each alignment, it is clear that the problem is untractacle with a brute-force method. Here is where dynamic programming enters the field.

Suppose we have an optimal alignment for two sequences $S$, of length $N$, and $T$, of length $M$, in which $S_i$ matches $T_j$. The alignment can be conceptually split into three parts:

1. Left subalignment: The alignment of the subsequences $(S_1 , \ldots , S_{i−1})$​ and $(T_1 , \ldots , T_{j −1} )$​.
2. Middle match/mismatch: The alignment of $S_i$​ with $T_j$​ (this could be a match or a mismatch).
3. Right subalignment: The alignment of the subsequences $(S_{i+1}, \ldots, S_N)$​ and $(T_{j+1} , \ldots, T_M)$​.

The overall alignment score is the sum of the scores from these three parts: the score of the left subalignment, the score for aligning $S_i$​ with $T_j$​, and the score of the right subalignment. Therefore, if the overall alignment is optimal, both the left and right subalignments must themselves be optimal. This follows from a cut-and-paste argument: imagine that one of the two subalignment was not optimal. This would mean there exists another alignment of these subsequences with a higher score. If such a better alignment for the subsequences existed, we could "cut" the suboptimal subalignment from the original alignment and "paste" in this better alignment. This would result in a new alignment for $S$ and $T$ with a higher total score than the supposed "optimal" alignment. This is a contradiction because the original alignment was assumed to be optimal. Since the same argument applies to both subalignments, the overall alignment's optimality depends on the optimality of both its left and right subalignments. Of course, this is true only if the scores are additive, so that the score of the overall alignment is the sum of the scores of the alignments of the subsequences. The implicit assumption is that the sub-problems of computing the optimal scoring alignments of the subsequences are independent.

Let $F_{i,j}$ be the score of the optimal alignment of $(S_1 , \ldots , S_i)$ and $(T_1 , \ldots , T_j)$. Since $i \in [0, N]$ and $j \in [0, M]$, the matrix $\hat F$ storing the solutions (*i.e.* optimal scores) of the subproblems has a size of $(M + 1) \times (N + 1)$.

We can compute the optimal solution for a subproblem by making a locally optimal choice based on the results from the smaller subproblems. Thus, we need to establish a recursive function that shows how the solution to a given problem depends on its subproblems and can be used to fill up the matrix $\hat F$. At each iteration we consider the four possibilities (insert, delete, substitute, match), and evaluate each of them based on the results we have computed for smaller subproblems.

We start by considering the linear gap penalty model, and define $d$, with $d < 0$, as the cost of a gap. We are now equipped to set the values of the elements of the matrix. Let's consider the first row: the value of $F_{0,j}$ is the cost of aligning a sequence of length $0$ (taken from $S$) to a sequence of length $j$ (taken from $T$), which can be obtained only by adding $j$ gaps, yielding $F_{0,j} = jd$. Likewise, for the first column we have $F_{i,0} = id$. Then, we traverse the matrix element by element. Let's consider a generic element $F_{ij}$: this is the cost of aligning the first $i$ characters of $S$ to the first $j$ characters of $T$. There are three ways we can obtain this alignment:

* a gap is added to $S$: the total cost is $F_{i-1, j} + d$;
* a gap is added to $T$: the total cost is $F_{i, j - 1} + d$;
* $S_i$ and $T_j$ are matched: the total cost is $F_{i - 1, j - 1} + s(S_i, T_j)$, where $s(x, y)$ is the cost of (mis)matching $x$ and $y$.

Leveraging the cut-and-paste argument sketched above, the optimal cost is given by the largest of the three values. Formally,

$$
F_{i,j} = \max
\begin{cases}
F_{i - 1, j} + d \\
F_{i, j - 1} + d \\
F_{i - 1, j - 1} + s(S_i, T_j).
\end{cases}
$$ (eq:needleman_wunsch)

This is a recursive relation: computing the value of any $F_{i,j}$ requires the knowledge of the values of its left, top, and top-left neighbours. Therefore, the fill-in phase amounts to traversing the table in row or column major order, or even diagonally from the top left cell to the bottom right cell. After traversing the matrix, the optimal score for the alignment is given by the bottom-right element, $F_{MN}$. In order to obtain the actual alignment we have to traceback through the choices made during the fill-in phase. It is helpful to maintain a pointer for each cell while filling up the table that shows which choice was made to get the score for that cell. Then the pointers can be followed backwards to reconstruct the optimal alignment. 

The complexity analysis of this algorithm is straightforward. Each update takes $\mathcal{O}(1)$ time, and since there are $MN$
elements in the matrix $\hat F$, the total running time is $\mathcal{O}(MN)$. Similarly, the total storage space is $\mathcal{O}(MN)$. For the more general case where the update rule is more complicated, the running time may be more expensive. For instance, if the update rule requires testing all sizes of gaps (*e.g.* the cost of a gap is not linear), then the running time would be $\mathcal{O}(MN(M + N)$).

:::{tip} A simple example
Consider the two DNA sequences $S = AGT$ and $T = AAGC$, a gap penalty $d = -2$, and $s(x, y) = \pm 1$, where the plus and minus signs are used for matches and mismatches, respectively.

```{figure} figures/needleman_wunsch_example.png
:name: fig:needleman_wunsch_example
:align: center

The dynamic programming table of the example during three stages of the fill-in phase: (a) at the beginning, (b) halfway through, and (c) at the end. The arrows point from each box to the box that has been used to compute its score. Panel (b) highlights one $(i, j)$ pair for which two of the cases of Eq. [](#eq:needleman_wunsch) give the same value. Each "branching" such as this one doubles the number of optimal alignments. Tracing back the arrows from the bottom-right corner (in green) to the top-left corner in grey) the optimal alignment(s) can be reconstructed.
```

[](#fig:needleman_wunsch_example) shows how the dynamic programming table of the problem looks when initialised, halfway through the fill-in phase, and after being fully traversed. Once the full table has been computed, the bottom-right box contains the optimal score, while the solutions (*i.e.* the optimal alignments) can be reconstructed by tracing back the matrix, following the arrows. Note that sometimes the maximum value computed from Eq. [](#eq:needleman_wunsch) is degenerate, *i.e.* the same value can be obtained by performing several operations (see [](#fig:needleman_wunsch_example)(b)). In this case there are multiple optimal alignments. For instance, for the simple example shown here there are two optimal alignments:

```
A-GT
|-|x
AAGC

-AGT
-||x
AAGC
```
:::

### Local alignment: the Smith-Waterman algorithm

The Needleman-Wunsch algorithm find the best possible alignment across the entire length of two sequences. It tries to align every character from the start to the end of the sequences, which means both sequences are considered in their entirety. This is called "global alignment", and it is most useful when the sequences being compared are of similar length and are expected to be homologous across their entire length.

Local alignment, on the other hand, focuses on finding the best alignment within a subset of the sequences. It identifies regions of similarity between the two sequences and aligns only those regions, ignoring the parts of the sequences that do not match well. Local alignment is particularly useful when comparing sequences that may only share a segment of similarity, such as when comparing domains within proteins, detecting conserved motifs, or identifying homologous regions in sequences that may not be overall similar[^local_alignment_DNA], which is why is very useful for the prediction of the 3D structure of proteins (or protein subdomains).

The most used method for local alignment is the Smith-Waterman algorithm, which is a modification of the Needleman-Wunsch algorithm. The key difference between the two lies in how the scoring matrices are constructed and scored. In Needleman-Wunsch, every cell in the matrix is filled to reflect the best global alignment, with the final alignment score found in the bottom-right corner of the matrix. Smith-Waterman, on the other hand, sets any negative scores to zero, which allows the algorithm to "reset" when the alignment quality dips. The highest score in the matrix indicates the end of the best local alignment, which is then traced back to identify the optimal aligned subsequence. This approach ensures that only the most relevant, highest-scoring local alignments are highlighted. Here is how the Smith-Waterman algorithm looks like in practice:

1. Initialisation: since a local alignment can start anywhere, the first row and column in the matrix are set to zeros, *i.e.* $F_{0,j} = jd$, $F_{i,0} = id$.
2. Iteration $\forall (i, j)$: this step is modified so that the score is never allowed to become negative but it is reset to zero. This is done by slightly modifying Eq. [](#eq:needleman_wunsch) as follows:
$$
F_{i,j} = \max
\begin{cases}
0\\
F_{i - 1, j} + d \\
F_{i, j - 1} + d \\
F_{i - 1, j - 1} + s(S_i, T_j).
\end{cases}
$$
3. Trace-back: starts from the position of the maximal number in the table and proceeds until a zero is encountered.

:::{tip} The Needleman-Wunsch algorithm
For reference, this is summary of the Needleman-Wunsch algorithm:

1. Initialisation: $F_{0,j} = jd$, $F_{i,0} = id$.
2. Iteration $\forall (i, j)$: Eq. [](#eq:needleman_wunsch).
3. Trace-back: starts from the bottom-right value and stops at the top-left corner.
:::

[^parsimony]: Note that this is not the only possible choice: we could choose a probabilistic method, for example using Hidden Markov Models (HMMs), that would assign a probability measure over the space of possible event paths and use other methods for evaluating alignments (*e.g.*, Bayesian methods).
[^frame-aware]: Indels (shorthand for "insertions/deletions") that cause frame-shifts in functional elements generally cause important phenotypic modifications.
[^local_alignment_DNA]: It is also valuable in cases where one sequence may be a subsequence of another, like when searching for a gene within a cromosome or a whole genome.

(sec:affine_gaps)=
### Affine gap penalty

For both the global and local alignment algorithms introduced we have used a linear gap penalty: the cost of adding a gap is constant, regardless of the nature of the aligned character, or of the length of the gap. From the biological point of view, this means that an indel of *any* length is considered as the result of independent insertions or deletions. However, in reality long indels can form in single evolutionary steps, and in these cases the linear gap model overestimates their cost. To overcome this issue, more complex gap penalty functions have been introduced. As mentioned before, a generic gap penalty function would result in an algorithmic complexity worse than $\mathcal{O}(N^2)$ (where for simplicity I'm considering two sequences of the same length). Let's see why. Any gap penalty function can be implemented in the Needleman-Wunsch or Smith-Waterman algorithms by changing the recursive rule, which can be done trivially for element $F_{i,j}$ by evaluating terms such as $\max_{0 \leq k \leq i} \lbrace F_{k,j} + d_{i - k} \rbrace$ and $\max_{0 \leq k \leq i} \lbrace F_{i,k} + d_{i - k} \rbrace$. However, this means that the update of every cell of the dynamic programming matrix would take $\mathcal{O}(N)$ instead of $\mathcal{O}(1)$, bringing the algorithmic complexity up to $\mathcal{O}(N^3)$ (for $N = M$ sequences).

However, the computational cost can be mitigated by using particular gap penalty functions. Here I will present the most common variant, which is known as the *affine* gap penalty. In this model, the cost of a gap of size $k$ is

$$
d_k = o + (k - 1) e,
$$

where $o$ and $e$ are the opening and extension penalties, respectively, and $o > e$. To incorporate affine-gap penalties into the Needleman-Wunsch or Smith-Waterman algorithms in an efficient manner, the primary adjustment involves tracking whether consecutive gaps occur in the alignment. This requires the alignment process to be split into three distinct cases: insertions, deletions, and matches/mismatches. Instead of using a single dynamic programming table as in the linear-gap penalty approach, three separate tables are used: $\hat I$ for insertions, $\hat D$ for deletions, and $\hat F$ for the overall score.

In this setup, each entry in the $\hat I$ table, denoted by $i_{i,j}$, stores the best alignment score when the last column includes an insertion (*i.e.* a gap in the first sequence), and each entry in the $\hat D $ table, $ d_{i,j} $, captures the best score for alignments ending in a deletion (*i.e.* a gap in the second sequence). As before, $F_{i,j} $ records the overall score. The recursive update rules of the three tables are

\begin{align}
I_{i,j} &= \max \begin{cases}
F_{i-1,j} + o\\
I_{i-1,j} + e
\end{cases}\\

D_{i,j} &= \max \begin{cases}
F_{i,j-1} + o\\
D_{i,j-1} + e
\end{cases}\\

F_{i,j} &= \max \begin{cases}
F_{i-1,j-1} + s(S_i, T_j)\\
I_{i, j}\\
D_{i,j}.
\end{cases}
\end{align}

Note that with this algorithm each cell update is $\mathcal{O}(1)$, and therefore the overall complexity remains the same ($\mathcal{O}(NM)$ or $\mathcal{O}(N^2)$ for same-length sequences).

:::{seealso} Python implementation
Head over [here](./notebooks/sequence_alignment.ipynb) for Jupyter notebook containing code implementing the Needleman-Wunsch and Smith-Waterman algorithms, with linear and affine gap penalty functions.
:::

(sec:substitution_matrices)=
### Substitution matrices

How is the substitution matrix $s(x, y)$ determined? One possibility is to leverage what we know about the biological processes that underlie mutations. For instance, in DNA the biological reasoning behind the scoring decision can be linked to the probabilities of bases being transcribed incorrectly during polymerization. We already know that of the four nucleotide bases, A and G are purines, while C and T are pyrimidines. Thus, DNA polymerase is much more likely to confuse two purines or two pyrimidines since they are similar in structure. As a result, a simple improvement over the uniform cost function used above is the following scoring matrix for matches and mismatches:

||A|G|T|C|
|-|-|-|-|-|
|A|+1|-0.5|-1|-1|
|G|-0.5|+1|-1|-1|
|T|-1|-1|+1|-0.5|
|C|-1|-1|-0.5|+1|

Here a mismatch between like-nucleotides (*e.g.* A and G) is less expensive than one between unlike nucleotides (*e.g.* A and C).

However, we can be more quantitative by taking a probabilistic approach. Here I will describe how the widely-used [BLOSUM matrices](https://en.wikipedia.org/wiki/BLOSUM) are built (see also the [original paper](doi:10.1073/pnas.89.22.10915)). We assume that alignment score reflects the probability that two similar sequences are [homologous](#sec:homology). Thus, given two sequences $S$ and $T$, we define two hypotheses:

1. The alignment between the two sequences is due to chance and the sequences are, in fact, unrelated.
2. The alignment is due to common ancestry and the sequences are actually related.

In order to compare two hypotheses, a good score is given by the logarithm of the ratio of their likelihoods (the so-called log-odds score). Therefore, for our case the alignment score is 

$$
A \equiv \log \frac{P(S, T|R)}{P(S, T|U)}.
$$

Under the rather strict assumption ("biologically dubious, but mathematically convenient", as aptly put in [](doi:10.1038/nbt0804-1035)) that pairs of aligned residues are statistically independent of each other, and thanks to the properties of logarithms and probabilities, the overall alignment score is the sum of the log-scores of the individual residue pairs. Considering 20 amino acids, there are 400 such log-scores, which form a 20x20 substitution (score) matrix. Each entry of the matrix takes the form

$$
s(x, y) = \frac{1}{\lambda} \log \frac{p_{xy}}{q_x q_y},
$$

where 

* $p_{xy}$ is the likelihood of the second hypothesis (the two residues are correlated);
* $q_x q_y$ is the likelihood of the first (null) hypothesis: the two residues are there by chance, their occurrence is unrelated and therefore the likelihood factorises in two terms that account for the average probability of observing those two residues in any protein;
* $\lambda$ is an overall scaling factor that is used to obtain values that can be rounded off to integers.

A positive score indicates that 

If $p_{xy} > q_x q_y$ (and therefore $s(x, y)$ is positive), the substitution $x \to y$, and therefore their alignment in homologous sequences, occurs more frequently than would be expected by chance, suggesting that this substitution is favored in evolution. These substitutions are called "conservative". As noted in [](doi:10.1038/nbt0804-1035), this definition is purely statistical and has nothing directly to do with amino acid structure or biochemistry. Likewise, substitutions with $p_{xy} < q_x q_y$, and therefore negative values of $s(x, y)$, are termed "nonconservative".

The procedure applied to compute the $p_{xy}$ and $q_x$ values starts with a collection of protein sequences.

1. The sequences are grouped into families based on their evolutionary relatedness. A common source for these families is the BLOCKS database, which contains aligned, ungapped regions of proteins that are highly conserved across members of the family. 
2. Within these protein families, highly conserved regions (blocks) are identified. These regions are short segments of amino acid sequences that are believed to be important for the protein's function and are conserved across different species. 
3. Once blocks are identified, the sequences within each block are clustered, choosing a threshold value $V$: two sequences that share more than this percentage of identical amino acids identity are clustered together, and only one representative from each cluster is used to avoid over-representation of very similar sequences.
4. Within each block, the actual counts of amino acid pairs are made.
    * The background frequencies $q_x$ are estimated by computing the overall frequency of each amino acid $x$ in the sequences, not considering any substitutions.
    * To evaluate the $p_{xy}$ terms, for each position in the aligned block, we count how many times one amino acid (say, Alanine) is substituted by another amino acid (say, Glycine) in the aligned positions across the different sequences. For instance, if there is an alignment of 5 sequences, for each single position we make $\binom{5}{2} = 10$ comparisons.

Of course, the final scores depend on $V$: if $V$ is large, protein blocks that are still rather similar will be considered belonging to different clusters, and therefore compared to each other to derive the scores; the resulting matrix will be more sensitive in identifying closely related sequences but less effective for more distantly related sequences (and vice versa for small values of $V$). Common values for $V$ are $45\%$, $62\%$, and $80\%$, which yield the matrices BLOSUM45, BLOSUM62, and BLOSUM80, with BLOSUM62 being the de-facto standard (and default) one. 

```{figure} figures/BLOSUM62.svg
:name: fig:BLOSUM62
:align: center

The BLOSUM62 matrix. The amino acids are grouped and coloured based on [Margaret Dayhoff's classification](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff#Table_of_Dayhoff's_encoding_of_amino_acids). Non-negative values are highlighted. Credits to [Ppgardne via Wikimedia Commons](https://upload.wikimedia.org/wikipedia/commons/f/f5/Blosum62-dayhoff-ordering.svg).
```

The BLOSUM62 matrix is shown in [](#fig:BLOSUM62). First of all, note that substitution matrices are symmetric, because the biological process of amino acid substitution is symmetric: as empirically observed, there is no preferred direction when one amino acid replaces another. Secondly, it is evident that conservative substitutions tend to be those between amino acids that are similar, as made evident by grouping the amino acids according to the classification introduced by [Margaret Dayhoff](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff#Table_of_Dayhoff's_encoding_of_amino_acids).

:::{tip} Margaret Dayhoff's classification
According to [Margaret Dayhoff](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff), one of the pioneers of biophysics and bioinformatics, amino acids can be classified in the following 6 groups, listed in the order with which they are shown in [](#fig:BLOSUM62):

|Amino acids|One-letter code|Property|Dayhoff|
|-|-|-|-|
Cysteine| C | Sulfur polymerization | a |
Glycine, Serine, Threonine, Alanine, Proline | G, S, T, A, P | Small | b |
Aspartic acid, Glutamic acid, Asparagine, Glutamine | D, E, N, Q | Acid and amide | c |
Arginine, Histidine, Lysine | R, H, K | Basic | d |
Leucine, Valine, Methionine, Isoleucine | L, V, M, I | Hydrophobic |e |
Tyrosine, Phenylalanine, Tryptophan | Y, F, W | Aromatic | f |
:::

Why the values of the matrix diagonal, representing the scores of "substituting" one amino acid with itself, are all different? 

> The rarer the amino acid is, the more surprising it would be to see two of them align together by chance. In the homologous alignment data that BLOSUM62 was trained on, leucine/leucine (L/L) pairs were in fact more common than tryptophan/tryptophan (W/W) pairs ($p_{LL} = 0.0371$, $p_{WW} = 0.0065$), but tryptophan is a much rarer amino acid ($f_L = 0.099$, $f_W = 0.013$). Run those numbers (with BLOSUM62's original $\lambda = 0.347$) and you get +3.8 for L/L and +10.5 for W/W, which were rounded to +4 and +11.
>
> -- [](doi:10.1038/nbt0804-1035)

### BLAST

The sheer volume of sequence data generated by modern-day high-throughput sequencing technologies presents a significant challenge. Databases now contain millions of nucleotide and protein sequences, each potentially spanning thousands of characters. When comparing a new sequence against these massive databases, traditional pairwise alignment methods like the ones we just discussed become computationally demanding. Indeed, performing a global or even local alignment between a query sequence and every sequence in a large database can require huge computational resources and time, especially as the number of sequences and their lengths continue to grow exponentially.

The need for a more efficient method resulted in the [most cited paper of the 1990s](doi:10.1016/S0022-2836(05)80360-2), where BLAST (Basic Local Alignment Search Tool), now a critical tool in bioinformatics, was introduced. BLAST operates by finding regions of local similarity between sequences, which is more computationally feasible and faster than aligning entire sequences globally. The algorithm requires a query sequence and a target database. First of all, the query sequence is broken down into smaller fragments (called words or $W$-mers). For each $W$-mer, a list of similar words is generated, and only those with a similarity measure that is higher than a threshold $T$ are retained add added to the final $K$-mer list. The similarity is evaluated by using a [substitution matrix](#sec:substitution_matrices) (BLOSUM62 is a common choice). Then, the target database is searched for matches to these words, extending the matches in both directions to find the best local alignments. BLAST assigns scores to these alignments based on the degree of similarity *via* a Smith-Waterman algorithm, with higher scores indicating closer matches.

:::{warning}
The Needleman-Wunsch and Smith-Waterman algorithms always find the optimal solution (*i.e.* the global minimum). By contrast, BLAST uses a heuristic approach, trading off some sensitivity for speed.
:::

:::{warning} TODO
ADD FIGURE
:::

The BLAST algorithm can be broken down into the following steps[^BLAST_wiki]

1. **Optional**: *remove low-complexity region or sequence repeats in the query sequence.* Here "Low-complexity region" means a region of a sequence composed of few kinds of elements. These regions might give high scores that confuse the program to find the actual significant sequences in the database, so they should be filtered out. The regions will be marked with an X (protein sequences) or N (nucleic acid sequences) and then be ignored by the BLAST program. To filter out the low-complexity regions, the SEG program is used for protein sequences and the program DUST is used for DNA sequences. On the other hand, the program XNU is used to mask off the tandem repeats in protein sequences.
2. *Make a $W$-letter word list of the query sequence.* The words of length $W$ ($W$-mers) in the query sequence are listed "sequentially", until the last letter of the query sequence is included. The method is illustrated in TODO. $W$ is usually 3 and 11 for a protein and a DNA sequence, respectively.
3. *List the possible matching words.* A [substitution matrix](#sec:substitution_matrices) (*e.g.* BLOSUM62) is used to match the words listed in step 2 with all the $20^W$ $W$-mers. For example, the score obtained by comparing PQG with PEG and PQA is respectively 15 and 12 with the BLOSUM62 matrix[^DNA_BLAST_words]. After that, a neighborhood word score threshold $T$ is used to reduce the number of possible matching words. The words whose scores are greater than the threshold *T* will remain in the possible matching words list, while those with lower scores will be discarded. For example, if $T = 13$ PEG is kept, but PQA is abandoned.
4.  *Organize the remaining high-scoring words into an efficient search tree.* This allows the program to rapidly compare the high-scoring words to the database sequences.
5.  *Repeat step 3 to 4 for each $W$-mer in the query sequence.*
6.  *Scan the database sequences for exact matches with the remaining high-scoring words.* The BLAST program scans the database sequences for the remaining high-scoring words, such as PEG. If an exact match is found, this match is used to seed a possible ungapped alignment between the query and database sequences.
7.  *Extend the exact matches to high-scoring segment pair (HSP).* The original version of BLAST stretches a longer alignment between the query and the database sequence in the left and right directions, from the position where the exact match occurred. The extension does not stop until the accumulated total score of the HSP begins to decrease. To save more time, a newer version of BLAST, called BLAST2 or gapped BLAST, has been developed. BLAST2 adopts a lower neighborhood word score threshold to maintain the same level of sensitivity for detecting sequence similarity. Therefore, the list of possible matching words list in step 3 becomes longer. Next, exact matched regions that are within distance $A$ from each other on the same diagonal are joined as a longer new region. Finally, the new regions are then extended by the same method as in the original version of BLAST, and the scores for each HSP of the extended regions are created by using a substitution matrix as before.
8.  *List all of the HSPs in the database whose score is high enough to be considered.* All the HSPs whose scores are greater than the empirically determined cutoff score $S$ are listed. By examining the distribution of the alignment scores modeled by comparing random sequences, a cutoff score $S$ can be determined such that its value is large enough to guarantee the significance of the remaining HSPs.
9.  *Evaluate the significance of the HSP score.* [It has been shown](doi:10.1073/pnas.87.6.2264) that the distribution of Smith-Waterman local alignment scores between two random sequences is
$$
p\left( S\ge x \right) =1-\exp \left( -KMN e^{-\lambda x } \right),
$$
where $M$ and $N$ are the length of the query and database sequences[^BLAST_MN], and the statistical parameters $\lambda$ and $K$ depend upon the substitution matrix, gap penalties, and sequence composition (the letter frequencies) and are estimated by fitting the distribution of the ungapped local alignment scores of the query sequence and of a lot (globally or locally) shuffled versions of a database sequence. Note that the validity of this distribution, known as the Gumbel extreme value distribution (EVD), has not been proven for local alignments containing gaps yet, but there is strong evidence that it works also for those cases. The expect score $E$ of a database match is the number of times that an unrelated database sequence would obtain a score $S$ higher than $x$ by chance. The expectation $E$ obtained in a search for a database of total length $N$
$$
E = K M N e^{-\lambda S}
$$
This expectation or expect value $E$ (often called $E$-score, $E$-value or $e$-value) assessing the significance of the HSP score for ungapped local alignment is reported in the BLAST results. The relation above is different if individual HSPs are combined, such as when producing gapped alignments (described below), due to the variation of the statistical parameters.
10. *Make two or more HSP regions into a longer alignment.* Sometimes, two or more HSP regions in one database sequence can be made into a longer alignment. This provides additional evidence of the relation between the query and database sequence. There are two methods, the Poisson method and the sum-of-scores method, to compare the significance of the newly combined HSP regions. Suppose that there are two combined HSP regions with the pairs of scores $(65, 40)$ and $(52, 45)$, respectively. The Poisson method gives more significance to the set with the maximal lower score $(45 > 40)$. However, the sum-of-scores method prefers the first set, because $65+40=105$ is greater than $52+45 = 97$. The original BLAST uses the Poisson method; BLAST2 uses the sum-of scores method.
11. *Show the gapped Smith-Waterman local alignments of the query and each of the matched database sequences.* The original BLAST algorithm only generates ungapped alignments including the initially found HSPs individually, even when there is more than one HSP found in one database sequence. By contrast, BLAST2 produces a single alignment with gaps that can include all of the initially found HSP regions. Note that the computation of the score and its corresponding $E$-value involves use of adequate gap penalties.

The $E$-value is the single most important parameter to rate the quality of the alignments reported by BLAST. For instance, if $E = 10$ for a particular alignment with score $S$, it means that there are $10$ alignments with score $\geq S$ that can happen by chance between any query sequence and the database used for the search. Therefore, this particular alignment is most likely not very significant. By contrast, values much smaller than $1$ (*e.g* $10^{-3}$ or even smaller) are likely to signal that the sequences are homologous. On most webserver, it is possible to filter out all matches that have an $E$-value larger than some threshold. For instance, on the [NCBI webserver](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome), the "expect threshold" defaults to $0.05$.

[^BLAST_wiki]: I took most of this description from [Wikipedia](https://en.wikipedia.org/wiki/BLAST_(biotechnology)#Algorithm), which in my opinion provides one the thorough explanation of the BLAST algorithm
[^DNA_BLAST_words]: For DNA words, common parameters are +5/-4 or +2/-3 for matches and mismatches.
[^BLAST_MN]: It is possible to derive expressions that use effective rather than true sequence lengths, to compensate for edge effects (an alignment that starts near the end of the query or database sequence is likely not to have enough sequence to build an optimal alignment).

Note that the sequence database used to search for matches is preprocessed first, which increases further the overall computational efficiency. With BLAST, the algorithmic complexity of searching the database for a sequence of length $N$ is only $\mathcal{O}(N)$. An online webserver[^BLAST_webserver] to run BLAST can be found [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome).

:::{seealso} Using BLAST
A short tutorial on how to use BLAST from the command line can be found [here](./notebooks/BLAST.ipynb).
:::

[^BLAST_webserver]: There are many!

(sec:MSA)=
### Multiple sequence alignment

Pairwise sequence alignment suffers from some issues that are intrinsic to it, and becomes glaring when trying to align sequences that are distantly related. In particular, all pairwise methods depend in some way or another on a number of parameters (scoring matrix, gap penalties, *etc.*), and it is hard to tell what is the "best" alignment if the method finds multiple alignments with the same score. Moreover, a pairwise alignment is not necessarily informative about the evolutionary relationship (and therefore about possibly conserved amino acids or whole motifs) of the sequences that are compared.

In order to overcome these issues, a number of multiple sequence alignment (MSA) methods have been developed. MSA extends the concept of pairwise protein alignment to simultaneously align three or more sequences, providing a broader view of evolutionary relationships, structural conservation, and functional regions among a group of proteins. Aligning multiple sequences make it possible to identify conserved amino acids that may be critical for protein function or stability and therefore are biologically significant. This, in turn, can help to predict structural features, functional motifs, and evolutionary patterns across species, providing a fundamental tool to build phylogenetic trees or inform the classification of proteins into families.

In principle, MSA can be carried out with the dynamic programming algorithms we already introduced. For instance, the recursive rule of the Needleman-Wunsch algorithm for globally aligning three sequences $S$, $T$, and $U$, which extends Eq. [](#eq:needleman_wunsch), is

$$
F_{i,j, k} = \max
\begin{cases}
F_{i - 1, j, k} + s(S_i, -, -) \\
F_{i, j - 1, k} + s(-, T_j, -) \\
F_{i, j - 1, k} + s(-, -, U_k) \\
F_{i - 1, j - 1, k} + s(S_i, T_j, -)\\
F_{i - 1, j, k - 1} + s(S_i, -, U_k)\\
F_{i, j - 1, k - 1} + s(-, T_j, U_k)\\
F_{i - 1, j - 1, k - 1} + s(S_i, T_j, U_k),
\end{cases}
$$ (eq:MSA)

where the dynamic programming "table" is now three-dimensional, and $s = s(a, b, c)$ is the cost function of aligning $a$, $b$, and $c$, which can be either residues or gaps. It is straightforward to see that the dimensionality of the table grows linearly with the number of sequences $M$, and therefore the algorithmic complexity is $\mathcal{O}(N^M)$. The exponential dependence on $M$ makes this approach practically unfeasible when working with real-world examples. 

Unfortunately, better-performing exact methods, *i.e.* methods that provide the global optimal solution by construction, do not exist (yet). As a result, we have to rely on heuristic methods, which only find local minima. One commonly used approach for multiple sequence alignment is called progressive multiple alignment. Here the prerequisite is that we need to know the evolutionary tree, often called the *guide tree*, that connects the sequences we wish to align, which is usually built by using some low-resolution similarity measure, often based on (global) pairwise alignment. Then, we start by aligning the two most closely related sequences in a pairwise fashion, creating what is known as the seed alignment. Next, we align the third closest sequence to this seed, replacing the previous alignment with the new one. This process continues, sequentially adding and aligning each sequence based on its proximity in the tree, until we reach the final alignment. Note that this is done using a "once a gap, always a gap" rule: gaps in an alignment are not modified during subsequence alignments. These methods have a computational complexity of $\mathcal{O}(M^2)$, which makes it possible to align thousands of sequences.

One of the most commonly used tools to perform MSAs is [Clustal Omega](doi:10.1038/msb.2011.75), which is a [command-line tool](http://www.clustal.org/omega/), but it is also available as a [webserver](https://www.ebi.ac.uk/jdispatcher/msa/clustalo). Clustal Omega builds the guide tree using an efficient algorithm (adapted from [](doi:10.1186/1748-7188-5-21)) which has an algorithmic complexity of $\mathcal{O}(M log M)$, making it possible to generate multiple alignments of hundreds of thousands sequences.

## Threading

```{tip}
Most of the text of this part has been adapted from [](doi:10.1385/1-59259-890-0:921) and @finkelstein2016protein.
```

If pair alignment tools find that a given query sequence is found to share more than 30% of its sequence with another, then it is common to think that a reasonable model for that sequence can be built. By contrast, an alignment yielding a similarity of 20% - 25% could be purely coincidental. In reality, things are more complicated, as it has been shown that proteins with rather high sequence identity could be very differently from a structural point of view ([](doi:10.1073/pnas.95.11.6073)). This and other results show that a single quantity, such as sequence identity, is not enough to determine the 3D similarity between two proteins, and more numbers (such as the length of the chain or of the well-aligned regions) are required to build reliable models. For instance, it makes sense that for a short 50-residue protein a 40% sequence identity would be required to generate a good match, while 25% may be enough for 250 residues. However, note that these numbers have a purely statistical value.

A better approach is to go beyond pairwise alignments by using sequence database searching programs such as BLAST which, as we have seen, provide E-values or similar quantities that estimate the reliability of a sequence match by looking at it in the context of the whole library of sequence scores. In addition, more sophisticated BLAST versions (such as PSI-BLAST) make it possible to obtain good matches with less than 20% sequence identity.

Thanks to these tools, and to the ever-growing number of sequences stored in databases, "simple" database searches are often enough to build a model of the query sequence out of a reliable homolog of known structure. However, if no matches are found, or if an independent confirmation is required, we need methods that do not rely on sequence homology. In this case we recast the problem of protein structure prediction as a problem of choice of the 3D structure best fitting the given sequence among many other possible folds. However, what can be the source of "possible" structures?

While an a priori classification is sometimes possible (see, *e.g.*, [](doi:10.1016/0022-2836(88)90366-X)) or [](doi:10.1016/0014-5793(85)80697-9)), a more practical answer is the [Protein Data Bank (PDB)](https://www.wwpdb.org/) where all the solved and publicly available 3D structures are collected. However, using structures stored in the PDB turns a "prediction" problem into a "recognition" one: a fold cannot be recognized if the PDB does not contain an already solved analog. This limits the power of recognition. Nevertheless, there is an important advantage associated to this procedure: if the protein fold is recognized among the PDB-stored structures, one can hope to recognize also the most interesting feature of the protein, namely its function—by analogy with that of an already studied protein.

Certainly, not all protein folding patterns have been collected in PDB yet; however, it most likely already includes the majority of all the folding patterns existing in nature. This hope, substantiated by [](doi:10.1038/357543a0), is based on the fact that the folds found in newly solved protein structures turn out to be similar to already known folds more and more frequently. Extrapolation shows that perhaps about 1500–2000 folding patterns of protein domains exist in genomes, and we currently know more than half of them (including the majority of the most common folds).

```{figure} figures/threading.png
:name: fig:threading
:align: center
:width: 600px

The basic idea of threading: a sequence is "threaded" through templates extracted from a database of known folds (*e.g.* the PDB), and the resulting structure is assigned an energy. The structure having the lowest energy is taken as the optimal candidate.
```

To recognize the fold of a chain having no visible homology with already solved proteins, one can use various superimpositions of the chain in question onto all examined (taken from an a priori classification or from PDB) 3D folds in search of the lowest-energy chain-with-fold alignment, as sketched in [](#fig:threading). This is called the *threading method*. When a chain is aligned with the given fold, it is threaded onto the fold's backbone until its energy (or rather, free energy) is minimized, including both local interactions and interactions between remote chain regions. The threading alignment allows "gaps" in the chain and in the fold's backbone (the latter are often allowed for irregular backbone regions only). Many different algorithms have been proposed for finding the correct threading of a sequence onto a structure, though many make use of dynamic programming in some form. Indeed, in principle, threading is similar to a homology search; the difference is that only sequences are aligned in a homology search, while threading aligns a "new" sequence with "old" folds.

Being physics-inspired, threading can be a powerful tool for structural biologists and biophysicists. It uses mostly statistics-derived pseudopotentials rather than actual energies, which are based on the contact statistics between amino acids as found in known protein structures[^like_substitution_matrices].

As with any other structure-prediction model, there are some issues with threading, the main ones being:

1. The conformations of the gapped regions remain unknown, together with all interactions in these regions.
2. Even the conformations of the aligned regions and their interactions are known only approximately, since the alignment does not include side-chain conformations (which may differ in "new" and "old" structures). Estimates show that threading takes into account, at best, half of the interactions operating in the protein chain, while the other half remain unknown to us. Thus, again, the protein structure is to be judged from only a part of the interactions occurring in this structure. Therefore, this can only be a probabilistic judgment.
3. The number of possible alignments is enormous (remember Levinthal's paradox?), which makes it very hard to sort out all possible threading alignments and single out the best one (or ones). Here there are many methods that can be employed, one of which being the statistical mechanics of one-dimensional systems, and the related dynamic programming techniques.

```{figure} figures/threading_example.png
:name: fig:threading_example
:align: center
:width: 600px

The 3D structures of (left) chicken histone H5 and (right) replication terminating protein (rtp). Note that the C-terminal helix of rtp is not shown, since it has no analoug in the histone. The root mean-squared deviation between 65 equivalent $C^\alpha$ positions in the two structures is 2.4 $\angstrom$. Taken from [](doi:10.1002/prot.340230311).
```

As an example, @finkelstein2016protein shows the structure prediction done for the replication terminating protein (rtp) by threading. Having threaded the rtp sequence onto all PDB-stored folds, [](doi:10.1002/prot.340230311) used threading with statistics-based pseudopotentials to show that the rtp fold must be similar to that of H5 histone. This a priori recognition, which is shown in [](#fig:threading_example), turned out to be correct.
 
However, it also turned out that the alignment provided by threading deviates from the true alignment obtained from superposed 3D structures of rtp and H5 histone. On the one hand, this shows that even a rather inaccurate picture of residue-to-residue contacts can lead to an approximately correct structure prediction. On the other hand, this shows once again that all the mentioned flaws (insufficiently precise interaction potentials, uncertainty in conformations of nonaligned regions, of side chains, *etc.*) make it possible to single out only a more-or-less narrow set of plausible folds rather than one unique correct fold. Indeed, the set of the "most plausible" folds can be singled out quite reliably, but it still remains unclear which of these is the best. The native structure is a member of the set of the plausible ones, it is more or less close to the most plausible (predicted) fold, but this is all that one can actually say even in the very best case.

Threading methods became a tool for a tentative recognition of protein folds from their sequences. The advantage of these methods is that they formulate a recipe: do this, this and this, and you will obtain a few plausible folds, one of which has a fairly high chance of being correct.

[^like_substitution_matrices]: Similar, in spirit, to the way [substitution matrices](#sec:substitution_matrices) are computed.

## AlphaFold

:::{tip}
The most technical part of this section has been taken/adapted from a [great blog post](https://www.blopig.com/blog/2021/07/alphafold-2-is-here-whats-behind-the-structure-prediction-miracle/) written by [Carlos Outeiral](https://carlos.outeiral.net/).
:::

AlphaFold2 (which I will refer to simply as "AlphaFold" from now on) is a recent deep learning model developed by DeepMind, presented in [](doi:10.1038/s41586-021-03819-2), designed to predict protein structures with remarkable accuracy. The work behind AlphaFold [has been awarded](https://www.nobelprize.org/prizes/chemistry/2024/press-release/) half of the 2024 Nobel Prize in Chemistry. The model's architecture is based on neural networks, and integrates several advanced techniques from machine learning and structural biology to predict the 3D structure of a protein out of its sequence (*i.e.*, the *folding problem*). I will briefly describe the internal architecture of Alphafold, and then show how to use it (in a slightly improved version called ColabFold, presented in [](doi:10.1038/s41592-022-01488-1)).

### Input and preprocessing

```{figure} figures/alphafold.png
:name: fig:alphafold
:align: center

The architecture of AlphaFold. Arrows show how the information flows among the components. The input (preprocessing) module, the evoform and the structure module are highlighted in red, blue and yellow, respectively. Adapted from [](doi:10.1038/s41586-021-03819-2).
```

The overall architecture, as presented in the original paper, is shown in [](#fig:alphafold). The first part, which handles the user input and is highlighted in red in [](#fig:alphafold), is a preprocessing pipeline that can be carried out independently of the rest as done, for instance, in [ColabFold](#sec:colabfold). First of all, the AlphaFold system uses the input amino acid sequence to query several databases of protein sequences, and constructs a [multiple sequence alignment](#sec:MSA) in order to determine the parts of the sequence that are more likely to mutate. The underlying idea is that, if two amino acids (possibly far apart along the sequence) are in close spatial contact, mutations in one of them will be closely followed by mutations of the other, so that the overall 3D structure is preserved. To make an extreme (and somewhat unrealistic) example, suppose we have a protein where an amino acid with negative charge (say, glutamate) is spatially close to an amino acid with positive charge (say, lysine), although they are both far away in the sequence. Most likely, the resulting electrostatic interaction stabilises the structure of the protein. If for evolutionary reasons the first AA mutates into a positively charged amino acid, the second AA will be under evolutionary pressure to mutate into a negatively charged amino acid in order to preserve the electrostatic attraction and therefore the contribution to the stability of the folded protein. The MSA makes it possible to detect this sort of correlations. 

```{figure} figures/myoglobin_templates.png
:name: fig:myoglobin_templates
:align: center
:width: 500

The 3D structure of human myoglobin (top left), african elephant myoglobin (top right, 80% sequence identity), blackfin tuna myoglobin (bottom right, 45% sequence identity) and pigeon myoglobin (bottom left, 25% sequence identity). Credits to [Carlos Outeiral](https://www.blopig.com/blog/2021/07/alphafold-2-is-here-whats-behind-the-structure-prediction-miracle/).
```

AlphaFold also tries to identify proteins that may have a similar structure to the input ("templates"), and constructs an initial representation of the structure called the "pair representation" which is, in essence, a model of which amino acids are likely to be in contact with each other. Finding templates follows a completely different, but closely related principle: proteins mutate and evolve, but their structures tend to remain similar despite the changes. In [](#fig:myoglobin_templates), for example, I display the structure of four different myoglobin proteins, corresponding to different organisms. You can appreciate that they all look pretty much the same, but if you were to look at the sequences, you would find enormous differences. The protein on the bottom right, for example, only has ~25% amino acids in common with the protein on the top left. 

### The Evoformer

For years, the use of MSAs to detect correlations between amino acids relied on statistical analysis, but its limited accuracy required substantial improvements. Early breakthroughs in coevolutionary analysis helped identify biases in the data and correct them using more advanced statistical techniques. These methods contributed to better but still imperfect predictions of protein structures.

In AlphaFold, the information about the MSA and the templates is jointly analysed by a special type of transformer, which is a deep learning method introduced in [](doi:10.48550/arXiv.1706.03762) that uses a mechanism called attention, allowing the model to focus on different parts of the input data, assigning varying importance to specific regions. This architecture became popular in tasks like natural language processing (NLP) due to its ability to model long-range dependencies in sequences. Instead of processing data sequentially, transformers look at all parts of the input simultaneously, which significantly speeds up learning and makes it more efficient. In bioinformatics, this approach has been adapted to analyze MSAs, improving how relationships between residues are understood. AlphaFold takes the application of transformers a step further with its Evoformer architecture (highlighted in blue in [](#fig:alphafold)), which processes both the MSA and the templates, allowing information to flow back and forth between the sequence and the structural representations. Unlike previous deep learning approaches, where geometric proximity was inferred only at the end, AlphaFold continuously updates and refines its structural predictions throughout the process.

Evoformer uses a two-tower transformer model, with one tower dedicated to the MSA (called the MSA transformer) and the other focusing on pairwise residue interactions (pair representation). At each iteration, these representations exchange information, improving both the understanding of the sequence and the predicted pair representation. This iterative refinement, performed 48 times in the published model, allows AlphaFold to adjust its predictions based on both the sequence data and evolving structural hypotheses.

For example, if the MSA transformer identifies a correlation between two residues, it hypothesizes that they are spatially close. This information is passed to the pair representation, which updates the structural model. The pair transformer can then detect further correlations between other residues, perhaps close to the first pair, refining the structure hypothesis. This back-and-forth exchange continues until the model reaches a confident structural prediction.

### The structure module

The information generated by the Evoformer is taken to the the structure module (highlighted in yellow in [](#fig:alphafold)). The idea behind the structure module is conceptually very simple, but its implementation is fairly complex.

The structure module considers the protein as a "residue gas": every AA is modelled as a triangle, representing the three atoms of the backbone. These triangles float around in space, and are moved by the network to form the structure. The transformations are parametrised as 4x4 affine matrices that handle rotations and translations. At the beginning of the structure module, all of the residues are placed at the origin of coordinates. At every step of the iterative process, AlphaFold produces a set of affine matrices that displace and rotate the residues in space. This representation does not reflect any physical or geometrical assumptions, and as a result the network has a tendency to generate structural violations. This is particularly visible in the supplementary videos of [](doi:10.1038/s41586-021-03819-2), that display some deeply unphysical snapshots (see [](#fig:alphafold_structure_module) for a striking example).

:::{figure} figures/alphafold_structure_module.mp4
:name: fig:alphafold_structure_module
:align: center

The evolving prediction of the CASP14 multi-domain target (863 residues) T1091. Individual domains' structure is determined early, while the domain packing evolves throughout the network. The network explores unphysical configurations throughout the process, resulting in the long "strings" that appear during the evolution. Supplementary Video 4 of [](doi:10.1038/s41586-021-03819-2), which can be downloaded [here](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03819-2/MediaObjects/41586_2021_3819_MOESM6_ESM.mp4).
:::

The secret sauce of the structure module is a new flavour of attention devised specifically for working with three-dimensional structures which is based on the very simple fact that the L2-norm of a vector is invariant with respect to translations and rotations[^IPA]. Why is this invariance such a big deal? You may understand this as a form of data augmentation: if the model knows that any possible rotation of translation of the data will lead to the same answer, it will need a lot less data to pull it away from wrong models and will therefore be able to learn much more.

Finally, it is also worth noting that the Structure Module also generates a model of the side chains. To simplify the system, their positions are parametrised by a list of dihedral angles $\chi_1$, $\chi_2$, *etc.* (depending on the AA), which are predicted in their normal form by the network, and implemented with standard geometric subroutines.

[^IPA]: An explanation of this Invariant Point Attention (IPA) is well-beyond my expertise. If you are interested, I refer to  section 1.8 of the Supplementary Information of the original paper.

### Recycling and loss function

AlphaFold leverages a so-called "recycling mechanism": after generating a final structure, it takes all the information (*i.e.* the MSA and pair representations, as well as the predicted structure) and passes it back to the beginning of the Evoformer blocks. In the original paper the whole pipeline is executed three times.

The quality of the training and of the final output of AlphaFold are determined by a specialized loss function called Frame Aligned Point Error (FAPE), a modified version of RMSD, to align atomic positions, which helps preventi incorrect protein chirality. However, this is only part of a more complex loss function, which is a weighted sum of various auxiliary losses. These include losses calculated from multiple iterations of the structure module and a distogram loss, which compares predicted 2D distance matrices with the true structure.

Another interesting aspect is MSA masking, inspired by self-supervised learning models, where some symbols in the MSA are masked and the model is asked to predict them. An additional trick employed is the so-called self-distillation: in this approach, they took a model trained exclusively on the PDB, and predicted the structures of ~300k diverse protein sequences. They then retrained the full model, incorporating a small random sample of these structures at every training cycle. The claim is that this operation allows the model to leverage the large amount of unlabelled data available in protein sequence repositories.

### Output

The main output of AlphaFold is the 3D coordinates of each atom of the final predicted structure. However, another very useful metric that is reported by the model is the predicted local-distance difference test (pLDDT), which provides a per-residue confidence score on a scale from 0 to 100, indicating how certain the model is about the position of each amino acid in the predicted 3D structure. A higher pLDDT score suggests greater accuracy, with scores above 90 being highly reliable, while scores below 70 indicate regions with lower confidence, possibly due to disorder or flexibility.

(sec:colabfold)=
### Using AlphaFold through ColabFold



# Nucleic acids

As for proteins, the structure of nucleic acids can be studied at different levels. The primary structure is the unidimensional list of nucleotides, which can be obtained fairly easily through sequencing. In DNA and RNA, the secondary structure is the list of base pairs that are formed between different parts of the strand(s), while the 3D structure is the description of how the polymer bends and twists. Here, contrary, to proteins, attractive interactions beyond base pairing (which includes stacking) are not very strong, and since nucleic acids are hydrophilic, there is no strong tendency of the nucleotides to tightly pack together. As a result, the free energy that stabilises a molecule's conformation mostly comes from the secondary (rather than tertiary, as in proteins) structure. Therefore, there are many useful things that we can learn about a certain RNA or DNA system by just knowing its secondary structure.

Here I will present methods and algorithms that have been developed for single chains of RNA, but they are mostly valid for (or extendable to) multi-stranded and DNA systems. Formally, given a chain of length $N$, a secondary structure can be represented as a graph where each nucleotide $i$ is a vertex, and each base pair $(i, j)$ is an edge connecting two vertices. The connectivity can be represented  by an adjancency matrix $\hat A$, where each element of the matrix $A_{i,j}$ is equal to $1$ if $i$ and $j$ are linked (either through the backbone or by hydrogen bonds) and $0$ otherwise. In the case of the single chain, nucleotides are ordered and the backbone is continuous, thus $a_{i,i+1} = 1$. Moreover, since each nucleotide can be involved in most one base pair, for each $i$ there is at most one $j$ different from $i - 1$ and $i + 1$ for which $a_{i,j} = 1$.

```{figure} figures/RNA_graph.png
:name: fig:RNA_graph
:align: center
:width: 600px

Possible secondary structures of an RNA strand of sequence `AUCAGCCGUAAGCGGUAACCUUGUUAGGU`, represented as graphs. The points are the nucleotides (*i.e.* the graph's vertices), while lines that connect them are the graph's edges. The black and coloured lines are the exterior and interior edges, respectively. In (a) no lines cross, while in (b) there are intersections between orange and green lines: a pseudoknot is present.
```

Note that secondary structure prediction becomes much easier if we assume that no pseudoknots can form. This is a strong assumption, but it is very convenient from an algorithmic point of view. In the graph/matrix representation, this condition is fulfilled if, for any pair of base pairs $(i, j)$ and $(k, l)$,

* if $i < k < j$ then $i < l < j$;
* if $k < i < l$ then $k < j < l$.

In the graph representation, pseudoknots are present if some edges intersect. Figure [](#fig:RNA_graph) shows the secondary structure (with and without knots) of a short RNA strand, represented as a graph.

Now we want to predict the secondary structure of the RNA, given its sequence. Since here we are dealing with one-dimensional sequences, I hope it will not be surprising to know that dynamic programming can be leveraged to find a solution to this problem, provided that

1. we use a scoring scheme whereby the free-energy contribution that each base pair (or, more generally, "local" secondary structure motif) has on the overall stability of the molecule is additive;
2. we assume that pseudoknots cannot form, so that the RNA can be split into two smaller ones which are independent.

Indeed, under these conditions the solution to the full problem can be built by solving subproblems, which can be done efficiently by applying dynamic programming. The next two sections will describe two algorithms that can be used to obtain the optimal structure, *i.e. the one with the minimum free-energy (MFE).

## Nussinov's algorithm

The first approach I will describe has been introduced by [](doi:10.1137/0135006). To keep it simple, I will use a version of the algorithm in which the (free) energy of an unpaired base is $0$, and that of a base pair is a fixed value, *i.e.* $-1$, regardless of the type of nucleotides involved.

```{figure} figures/nussinov.png
:name: fig:nussinov
:align: center
:width: 600px

The graph representation of Nussinov's recursion relation. Given a subsequence $S_{ij}$, $i$ is either unpaired (first term after the equal sign), or is paired to some other nucleotide $i < k \leq j$.
```

The intuition behind Nussinov's algorithm, shown pictorially in [](#fig:nussinov), is the following: for any subsequence $[i, j]$, $S_{i,j}$, the $i$-th base can either remain unpaired or be paired with some $k$-th base where $i < k \leq j$. If the $i$-th base is unpaired, the (free) energy of $S_{i,j}$, denoted as $F_{i,j}$, simply reduces to the energy of the subsequence $S_{i+1,j}$, or $F_{i+1,j}$. This forms the first term of the Nussinov recurrence relation.

On the other hand, if the $i$-th base pairs with the $k$-th base, then $F_{i,j}$ consists of the energy contribution of this pairing, denoted as $s_{i,k}$, plus the energies of the two resulting subproblems: the subsequence $S_{i+1,k-1}$, represented by $F_{i+1,k-1}$, and the subsequence $S_{k+1,j}$, represented by $F_{k+1,j}$. Finding the optimal $k$ that minimizes this value gives us the second term in the Nussinov recurrence relation.

:::{note} The original algorithm
If you read the original paper, you will notice that the algorithm was conceived a bit differently. For instance, the optimal secondary structure was defined as the structure maximising the number of base pairs rather than the one minimising the free energy. Here I choose to do the latter to make the algorithm a bit more general, since in principle the energy contribution of a base pair can be made to depend on the types of the two nucleotides.
:::

Therefore, the optimal energy of the subsequence $S_{i,j}$ is the minimum of the energy when the $i$-th base is unpaired and the energy when it is paired with the optimal $k$-th base, which can be written explicitly as follows:

$$
F_{i,j} = \min \left\lbrace F_{i+1,j}, \min_{k} \lbrace F_{i+1,k-1} + F_{k+1, j} + s_{i,k} \rbrace  \right\rbrace.
$$ (eq:nussinov)

In this case the dynamic programming table $\hat F$ is of size $N \times N$, and is initialised so that its diagonal entries are set to zero, since a nucleotide cannot bind to itself. The other entries are set iteratively by starting from the bottom-right entry, where $i = j = N - 1$, so that the matrix is progressively filled up from left to right and bottom to top. The final score, representing the optimal solution for the entire sequence, is found in the upper right corner of the matrix, corresponding to the subsequence $S_{0,N-1}$. Since each $(i,j)$ entry requires an $\mathcal{O}(N)$ minimisation, and there are $\mathcal{O}(N^2)$ entries, the total algorithmic complexity is $\mathcal{O}(N^3)$.

As always with dynamic programming methods, in addition to the optimal score we can obtain the secondary structure by using a traceback matrix $\hat K$, which is initialised during the fill-in phase. In particular, given the optimal secondary structure of the subsequence $S_{i,j}$, $K_{ij} = 0$ if $i$ is unpaired, while $K_{ij} = k$, where $k$ is the nucleotide that is paired to $i$, otherwise. A possible recursive traceback function is

:::{code} plaintext
DEFINE K as the traceback matrix K
DEFINE P as the list of base pairs

FUNCTION traceback(i, j)
    k = K(i, j)
    IF i > j
        RETURN
    IF k == 0
        traceback(i + 1, j)
    ELSE
        ADD (i, k) TO P
        traceback(i + 1, k - 1)
        traceback(k + 1, j)
:::

Note that the complexity of the traceback algorithm is $O(N^2)$, since we move through a $N \times N$ matrix and perform a $\mathcal{O}(1)$ operation on each entry we visit.

The Nussinov's algorithm is quite simplistic and comes with several limitations. In its naive implementation, the algorithm does not account for some important factors in RNA folding. Most notably, it does not take into account that, as we already know, stacking interactions between neighboring base pairs are crucial for RNA stability, at least as much as hydrogen bonding. To address this and other limitations, it is fundamental to incorporate biophysical factors into the prediction model. One improvement is to assign energies to structural elements rather than to individual base pairs, so that the total energy of the RNA structure becomes the sum of the energies of these substructures.

(sec:zuker)=
## Zuker's algorithm

Most of the issues mentioned above can be overcome with the algorithm devised by [Zuker and Stiegler](doi:10.1093/nar/9.1.133). The algorithm is based on the idea that RNA folds into a structure that minimizes its free energy, taking into account not only base pairing, but also other factors. In particular, the algorithm assumes that the RNA folding process can be described by an energy model where base pairs, stacking interactions, and loops (hairpin loops, bulge loops, interior loops, and multi-branch loops) contribute to the total free energy of the molecule. The [Turner](doi:10.1093/nar/gkl472) energy model is commonly used in implementations, which provides empirical parameters for these contributions.

```{figure} figures/zuker_graph.png
:name: fig:zuker_graph
:align: center
:width: 600px

Two schematic representations of the secondary structure of a simple RNA molecule. (a) The conventional representation. (b) The same structure, represented as an undirected graph with exterior and interior edges. The legend in the bottom right applies to both panels. Note that here "bifurcation loop" is used in place of multibranched loop. Taken from [](doi:10.1093/nar/9.1.133).
```

The key idea can be fully appreciated by representing an RNA secondary structure as an [undirected graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)#Undirected_graph), drawn as a semicircle with chords connecting some edges. In the following description I borrow heavily from the original paper, which still is (after more than 40 years) one of the best resources on the subject[^read_original_papers]. Here are a few definitions:

* the $N$ nucleotides are the vertices of the graph;
* the $N - 1$ arcs that connect each pair of $i$ and $i + 1$ nucleotides are the *exterior edges* and represent the backbone bonds;
* chords connecting vertices represent base pairs, and therefore connect only complementary nucleotides (C-G, A-U or G-U); these chords are called *interior edges* and have the following properties:
    1. interior edges do not touch, since each nucleotide can be involved in no more than one base pair;
    2. interior edges do not cross, since pseudoknots are not allowed.
* a *face* is a region of the graph that is bounded on all sides by edges.

These definitions are pictorially represented in [](#fig:zuker_graph). The difference with Nussinov's algorithm is that the (free) energy of a structure is associated not with the base pairs, but with the regions that are bounded by the bonds, *i.e.* with the graph faces associated to that structure. Indeed, the possible secondary structure motifs can be classified according to how they can be represented in terms of faces:

* A face bounded by a single interior edge is a hairpin loop, whose length is the number of exterior edges (the number of nucleotides composing the loop is this number minus one);
* A face bounded by two interior edges is further classified into three groups:
    1. if the interior edges are separated by two single exterior edges on both sides, the face is a stacking region (or loop);
    2. if the face is bounded by a single exterior edge on one side, and multiple exterior edges on the other, the face is a bulge loop;
    3. if there are multiple exterior edges on both sides, the face is an interior loop;
* A face that has $k$ ($k > 2$) interior edges is a $k$-multiloop (or a multibranched loop of order $k$).

The energy of a given structure is the sum of the energies associated to each of its faces, $E_P = \sum_F E_F$. This is a powerful way of setting the statistical weight of each motif by using experimentally-determined values: for instance, since hairpins with loops shorter than $3$ are not possible, we can assign the energy $E_F = \infty$ to hairpin loops having fewer than four exterior edges. Any other specific effect (*e.g.* complicated sequence-dependent effect) can be naturally incorporated in this framework.

Thanks to this decomposition, if we know (or can estimate) the energy contributions associated to each secondary structure motif, we can compute the total energy of any given structure. However, enumerating all possible structures to find the one with the lowest (free) energy would be computationally impossible even for rather small sequences. The problem can be overcome by using dynamic programming, with an approach similar (but different) to that of Nussinov.

As before, $S_{i,j}$ is the subsequence $[i, j]$. For each $S_{i,j}$ we define

* $W_{i,j}$ as the minimum free energy of the subsequence;
* $V_{i,j}$ as the minimum free energy of all structures formed by $S_{i,j}$, in which $S_i$ and $S_j$ are paired with each other; if $S_i$ and $S_j$ cannot pair with each other, $V_{i,j} = \infty$.

The constraint on the minimum length of hairpin loops can be enforced by setting $W_{i,j} = V_{i,j} = \infty$ if $j - i \leq 4$. By constrast, if $j - i = d > 4$, $W_{i,j}$ and $V_{i,j}$ can be written in terms of $W_{k,l}$ and $V_{k,l}$ for various pairs $(k, l)$ for which $l - k < d$, *i.e.* in terms of the minimum energies of smaller substructures. The recursive relations for the two matrices, $\hat V$ and $\hat W$, both of size $N \times N$, can be found by considering the possible structures formed by $S_{i,j}$.

```{figure} figures/zuker_recursive_1.png
:name: fig:zuker_recursive_1
:align: center
:width: 600px

Possible substructures for the subsequence $S_{ij}$, constrained to the presence of a $(i, j)$ edge. Adapted from [](doi:10.1093/nar/9.1.133).
```

We start with $\hat V$. First of all, note that this matrix considers only those substructures that have $S_i$ and $S_j$ paired or, in other words, substructures whose graph representation has the $(i, j)$ interior edge. The presence of this additional edge means that the resulting optimal substructure also has an additional face compared to any other substructure we already evaluated. In particular, as shown in [](#fig:zuker_recursive_1),

1. $(i, j)$ can close a hairpin, thus contributing with an energy $F_H(i, j)$;
2. $(i, j)$ can close a face containing exactly two interior edges, with the other edge being $(k, l)$, with $i < k < l < j$. The face contributes an energy $F_L(i, j, k, l)$ that depends on the face which, as before, can be of three types:
    1. $k = i + 1$ and $l = j - 1$: the face is a stacking region;
    2. $k = i + 1$ or $l = j - 1$ (but not both): the face is a bulge region;
    3. $k > i + 1$ and $l < j - 1$: the face is an interior loop.
3. $(i, j)$ can close a face containing $k$ interior edges, with the other $k - 1$ edges being $(i_1, j_1), \ldots , (i_{k - 1}, j_{k - 1})$. This is a $k$-multiloop, and its energy contribution is $F_M(i, j, i_1, j_1, \ldots, i_{k - 1}, j_{k - 1})$.

In all these cases, the $F_H$, $F_L$ and $F_M$ penalty functions are model parameters, and can be estimated by experiment (see the section on [nearest-neighbour models](#sec:NN_models)). Of course, there are multiple substructures in which the $(i, j)$ edge is present. We make sure to select the one with the lowest free energy by first obtaining the optimal substructures corresponding to cases 2. and 3., and then to select the optimal substructure among the three possibilities listed above. This procedure translates to the following recursive relation:

$$
V_{i,j} = \min \begin{cases}
F_H(i, j) \\
F_S(i, j) + V_{i + 1, j - 1}\\
\min_{i < k < l < j} \lbrace F_L(i, j, k, l) + V_{k, l} \rbrace\\
\begin{aligned}
\min_{k, i < i_1 < j_1 < \ldots \ < i_{k-1} < j_{k-1} < j} \lbrace & F_M(i, j, i_1, j_1, \ldots, i_{k - 1}, j_{k - 1}) + \\
& + \sum_{1 \leq l \leq k} V_{i_l, j_l} \rbrace,
\end{aligned}
\end{cases}
$$ (eq:V_recursion)

where $F_S(i, j) \equiv F_L(i, j, i + 1, j - 1)$ is the energy associated to the possibility that the new face is a stacking region, for which we do not have to carry out any minimisation, since we already know that $k = i + 1$ and $l = j - 1$.

:::{note} A simplified model
Strangely enough, the [original paper](doi:10.1093/nar/9.1.133) contains an excellent description of the simplified model where the energy cost of multibranched loops is zero. However, the authors explicitly state that they use a version of the algorithm where this cost is different from zero, but do not describe it.

In the simplified model described in [](doi:10.1093/nar/9.1.133), the contribution due to a multibranched loop is not given by a specific energy function $F_M$, but it is written in terms of the energy of two substructures (which is why Zuker et al use the term "bifurcation loop"). In particular, the contribution of a $k$-multiloop closed by the $(i, j)$ edge is written as

$$
W_{i + 1, k} + W_{k + 1, j - 1}
$$

for some nucleotide $k$ satisfying $i + 1 < k < j - 1$. This is a strong approximation, since it does not assign any penalty to the presence of a multiloop *per se*: its contribution is only due to the energy of the parts that compose it. The updated recursion for $\hat V$ thus reads

$$
V_{i,j} = \min \begin{cases}
F_H(i, j) \\
F_S(i, j) + V_{i + 1, j - 1}\\
\min_{i < h < k < j} \lbrace F_L(i, j, h, k) + V_{h, k} \rbrace\\
\min_{i + 1 < k < j - 2} \lbrace W_{i + 1, k} + W_{k + 1, j - 1} \rbrace
\end{cases}
$$ (eq:V_recursion_simplified)
:::

```{figure} figures/zuker_recursive_2.png
:name: fig:zuker_recursive_2
:align: center
:width: 600px

Possible substructures for the subsequence $S_{ij}$, this time with no constraints. Adapted from [](doi:10.1093/nar/9.1.133).
```

Now we can write down the recursive relation for $\hat W$. There are once again three possibilities, shown schematically in [](#fig:zuker_recursive_2):

1. $S_i$ and $S_j$ (one of them or both) do not form base pairs in the substructure. Since dangling ends do not contribute to the overall energy, $W_{i,j} = W_{i + 1, j}$ or $W_{i, j} = W_{i, j - 1}$.
2. $S_i$ and $S_j$ form a base pair. In this case $W_{i, j} = V_{i, j}$, which we have already computed.
3. $S_i$ and $S_j$ base pair, but not with each other. In other words, there exist interior edges $(i, k)$ and $(l, j)$ (with $i < k < l < j$). In this case there is an "open bifurcation", as the energy can be written in terms of two substructures: $W(i, j) = W(i, k) + W(k + 1, j) = W(i, l - 1) + W(l, j)$.

The last case requires a minimisation over $k$ (or, equivalently $l$), which translates to the following recursive relation:

$$
W_{i,j} = \min \begin{cases}
W_{i + 1, j}\\
W_{i, j - 1}\\
V_{i, j}\\
\min_{i < k < j} \lbrace W_{i, k} + W_{k + 1, j} \rbrace.
\end{cases}
$$ (eq:W_recursion)

:::{note} A simpler (but less obvious) recursive relation for $\hat W$
If the the $W_{i,j}$ entries with $j - i \leq 4$ are initialised to zero rather than $\infty$, Eq. [](#eq:W_recursion) can be rewritten as

$$
W_{i,j} = \min \begin{cases}
W_{i, j - 1}\\
\min_{i \leq k < j - 4} \lbrace W_{i, k - 1} + V_{k, j} \rbrace.
\end{cases}
$$ (eq:W_recursion_simpler)

This can be proven by rewriting the last case of Eq. [](#eq:W_recursion) so that the first term of the $\min$ argument is $W_{i, k - 1}$ rather than $W_{i, k}$, and the minimisation is carried out over the same values of $k$ of Eq. [](#eq:W_recursion_simpler), yielding

$$
\min_{i \leq k < j - 4} \lbrace W_{i, k - 1} + W_{k, j} \rbrace.
$$

First of all, note that, by construction, $j$ is always paired with $k$, so that $W_{k, j} = V_{k, j}$. Moreover, the new relation contains two additional terms compared to the original one: $W_{i + 1, j}$ and $(W_{i, i - 1} + V_{i, j}) = V_{i, j}$, which are equal to the first and third case of Eq. [](#eq:W_recursion), respectively.
:::

The recursive relations for $\hat V$ and $\hat W$ can be used in a dynamic programming code to obtain the optimal secondary structure of the full sequence $S$. Since there is the constraint on the minimal loop length, after initialisation we start filling the matrices from the $(N - 5, N - 1)$ entry, and continue from bottom to top and from left to right. As in Nussinov's algorithm, the energy of the optimal structure of the entire sequence will be stored at the end of the fill-in phase in the $W_{N - 1, 0}$ entry, and the optimal secondary structure itself can be retrieved by using traceback matrices.

What is the algorithmic complexity of the method? In the last case of Eq. [](#eq:V_recursion), we minimise on the order of the multiloop, and on all its interior edges. This operation has an exponential complexity, and therefore completely kills the computational efficiency. There exist several approximations that can be leveraged to obtain a polynomial complexity. The most extreme one is given by Eq. [](#eq:V_recursion_simplified), which completely neglects the energetic penalty of the multiloop. A more realistic approximation is the following:

$$
F_M(i, j, i_1, j_1, \ldots, i_{k - 1}, j_{k - 1}) \approx a + b k + c k',
$$

where $k'$ is the number of unpaired bases within the multiloop, and $a$, $b$ and $c$ are parameters, to be estimated experimentally. This approach makes use of another auxiliary matrix, $\hat M$, whose generic entry $M_{i,j}$ is the optimal energy of $S_{ij}$, with the constraint that $S_i$ and $S_j$ are part of a multiloop. The recursive relation takes into account the possibilities that $i$ and/or $j$ are unpaired, which contributes $c$, are paired with each other, which contributes $b$, or are part of two multi-loop substructures, and therefore

$$
M_{i,j} = \min \begin{cases}
M_{i, j - 1} + c\\
M_{i + 1, j} + c\\
V_{i, j} + b\\
\min_{i < k < j} \lbrace M_{i, k} + M_{k + 1, j} \rbrace.
\end{cases}
$$ (eq:M_recursion)

We can now rewrite the last case of Eq. [](#eq:V_recursion) by making use of the new matrix, obtaining

$$
V_{i,j} = \min \begin{cases}
F_H(i, j) \\
F_S(i, j) + V_{i + 1, j - 1}\\
\min_{i < k < l < j} \lbrace F_L(i, j, k, l) + V_{k, l} \rbrace\\
\min_{i < k < j} \lbrace M_{i + 1, k} + M_{k + 1, j - 1} + a \rbrace.
\end{cases}
$$ (eq:V_recursion_M)

The algorithmic complexity of the four cases are $\mathcal{O}(1)$, $\mathcal{O}(1)$, $\mathcal{O}(N^2)$ and $\mathcal{O}(N)$. Since there are $\sim N^2$ entries, the overall algorithmic complexity is $\mathcal{O}(N^4)$, and the required storage space is $\mathcal{O}(N^2)$, since only $N \times N$ matrices are required. The computational efficiency can be improved by limiting the size of a bulge or interior loop to some value (often taken to be $30$), which brings the complexity of the third case of Eq. [](#eq:V_recursion_M) down to $\mathcal{O}(N^2)$, and the overall algorithmic complexity down to $\mathcal{O}(N^3)$.

A very nice open-source implementation of the Zuker's algorithm can be found [here](https://github.com/Lattice-Automation/seqfold).

## Beyond the MFE: the McCaskill algorithm

The derivations in this section are taken from [](doi:10.1371/journal.pcbi.1006341).

The Nussinov's and Zuker's algorithm find the optimal secondary structure of an RNA strand, defined as the structure that minimise the overall (free) energy. However, the conformations of macromolecules that are in thermal equilibrium with their environment are not fixed: in principle *any* allowed structure can be visited, given enough time (remember the concept of ergodicity and phase space?). In particular, if we define $\mathcal{P}$ as the *structural ensemble*, *i.e.* the ensemble of allowed structures, the probability that a macromolecule has a specific secondary structure $P$ is given by

$$
p(P) = \frac{e^{-\beta F_P}}{\sum_{P' \in \mathcal{P}} e^{-\beta F_{P'}}} \equiv \frac{e^{-\beta F_P}}{Q},
$$

where $F_P$ is the energy of secondary structure $P$, and $Q$ is the partition function. We now look for a recursive relation to compute the partition function of a sequence $S$. Given a subsequence $S_{ij}$, $i$ is either unpaired, or it is paired with the $k$-th nucleotide, where $i < k \leq j$. For the first case, the partition function is that of the smaller subsequence $[i + 1, j]$, $Q_{i + 1, j}$. For the latter case, the $(i, k)$ edge splits the subsequence in two independent subproblems, and the overall partition function is given by the product of their partition functions, $Q_{k + 1, j}$ and $Q_{i + 1, k - 1} q_{i, k}$. Look at [](#fig:nussinov) and you will realise that it is the same splitting! The total partition function will be a sum over all these cases, *viz*

$$
Q_{i,j} = Q_{i + 1, j} + \sum_{i < k \leq j} Q_{k + 1, j} Q_{i + 1, k - 1} q_{i, k},
$$ (eq:mccaskill)

where $q_{kj} \equiv e^{-\beta \Delta G_{k,j}}$ is the statistical weight of the base pair formed by $S_k$ and $S_j$. This relation makes it possible to write down a dynamic programming code that fill the $\hat Q$ matrix with a complexity $\mathcal{O}(N^3)$. As before, the partition function of the sequence is found in the $Q_{0, N-1}$ entry.

:::{note} The connection with Nussinov's algorithm
Here I'm implicitly using a simplified version of the McCaskill algorithm, where the overall energy is only due to base pairing. If you look closely, you will notice that there is a direct connection between Eq. [](#eq:mccaskill) and Eq. [](#eq:nussinov), which is summarised in the following table (adapted from the Appendix C of @finkelstein2016protein):

| MFE | Partition function |
| --- | --- |
| Energy $E$ | Boltzmann factor $e^{-\beta E}$ |
| Minimisation | Summation |
| Summation | Multiplication |

Therefore, it is possible to turn any MFE algorithm into an algorithm that can compute partition functions, provided that there is no overlap between the cases appearing in the recursive relations of the former, *i.e.* that some structures are counted more than once. Indeed, ambiguous decompositions (as they are called) are not a problem when carrying out a minimisation, but they are when doing a summation...
:::

How can we use the information contained in the partition function to obtain the equilibrium (structural) ensemble of the strand? The most granular information we can compute is the probability that two nucleotides $i$ and $j$ are base paired, $p(i, j)$. However, in order to do so we first have to introduce the auxiliary matrix $\hat Q^{\rm bp}$, whose generic entry $Q^{\rm bp}_{i,j}$ stores the partition function of the subsequence $S_{i,j}$, with the constraint that nucleotides $i$ and $j$ form the $(i, j)$ base pair. Its recursive relation is

$$
Q^{\rm bp}_{i,j} = \begin{cases}
Q_{i + 1, j - 1} q_{i,j} & \text{if } i \text{ and } j \text{ are complementary}\\
0 & \text{otherwise} \\
\end{cases}
$$ (eq:mccaskill_bp)

so that the recursive relation of $\hat Q$ (Eq. [](#eq:mccaskill)) becomes

$$
Q_{i,j} = Q_{i + 1, j} + \sum_{i < k \leq j} Q_{k + 1, j} Q^{\rm bp}_{i, k}.
$$

Note that with this formalism it is easier to employ more realistic energy functions by adding contributions (due to interior and multiloops, for instance) to Eq. [](#eq:mccaskill_bp), in analogy with the strategy devised by Zuker. See [here](doi:10.1186/1748-7188-6-3) for additional details.

We can write down the probability $p(i, j)$ recursively, as the sum of two contributions:

1. The probability that the $(i, j)$ edge exists and it is external, *i.e.* that there are no edges $(k, l)$ for which $k < i < j < l$. This is just
$$
p_{\rm ext}(i, j) = \frac{Q_{0, i - 1} Q^{\rm bp}_{i,j} Q_{j + 1, N - 1}}{Q_{0, N-1}}.
$$
2. The probability that the $(i, j)$ edge exists and it is enclosed by other base pairs. This is a sum over all the possible enclosing base pairs $(h, k)$, with $h < i < j < k$, of the probability that the $(h, k)$ edge exists, $p(h, k)$, times the probability that the $(i, j)$ edge is the exterior edge of the $S_{h + 1, k - 1}$ subsequence, which is given by[^denominator_difference]
$$
\frac{Q_{h + 1, i - 1} Q^{\rm bp}_{i,j} Q_{j + 1, k - 1}}{Q_{h + 1, k - 1}}.
$$
All in all, this contribution amounts to
$$
P_{\rm int}(i, j) = \sum_{h < i < j < k} p(h, k) \frac{Q_{h + 1, i - 1} Q^{\rm bp}_{i,j} Q_{j + 1, k - 1}}{Q_{h + 1, k - 1}}.
$$

The total probability that $i$ and $j$ form a base pair is thus

$$
p(i, j) = \frac{Q_{0, i - 1} Q^{\rm bp}_{i,j} Q_{j + 1, N - 1}}{Q_{0, N-1}} + \sum_{h < i < j < k} p(h, k) \frac{Q_{h + 1, i - 1} Q^{\rm bp}_{i,j} Q_{j + 1, k - 1}}{Q_{h + 1, k - 1}}.
$$

Note that the algorithmic complexity of this computation is $\mathcal{O}(N^4)$, since there is the double sum over $h$ and $k$, which has to be carried out for every $i, j$ pair. However, the complexity can be brought down to $\mathcal{O}(N^3)$ by introducing another auxiliary matrix (see *e.g.* [](doi:10.1186/1748-7188-6-3)).

:::{seealso} Python implementation
Head over [here](./notebooks/RNA_folding.ipynb) for a Jupyter notebook containing code implementing the algorithms discussed here. My version of Zuker's algorithm is highly simplified so that the code is as short and readable as possible.
:::

[^read_original_papers]: It is very common for papers that contributed greatly to a field to be also very well thought out and written. This is only partially true for the paper in question (as mentioned below).
[^denominator_difference]: Eq. (10) of [](doi:10.1371/journal.pcbi.1006341) is written slightly differently, since they use the fact that $Q^{\rm bp}_{h, k} = Q_{h + 1, k - 1} q_{hk}$

## Third-party software

You can find several codes online implementing the algorithms we just discussed. However, many of them are a bit rough around the edges and/or incomplete (or wrong). Some interesting and well-done open source implementations are available in different languages. [seqfold](https://github.com/Lattice-Automation/seqfold) implements the Zuker algorithm for both DNA and RNA, and it is a rather readable code written in Python. At [this website](http://rna.informatik.uni-freiburg.de/Teaching/) there are several pages where you can test "RNA algorithms" (several of which we already introduced). Interestingly, the authors (from Freiburg University) also released an excellent [paper](doi:10.1371/journal.pcbi.1006341) detailing their implementations, as well as their [code](https://github.com/BackofenLab/RNA-Playground).

The current state-of-the-art software codes for RNA folding are RNAfold, bundled in the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) package, and mFold, which was originally developed by M. Zuker and is been since merged with other packages in the [UNAFold](http://www.unafold.org/) software. Both softwares can be used online as webservers, but can also downloaded and installed locally. However, ViennaRNA can be downloaded freely, while UNAFold requires a license (although older versions of mFold can be freely downloaded [here](http://www.unafold.org/download/mfold-3.6.tar.gz)).

The best online webserver is surely the one based on ViennaRNA, which can be found [here](http://rna.tbi.univie.ac.at/). It provides tools to do almost everything related to RNA secondary structures, and I invite you to use it (and possibly read the accompanying papers). The tool that implements the algorithms of Zuker and McCaskill is [RNAfold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi). Try it, and compare its predictions to those of the [simple codes](./notebooks/RNA_folding.ipynb) we discussed in class.

### NUPACK

An online webserver, which is also a Python library, that has a slightly different scope but uses algorithms that are very much inline with those we just studied, is NUPACK. The most recent NUPACK version is described [here](doi:10.26434/chemrxiv-2022-xv98l). NUPACK has been free-to-use for many years, but it is apparently on the verge of requiring a paid subscription to use its webservices. By contrast, for now the Python library, which has an open-source-like license for non-commercial academic use, is free-to-download [here](https://nupack.org/download/software), provided that you are a registered user.

I want to briefly describe what NUPACK does, and why it is one of the most used software by the DNA/RNA nanotechnology community. At its core, NUPACK integrates statistical thermodynamics and advanced algorithms to model the folding, hybridization, and thermodynamic stability of nucleic acids. The main appeal of NUPACK compared to other packages lies in its ability to model multi-strand nucleic acid systems. In contrast to simpler single-strand predictions, many DNA nanotechnology applications involve multiple strands interacting to form complex assemblies. NUPACK handles these multi-stranded configurations by extending the core ideas found in the classical algorithms we studied, generalizing them to more complicated systems where the combinatorial possibilities of base-pairing grow exponentially not only with the number of nucleotides, but also with the number of strands.

In practical terms, NUPACK enables researchers to predict the secondary structures that DNA or RNA strands will adopt under given conditions. NUPACK's algorithms minimize free energy to predict the most stable structure, but also consider an ensemble of possible structures, going beyond the minimum free energy configuration. Indeed, is uses McCaskill-like algorithms to compute the partition function, providing information about the entire landscape of possible configurations. Going beyond folding, NUPACK can also design sequences for specific structures: In nanotechnology, designing DNA or RNA that reliably folds into a desired structure is essential for creating functional nanodevices or molecular machines. NUPACK optimizes nucleotide sequences to favor specific target structures while minimizing the formation of undesired alternatives, often through iterative sequence design techniques rooted in statistical optimization.

:::{seealso} Python implementation
Head over [here](./notebooks/NUPACK.ipynb) for some examples on how to use the NUPACK library.
:::

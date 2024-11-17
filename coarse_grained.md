---
title: Coarse-grained models
exports:
   - format: pdf
---
 
# Accuracy

:::{tip}
This part has been adapted (taken almost verbatim in some parts) from [](doi:10.1021/acs.jctc.2c00643). However, in that review accuracy is discussed only in the context of bottom-up CG models. Here I have tried to provide a more inclusive discussion.
:::

How accurate a particular CG model is? Answering this question is hard not because it requires costly calculations, but because there are no metrics that makes it possible to provide an unambiguous answer. What do we mean with "accurate"? Following [](doi:10.1021/acs.jctc.2c00643), I define three separate yet related measures of CG model fidelity: consistency, representability, and transferability. The use and importance of each of these measures are dependent upon the scientific question of interest. It is therefore worthwhile to discuss each of these metrics with the understanding that all three contribute to the overall accuracy of a CG model. Although the original review discussing this topic focusses on bottom-up coarse-graining, I think that some of these metrics can, at least to some extent, extended to top-down models. Indeed, while it is true that top-down models are, by definition, not consistent with statistical mechanics, in the sense that they do not provide a direct connection between the FG and CG worlds but are instead primarily models in the larger sense of the word, they can nevertheless reproduce, sometimes quantitatively, some "macroscopic" properties of the coarse-grained system (*e.g.* thermodynamics).

# Bottom-up

Bottom-up coarse-graining is a robust strategy in molecular modeling that allows detailed atomistic information to be condensed into simplified representations. This approach systematically derives effective potentials and interaction parameters from molecular-level data, reducing complex systems to manageable models while retaining critical physical characteristics. By representing groups of atoms or molecules as single units, bottom-up coarse-graining enables the study of larger systems or longer timescales than would be feasible with fully atomistic simulations, making it especially valuable for exploring mesoscale phenomena.

In practice, bottom-up coarse-graining methods often use information from all-atom simulations to inform the development of coarse-grained potentials. One widely used approach is the iterative Boltzmann inversion, which uses distributions from high-resolution simulations to create effective potentials for coarse-grained models. Other techniques, like the force-matching or multiscale coarse-graining methods, involve more direct optimization to match forces and structural correlations between atomistic and coarse-grained representations. The goal is to achieve a model that can replicate the thermodynamic and structural properties of the original system accurately.

Bottom-up coarse-graining procedures often retain the molecular specificity and accuracy needed for a range of applications, from materials science to biomolecular simulations. However, this accuracy comes at a computational cost, as it requires initial high-resolution data and complex optimization procedures. Nonetheless, by sacrificing only unnecessary degrees of freedom, bottom-up coarse-graining enables researchers to strike a balance between simulation speed and model fidelity, making it an invaluable approach for studying mesoscale phenomena while staying rooted in molecular detail.

A key challenge in bottom-up coarse-graining is that while the coarse-grained model allows for efficient simulation of broader-scale dynamics, it inherently sacrifices fine-grained structural detail. This is where backmapping becomes essential. Backmapping is the process of reconstructing detailed atomic configurations from a coarse-grained representation, effectively transforming a simplified model back into a high-resolution form. This procedure is crucial for bridging scales, as it enables researchers to interpret coarse-grained simulations in terms of atomistic detail when necessary, such as when making predictions about properties that depend on atomic-scale structures.

To achieve accurate backmapping, algorithms often combine statistical inference with additional energetic minimization to produce atomistic configurations that are consistent with the original coarse-grained trajectory. By leveraging data from the coarse-grained simulation as initial constraints, these methods reconstruct missing details in a way that reflects the system's underlying physical behavior. Backmapping is essential for applications where final atomic-scale configurations are needed, such as in the prediction of molecular interactions or the preparation of input structures for further atomistic simulations. Through bottom-up coarse-graining and efficient backmapping, researchers can navigate seamlessly between scales, balancing computational efficiency with molecular accuracy.

# Top-down

## oxDNA/oxRNA

## Vertex models

## A hybrid approach: the Martini force field

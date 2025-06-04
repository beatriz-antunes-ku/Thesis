# Bioinformatics Master Thesis

__Title:__ "Benchmarking prediction of mutational effects on protein-protein interactions on MAVE data"

__Abstract:__

Understanding the impact of missense mutations on protein stability and interactions is
key to revealing their functional effects and potential disease relevance. Distinguishing
variants that disrupt protein function from neutral variants remains a challenge. In
this study, the influence of mutations on protein-protein interactions is analysed using
MAVE data from Weng et al. (2023), which focuses on KRAS and six interaction
partners. AlphaFold3 is used to predict the structure of this protein in monomeric form
as well as a complex. These predictions are then compared to existing crystal structures
to assess the accuracy of computationally predicted models. MAVE experiments are
costly and challenging to perform, leading to limited availability of data. While
computational models offer a viable solution, it is essential to benchmark them with
experimental data to assess the reliability of the results. This would facilitate a more
comprehensive understanding of the mutational effects on a wider range of proteins,
thereby contributing to the advancement of drug development. Rosetta is used to
predict the change in stability upon mutation, and it is executed on two separate
occasions: to estimate the stability of either the monomer or complex, and then to
measure solely the interaction stability. Rosetta scores are compared to MAVE ones
that probe cellular abundance and interaction to evaluate their ability in distinguishing
accurately between wild-type-like and deleterious mutations. The performance of
ESM-1v is evaluated as a comparative alternative to Rosetta for predicting the effects
of mutations. AlphaFold structures are a reliable replacement to crystal structures
as inputs for Rosetta. The interface scoring captures protocol better the effects of
mutations compared to interaction MAVE. This comprehensive approach provides
a detailed evaluation of the strengths and limitations of computational models in
predicting mutation effects on protein-protein interactions.

## Shiny app: [Thiery et al., 2022](https://www.biorxiv.org/content/10.1101/2022.02.15.480567v1)

</br>

### Background

The vertebrate ‘neural plate border’ is a transient territory at the edge of the neural plate containing precursors for all ectodermal derivatives: the neural plate; neural crest; placodes; and epidermis. In our study, we used single cell RNA sequencing to characterise the cellular heterogeneity at the chick neural plate border, as these cell lineages segregate. We classify our cells in an unbiased manner to decipher when key cell states first arise, before modeling transcriptional dynamics as the different ectodermal lineage segregate.

We find that genes previously classified as neural plate border specifiers typically exhibit dynamic expression patterns and are biased towards either placodal or neural crest fates. Instead, given its transient and heterogeneous characteristics, the neural plate border should be defined based on the co-expression of alternative placodal and neural crest gene modules. Finally, we propose a ‘gradient border’ model for cell fate choice at the neural plate border, with the probability of cell fate allocation closely tied to the spatiotemporal positioning of cells.

---
</br>

### Shiny

In this Shiny app, we are providing you with full access to our data which you can explore until your heart is content. In each section, you have the option to select: the full dataset; each of the individual developmental stages sampled; or the neural plate border subset (placodal, neural crest or neural plate border cell states from HH5 onward).

</br>

#### Functionality:

>##### Feature Plots

>Here you can visualise the expression of genes which are variable across at least one of our data subsets.

</br>

>##### Lineage Dynamics

>[scVelo](https://scvelo.readthedocs.io/) and [CellRank](https://cellrank.readthedocs.io/en/stable/index.html) were used to predict the directionality of transcriptional change over time. In doing so we obtained measurements of a given cells latent time and the probability that a cell will differentiate into each of the main ectodermal lineages.

>To model transcriptional dynamics during lineage segregation, we generate generalised additive models (GAMs) of a given genes expression over latent time, weighted by our predicted lineage probabilities. In this section you can visualise the expression of genes of interest across our full dataset, or across our neural plate border subset if you prefer a higher resolution view of placodal/neural crest segregation.

</br>

>##### Co-expression analysis

>One of the key features of the neural plate border is that it does not exclusively express any given set of genes. Instead, we find that previously termed 'neural plate border specifiers' are biased towards either placodal or neural crest lineages. Instead we suggest that the neural plate border is in an unstable state which is better classified based on the co-expression of markers that are later restricted towards alternative fates.

>Here you can visualise the co-expression of any two genes of interest in individual cells. Given that pairs of genes are co-expressed at different levels, we have added the option to threshold out cells with low co-expression of your genes of interest. This allows you to easily see which are the cells that co-expressing genes at the highest levels.

---
</br>

### Reproducibility

All of our analysis is ran using custom Nextflow pipelines and associated Docker containers. We have made use of containers built by the Rocker team - alongside reproducing our analyses, these containers also allow you to interactively explore our data in RStudio server whilst still using the computing environment in which our analysis was carried out.

For full documentation on how to reproduce our analysis visit our [GitHub page](https://github.com/alexthiery/10x_neural_plate_border).


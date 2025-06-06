# clusterMLD: Clustering Multivariate Longitudinal Data

This package is an extension of the original clusterMLD package available at <https://github.com/junyzhou10/clusterMLD>.

The original use -as described below- is untouched.

> Support to efficiently cluster multivariate longitudinal data with sparse (very little observed points for each subject) and irregular (the observation occassions are not aligned across subjects) observations.

> Longitudinal data with multiple outcomes is supported as well, where some of them could be pure noise (non-distinguishable).

> Support the case with potentially unbalanced cluster size, i.e., the number of subjects in some clusters are way outnumbered by the others.

> In sum, the package is capable of clustering **sparse, irregular, unbalanced, and multivariate continuous** longitudinal data

Our main contribution is the function **simLongData()** which generates new samples of the `Longdat` and `Longdat2`. We "reverse engineered" the simulation from the detailed description in the paper.

There is also a new utility plotting function `plot_multivariate_trajectories()` as well as the function `generate_random_curve()`.

## Installation

To install the package:

``` r
devtools::install_github("markusloecher/clusterMLD")
```

## Usage

Use main function `LongDataCluster(x, Y, id, ...)` to cluster longitudinal data in long format. Parallel computing is supported by specifying `parallel = TRUE`.

`DendroPlot(Cluster.object)` yields corresponding dendrogram, where Cluster.object is the output from `LongDataCluster()`

`MeanPlot(Cluster.object)` yields corresponding mean curves of each detected cluster.

## Exploration

The original author created an R shiny app for a better illustration/visualization of the package. Please refer to [clusterMLD_ShinyApp](https://github.com/junyzhou10/clusterMLD_ShinyApp) for more details.

The App is published at <https://junyzhou.shinyapps.io/clusterMLD_ShinyApp/>, with two toy examples uploaded already. Please note that it will take more time to run parallel online than local.

## Limitations

The current version mainly support for the continuous outcomes. Algorithm for other outcome types, such as binary and time-to-event, are now under development.

## Reference

Junyi Zhou, Ying Zhang & Wanzhu Tu (2022) clusterMLD: An Efficient Hierarchical Clustering Method for Multivariate Longitudinal Data, Journal of Computational and Graphical Statistics, [DOI: 10.1080/10618600.2022.2149540](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2149540)

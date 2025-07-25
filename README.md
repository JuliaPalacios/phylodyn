[![CRAN_Status_Badge](http://www.r-pkg.org/badges/last-release/phylodyn)](http://cran.r-project.org/package=phylodyn)

phylodyn
========

This is a branch forked from mdkarcher/phylodyn. The purpose of this branch of `phylodyn` is to incorporate the sampling of genealogies and facilitate phylodynamic inference and analysis from data directly. This branch includes some Python code and installation may be problematic. We are working to resolve this.

## Installation

1. Install (if necessary) package dependencies and helpers `ape`, `spam` and `devtools` using `install.packages`.

2. Install `INLA` using `install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)` 
or check [r-inla.org](http://www.r-inla.org/download) for the most up-to-date installation instructions.

3. Load `devtools` using `library(devtools)`.

4. Install `phylodyn` using

    a. `install_github("JuliaPalacios/phylodyn")`, or

    b. `install_github("JuliaPalacios/phylodyn", build_vignettes = TRUE)` if you want some illustrative vignettes (note: using `build_vignettes = TRUE` will make the install take longer).

## Vignettes

1. [SimpleBNPR](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/SimpleBNPR.Rmd): A short example showing how to use BNPR and BNPR-PS on simulated data, illustraring methodology in [2] and [5].

2. [NewYorkInfluenza](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/NewYorkInfluenza.Rmd): A case study analyzing influenza data from New York, reproducing analysis in [5] on data from [1].

3. [RegionalInfluenza](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/RegionalInfluenza.Rmd): A case study analyzing influenza data from nine geographic regions, reproducing analsyis in [5] on data from [3].

4. [RegionalSeasonality](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/RegionalSeasonality.Rmd): A case study analyzing influenza seasonality from nine geographic regions, reproducing analsyis in [5] on data from [3].

5. [SimplePhyloinfer](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/SimplePhyloinfer.Rmd): A short example comparing BNPR with a split HMC MCMC sampler approach, illustrating methodology in [4].

6. [LongPhyloinfer](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/SimplePhyloinfer.Rmd): A longer example comparing BNPR with multiple MCMC samplers, including split HMC as in SimplePhyloinfer, illustrating methodology in [4].

7. [LocalGenealogies](https://github.com/mdkarcher/phylodyn/blob/master/vignettes/LocalGenealogies.Rmd): A short example of MCMC-based inference of effective population size trajectories from a sequence of local genealogies. Genealogies are assumed to be a realization of the Sequentially Markov Coalescent (SMC') model. The methodology is developed in [6]

8. [KingmanTajimaCounting](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/CountSimulatedCoalescentTrees.Rmd): Example that explains how to count tree spaces conditionally on the data. The methodology is developed in [8]

9. [DistanceTrees](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Distance_RankedGenealogies.Rmd) [10]. 

10. [Tajima estimation](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Tajima_estimation.Rmd) Phylodynamic inference by sampling Tajima tree. Methods developed in [9] and [12]. 

11. [Tajima plotting data](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Tajima_prepare&plotdata.Rmd) Plot the perfect phylogeny of the data.

12. [Adaptive preferential sampling](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Adaptive_prefsamp_INLA.Rmd) Show how to use the adaptive preferential sampling using the INLA approximation developed in [13].

13. [Simulation](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Simulation.Rmd): A short example of simulation of genealogies under isochronous and heterochronous sampling schemes, illustraring methodology in [2] and [5].

14. [Selection](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Parameteric_growth_comparison.Rmd): An application for the hierarchical model described in an upcoming draft to detect different growth rates of the effective population size.

15. [Lineage Tracing Simulation](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/SLT_Simulation.Rmd): Simulation of single cell lineage tracing phylogenies.

16. [Lineage Tracing Estimation]: An application of coalescent single lineage tracing inference. 
    
17. [Multifurcating](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/Multi_Resolution_EPS.Rmd): Simulation of mutifurcating phylogenies from the Lambda coalescent and estimation of Lambda coalescent parameters from this model.
    
18. [Labeled distance](https://github.com/JuliaPalacios/phylodyn/blob/master/vignettes/label_distance.Rmd): Distance between ranked labeled trees.

19. Bounded coalescent. 

## Datasets

Datasets below can be found at: https://github.com/mdkarcher/PhyloData/

1. **New York influenza** BEAST XML for inferring genealogy using sequence data from [1].
    * NewYork.xml

2. **Regional influenza** BEAST XML for inferring genealogy using sequence data from [3].
    * Europe.xml
    * India.xml
    * JapanKorea.xml
    * NorthChina.xml
    * Oceania.xml
    * SouthAmerica.xml
    * SouthChina.xml
    * SoutheastAsia.xml
    * USACanada.xml

## References

1. A. Rambaut, O. G. Pybus, M. I. Nelson, C. Viboud, J. K. Taubenberger, E. C. Holmes
[The genomic and epidemiological dynamics of human influenza A
virus](http://www.nature.com/nature/journal/v453/n7195/full/nature06945.html).
*Nature*, 453(7195): 615–619, 2008.

2. J. A. Palacios and V. N. Minin.
[Integrated nested Laplace approximation for Bayesian nonparametric phylodynamics](http://www.auai.org/uai2012/papers/310.pdf).
In *Proceedings of the Twenty-Eighth International Conference on Uncertainty in Artificial Intelligence*, pages 726–735, 2012.

3. D. Zinder, T. Bedford, E. B. Baskerville, R. J. Woods, M. Roy, M. Pascual.
[Seasonality in the migration and establishment of H3N2 Influenza lineages with epidemic growth and decline](http://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-014-0272-2).
*BMC Evolutionary Biology*, 14(1): 272, 2014.

4. S. Lan, J. A. Palacios, M. Karcher, V. N. Minin, and B. Shahbaba.
[An Efficient Bayesian Inference Framework for Coalescent-Based Nonparametric Phylodynamics](http://bioinformatics.oxfordjournals.org/content/31/20/3282),
*Bioinformatics*, 31(20): 3282-3289, 2015.

5. M. D. Karcher, J. A. Palacios, T. Bedford, M. A. Suchard, and V. N. Minin.
[Quantifying and mitigating the effect of preferential sampling on phylodynamic inference](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004789).
*PLOS Computational Biology*, 12:e1004789, 2016.

6. J.A. Palacios, J. Wakeley,  and S. Ramachandran. [Bayesian nonparametric inference of population size changes from sequential genealogies.](http://www.genetics.org/content/early/2015/07/28/genetics.115.177980) *Genetics* Vol. 201:281-304, 2015.

7. M. Karcher M, J.A. Palacios, S. Lan, V.N. Minin. 
[phylodyn: an R package for phylodynamic simulation and inference](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12630/full), 
Molecular Ecology Resources, 17, 96-100, 2017.

8. L. Cappello and J.A. Palacios. [Sequential Importance Sampling for Multi-Resolution Kingman-Tajima Coalescent Counting](https://arxiv.org/abs/1902.05527), *Annals of Applied Statistics*, in press

9. J.A. Palacios, A. Veber, L. Cappello, Z. Wang, J. Wakeley, S. Ramachandran. [Bayesian Estimation of Population Size Changes by Sampling Tajima's Trees](https://www.genetics.org/content/early/2019/09/11/genetics.119.302373.article-info?versioned=true), *Genetics*, 214 (3) 967-986

10. J. Kim, N.A. Rosenberg, J.A. Palacios. [Distance Metrics for Ranked Evolutionary Trees](https://doi.org/10.1073/pnas.1922851117), *Proceedings of the National Academy of Sciences*, 117(46):28876-28886, 2020.

11. O. Maliet, F. Gascuel, A. Lambert. [Ranked Tree Shapes, Nonrandom Extinctions, and the Loss of Phylogenetic Diversity](https://doi.org/10.1093/sysbio/syy030), *Systematic Biology*, 67(6):1025–1040, 2018.


12. L. Cappello, A. Veber., J.A. Palacios. [The Tajima Heterochronous n-Coalescent: Inference from Heterochronously Sampled Molecular Data](https://www.tandfonline.com/doi/full/10.1080/01621459.2024.2330732), *JASA Applications and Case Studies*, 2024.

13. L. Cappello, J.A. Palacios. [Adaptive Preferential Sampling in Phylodynamics](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010897), *PLOS Comp. Biology*, 2023.

14. J. Zhang, J.A. Palacios. [Multiple merger coalescent inference of effective population size](https://royalsocietypublishing.org/doi/epdf/10.1098/rstb.2023.0306), Philosophical Transactions of the Royal Society B. 

    
## Contributors

Julie Zhang 

Lorenzo Cappello

Vladimir Minin

Michael Karcher

Jaehee Kim

Samyak Rajanala

Shiwei Lan

Wanjing Zhang

Julia Palacios






<!-- README.md is generated from README.Rmd. Please edit .Rmd file -->

<figure>
<img src="./fig/QBRC.jpg" alt="QBRC_logo" />
<figcaption aria-hidden="true">QBRC_logo</figcaption>
</figure>

Please visit our website for more bioinformatics tools:
<https://qbrc.swmed.edu/labs/wanglab>

# SClineager

This is a README file of the R package *SClineager*. In our paper
[available in Cell
Reports](https://doi.org/10.1016/j.celrep.2020.108589), we develop a
Bayesian hierarchical model that performs lineage tracing of single
cells based on genetic markers. Single cell variant calling has two
inherent issues: (1) low coverage on many positions in many cells and
(2) allelic bias due to true monoallelic expression in single cells or
due to sampling bias. This algorithm infers genetic trajectories of
cells by taking these two issues into account. More details about the
structure of the data can be found in the example dataset that goes
along with this R package. The details of the Bayesian model can be
found in our upcoming paper. Detailed usage instructions can be found in
the function manual. Here we provide a basic workflow.

## Configuration

- For this exposition, *Rcpp* ver 1.0.2, *MCMCpack* ver 1.4.5, *vioplot*
  ver 0.3.4, *gplots* ver 3.0.3 are used.

- *SClineager* ver 1.30 is performed in R ver 4.0.3 interfaced with a
  MacBook Pro with a 2.3 GHz Quad Core Intel Core i7 and 32 GB RAM.

## Installation of the package

To install our package, you may simply execute the following codes. The
installation will take about a minute.

``` r
# install.packages("devtools")
devtools::install_github("inmybrain/SClineager", subdir = "SClineager") # don't forget to specify subdir!
```

If you come across a problem like
[this](https://github.com/r-lib/remotes/issues/130), please refer to
[this
answer](https://github.com/r-lib/remotes/issues/130#issuecomment-423830669)
in that issue.

<!-- Or you can install the source file using the command line after downloading it from [here](XXX) (NOT AVAILABLE NOW); -->
<!-- ```{bash, eval = FALSE} -->
<!-- R CMD INSTALL BayesianMIR_1.0.tar.gz -->
<!-- ``` -->

### Dependencies

- To install *SClineager*, you will need the following R packages (and
  their dependencies): *Rcpp* (\>= 1.0.2), *MCMCpack*, *vioplot*,
  *gplots*. This information is specified in DESCRIPTION file in
  ‘’SClineager’’ folder

- No specific operating system is required.

<!-- ## Download the example dataset -->
<!-- The dataset used in the following illustration of _SClineager_ can be downloaded in [here](?). -->

## A basic example of using the package

The names of folders in this example can be anything. Just change the
codes to read from the correct directories. The names and formats of the
.txt files have to be the same. First, your working directory should be
at *test* that contains two sub-folders *mutations* and *processed*.
Please see the following directory structure used in our example. The
*mutations* folder contains the mutation data of different cells (such
as P301_9_5) of different samples (such as 301). The *processed* folder
will need to contain empty sub-folders corresponding to the sample
names.

``` bash
── test
    ├── mutations
    │   ├── 301
    │   │   ├── P301_9_5
    │   │   │   ├── coverage.txt
    │   │   │   └── germline_mutations_hg38.txt
    │   │   ├── P301_9_72
    │   │   │   ├── coverage.txt
    │   │   │   └── germline_mutations_hg38.txt
    │   │   └── P301_9_82
    │   │       ├── coverage.txt
    │   │       └── germline_mutations_hg38.txt
    │   └── 304
    │       ├── P304_8_188
    │       │   ├── coverage.txt
    │       │   └── germline_mutations_hg38.txt
    │       ├── P304_9_40
    │       │   ├── coverage.txt
    │       │   └── germline_mutations_hg38.txt
    │       └── P304_9_63
    │           ├── coverage.txt
    │           └── germline_mutations_hg38.txt
    └── processed
        ├── 301
        └── 304
```

### Read data

We read each mutation information using `read_sclineager`, which takes
about 30 seconds.

``` r
coverage_cutoff <- 3
coverage_percentage <- 0.2
cell_percentage <- 0.2
artefact_percentage <- 0.03

for (folder in c(301, 304)){
  print(folder)
  runinfo <- data.frame(
    Cell = list.files(paste("./mutations/", folder, sep = "")),
    Path = list.files(paste("./mutations/", folder, sep = ""), full = T),
    stringsAsFactors = F
  )
  out_folder <- paste("./processed/", folder, sep = "")
  preprocess_genetics <- read_sclineager(runinfo,
                                        coverage_cutoff,
                                        coverage_percentage,
                                        cell_percentage,
                                        out_folder,
                                        artefact_percentage)
}
```

After running the above code, we get for each mutation a RData file that
will be used in the SClineager model and a pdf file that visualizes data
information. `cleaned.RData` contains an R object named `results`, which
has the following information:

- `mutations_mat` : raw matrix of variants and VAFs
- `runinfo` : same as the `runinfo` parameter in the input, but contains
  only cells found in `mutations_mat` in the same order
- `annotation` : annotation information of the variants in
  `mutations_mat`
- `coverage_mat` : coverage data of the same cells and same variants as
  in `mutations_mat`

``` bash
── processed
    ├── 301
    │   ├── cleaned.RData
    │   └── summary.pdf
    └── 304
        ├── cleaned.RData
        └── summary.pdf
```

### Exploratory analysis

Basic data visualization is performed through `explore_sclineager`,
which takes less than a second.

``` r
folders <- list.files("./processed", full.names = T)
explore_sclineager(folders, coverage_cutoff, "./mCelSeq2_exploratory.pdf")
```

### Run SClineager

The main function `run_sclineager` performs the MCMC sampling and
returns estimation results. The analysis has to be performed for each
sample independently. The whole procedures take less than a few seconds.
The results will be saved as `results.RData` in `file_out`, which
contains an R object named `results` with the following attributes:

- `genotype_mat` : inferred VAF matrix
- `genotype_mat_orig` : the raw VAF matrix
- `sigma` : inferred covariance matrix
- `genotype_mat_all` and `sigma_all`: sampled parameters at each
  iteration
- `runinfo`, `annotation` : same as above

``` r
for (folder in c(301, 304)){
  print(folder)

  file_out <- paste("./processed/", folder, sep = "") # a folder to hold temporary results
  categories <- c("chrM",
                 "synonymous SNV",
                 "nonsynonymous SNV",
                 "UTR5",
                 "UTR3",
                 "nonframeshift substitution",
                 "frameshift substitution")
  max_iter <- 100
  file_in0 <- paste("./processed/", folder, "/cleaned.RData", sep = "")
  mask_genes <- c("HLA-A", "HLA-B", "HLA-C")
  vaf_offset <- 0.01
  dfreedom <- 100
  skip_common <- T
  control <- list(a = 0.5,
                 b = 0.5,
                 c = 1,
                 d = 20,
                 e = 0.4)
  
  sclineager_results <- run_sclineager(
    file_in = file_in0,
    folder = file_out,
    categories = categories,
    max_iter = max_iter,
    keep_genes = "all",
    mask_genes = mask_genes,
    vaf_offset = vaf_offset,
    dfreedom = dfreedom,
    skip_common = skip_common,
    psi = NULL,
    control = control,
    save = FALSE
  )
}
```

As a result, `results.RData` and `imputation_results.pdf` are generated
under each mutation folder.

``` bash
── processed
    ├── 301
    │   ├── cleaned.RData
    │   ├── imputation_results.pdf
    │   ├── results.RData
    │   └── summary.pdf
    └── 304
        ├── cleaned.RData
        ├── imputation_results.pdf
        ├── results.RData
        └── summary.pdf
```

<!-- ### Visualization -->
<!-- Using the fitted model, a scatter plot for multiple instance regression can be provided as follows. -->
<!-- ```{r, eval = TRUE} -->
<!-- MIScatterPlot(tidydata = tidydata,  -->
<!--               bag_size = 5, -->
<!--               true_primary = lapply(1:tidydata$nsample, function(x) rep(c(T,F), c(1, ninst - 1))),  -->
<!--               pred_primary = lapply(split(BMIR_fit$pip[,1], tidydata$membership), function(x) rank(-x, ties.method = "min") <= 1) -->
<!-- ) -->
<!-- ``` -->
<!-- Using slightl modified `ggmcmc::ggs_density` function, we can have the Bayesian inference. -->
<!-- ```{r, eval = TRUE} -->
<!-- # install.packages("ggmcmc") -->
<!-- library("ggmcmc") -->
<!-- ggs_density <- function (D, ncol, family = NA, rug = FALSE, hpd = FALSE, greek = FALSE)  -->
<!--   ## - ncol is added! -->
<!--   ## - ci -> ggmcmc::ci -->
<!--   ## - [Low, High] interval is commented -->
<!-- { -->
<!--   if (!is.na(family)) { -->
<!--     D <- get_family(D, family = family) -->
<!--   } -->
<!--   if (attributes(D)$nChains <= 1) { -->
<!--     f <- ggplot(D, aes(x = value)) -->
<!--   } -->
<!--   else { -->
<!--     f <- ggplot(D, aes(x = value, colour = as.factor(Chain),  -->
<!--                        fill = as.factor(Chain))) -->
<!--   } -->
<!--   f <- f + geom_density(alpha = 0.3) + scale_fill_discrete(name = "Chain") +  -->
<!--     scale_colour_discrete(name = "Chain") -->
<!--   if (!greek) { -->
<!--     f <- f + facet_wrap(~Parameter, ncol = ncol, scales = "free") -->
<!--   } -->
<!--   else { -->
<!--     f <- f + facet_wrap(~Parameter, ncol = ncol, scales = "free",  -->
<!--                         labeller = label_parsed) -->
<!--   } -->
<!--   if (rug)  -->
<!--     f <- f + geom_rug(alpha = 0.1) -->
<!--   if (hpd) { -->
<!--     ciD <- ggmcmc::ci(D) -->
<!--     f <- f + geom_segment(data = ciD, size = 2, color = "blue", inherit.aes = FALSE,  -->
<!--                           aes(x = low, xend = high, y = 0, yend = 0))  -->
<!--     # +geom_segment( -->
<!--     #   data = ciD, -->
<!--     #   size = 1, -->
<!--     #   inherit.aes = FALSE, -->
<!--     #   aes( -->
<!--     #     x = Low, -->
<!--     #     xend = High, -->
<!--     #     y = 0, -->
<!--     #     yend = 0 -->
<!--     #   ) -->
<!--     # ) -->
<!--   } -->
<!--   return(f) -->
<!-- } -->
<!-- ggs_mcmc <- ggmcmc::ggs(BMIR_fit$mcmclist) -->
<!-- ggs_mcmc$Parameter <- factor(ggs_mcmc$Parameter, labels = c(paste0("coef", 1:(nfeature + 1)), "sig2_error")) -->
<!-- ggs_density(ggs_mcmc %>%  -->
<!--               filter(Iteration > ntotal * 1 / 4),  -->
<!--             ncol = 2, -->
<!--             hpd = TRUE) +  -->
<!--   geom_vline(data = data.frame(Parameter = c(paste0("coef", 1:(nfeature + 1)), "sig2_error"), -->
<!--                                true_val = c(rep(2, 1 + nfeature), 1)), -->
<!--              aes(xintercept = true_val), color = "red") + -->
<!--   labs(x = "Value", y = "Density") +  -->
<!--   theme(axis.text.y = element_blank(), -->
<!--         axis.ticks.y = element_blank()) -->
<!-- ``` -->
<!-- ### Prediction in new bags -->
<!-- When new bags (i.e. without labels) are given, we can predict both labels and primary instances using `predict.BMIR`. -->
<!-- ```{r, eval = TRUE} -->
<!-- pred_fit <- predict.BMIR(BMIRchain = BMIR_fit$mcmclist$Chain1,  -->
<!--                          pip = BMIR_fit$pip[,1],  -->
<!--                          tidydata = tidydata,  -->
<!--                          newtidydata = newtidydata,  -->
<!--                          k = 1) -->
<!-- ``` -->
<!-- Let us see how prediction works. -->
<!-- ```{r, eval = TRUE} -->
<!-- ggplot(data = data.frame(pred = pred_fit$newtidydata$label,  -->
<!--                          true = label[-(1:100)]),  -->
<!--        mapping = aes(x = pred, y = true)) +  -->
<!--   geom_point() + geom_abline(intercept = 0, slope = 1, color = "red") -->
<!-- ``` -->

## Another toy example with simulated data

You can reach one of toy datasets used in our simulation study through
[this link](https://github.com/inmybrain/SClineager/tree/master/data),
which points to the folder `data` of this repository. We provide this
small toy dataset because the runtime will be very short.

The main function of *SClineager* is `sclineager_internal`, which runs
the Bayesian sampling with mutations and coverage inputs. The sampling
with this dataset takes about 5-6 seconds.

``` r
load("./data/toy_data.RData")

fit_scl <-
  sclineager_internal(
    mutations_mat = res_data$mutations_mat_obs,
    coverage_mat = res_data$coverage_mat_obs,
    max_iter = 100,
    vaf_offset = 0.01,
    dfreedom = ncol(res_data$mutations_mat_obs),
    psi = diag(1, ncol(res_data$mutations_mat_obs)),
    save = F
  )
```

The following code block is to generate heatmaps for raw and fitted
VAFs.

``` r
library("reshape2") # convert a wide matrix to a long matrix
library("lemon") # to make a shared legend

melted_VAF <- melt(t(fit_scl$genotype_mat))
melted_VAF$Var1 <- factor(x = melted_VAF$Var1,
                          levels = (unique(melted_VAF$Var1))[order_VAF],
                          ordered = TRUE)

gg_fig <-
  ggplot(data = melted_VAF %>%
           filter(Var2 != 71)) +
  geom_tile(aes(x = Var2, y = Var1, fill = value), color = "white") +
  scale_fill_gradient2(
    low = "red",
    high = "white",
    mid = "yellow",
    midpoint = 0.5,
    limit = c(0, 1),
    # space = "Lab",
    name = "Variant Allele Frequency"
  ) +
  labs(x = "Variant", y = "Cell") +
  user_theme +
  theme(axis.text.y = element_blank()) +
  theme(# Change legend key size and key width
    legend.key.size = unit(1.5, "cm"),
    legend.key.height = unit(0.5, "cm"))

### raw
melted_VAF_raw <- melt(t(res_data$mutations_mat_obs))
melted_VAF_raw$Var1 <- factor(x = melted_VAF_raw$Var1, 
                              levels = (unique(melted_VAF_raw$Var1))[order_VAF],
                              ordered = TRUE)
```

The results are compared below.

``` r
grid_arrange_shared_legend(gg_fig + ggtitle("Estimated"), 
                           (gg_fig %+% melted_VAF_raw) + ggtitle("Raw"), 
                           nrow = 1, 
                           position = "top")
```

<figure>
<img src="./fig/toy_VAF.png" alt="heatmap_VAF" />
<figcaption aria-hidden="true">heatmap_VAF</figcaption>
</figure>

## CTCL dataset

We also provide this CTCL dataset from the Mimitou et al. publication
[\[1\]](#1), which we processed and used in our manuscript.

<!-- ## Runtime and memory usage -->
<!-- We perform a series of simulations with different number of cells and variants in these ranges (50, 100, 150, 200 cells and 80, 200, 400, 800, 1200 variants), and measure the runtime and memory used for runnning the internal function `sclineager_internal`. The  length of MCMC chains is 10000. Time is measured in mins. Each combination of the simulation is repeated 100 times to compute an average and a standard deviation.  -->
<!-- |  Time   |J=80             |J=200            |J=400            |J=800             |J=1200             | -->
<!-- |:-----|:----------------|:----------------|:----------------|:-----------------|:------------------| -->
<!-- |I=50  |30.84 (3.159)    |74.23 (2.830)    |77.21 (3.972)    |69.46 (1.642)     |89.97 (1.506)      | -->
<!-- |I=100 |56.95 (0.981)    |140.88 (5.643)   |187.59 (2.263)   |252.49 (2.179)    |382.87 (4.168)     | -->
<!-- |I=150 |98.77 (1.485)    |287.23 (6.868)   |457.03 (5.784)   |726.07 (8.905)    |1125.95 (7.423)    | -->
<!-- |I=200 |184.39 (1.806)   |533.56 (7.774)   |942.97 (12.081)  |1622.80 (16.421)  |2488.86 (32.075)   | -->
<!-- Table: Runtime information for _SClineager_. -->
<!-- Since the memory usage is constant across different repetitions, we show it without a standard deviation. -->
<!-- Memory usage is measured in MB. -->
<!-- | Memory |J=50  |J=100 |J=150 |J=200 |J=250 | -->
<!-- |:------|:-----|:-----|:-----|:-----|:-----| -->
<!-- |I=80   |0.16  |0.22  |0.29  |0.35  |0.42  | -->
<!-- |I=200  |0.53  |0.69  |0.85  |1.01  |1.17  | -->
<!-- |I=400  |1.65  |1.97  |2.29  |2.61  |2.93  | -->
<!-- |I=800  |5.81  |6.45  |7.09  |7.73  |8.37  | -->
<!-- |I=1200 |12.53 |13.49 |14.45 |15.41 |16.37 | -->
<!-- Table: Memory usage for _SClineager_. -->

## References

<a id="1">\[1\]</a> Mimitou, E. P. et al. (2019) Multiplexed detection
of proteins, transcriptomes, clonotypes and CRISPR perturbations in
single cells. Nat. Methods 16, 409–412.

## Notes

- The variant calling of the single cell sequencing data should be
  performed by our variant calling pipeline:
  <https://github.com/tianshilu/QBRC-Somatic-Pipeline> (**temporarily
  closed as of 20.09.24**). One should use the “tumor-only” mode for
  calling mutations, and set “keep_coverage” (keep coverage data) to 1.

- Also, please change the Perl script somatic.pl so that `$skip_recal=1`
  (skip base recalibration), and `$lofreq=1` (use lofreq for scATAC-seq)
  or =0 (use strelka for scRNA-Seq). Alternatively, the users can
  formulate their own variant calling results into the format of our
  pipeline’s results. Please refer to the example datasets incorporated
  in this R package for our format.

- The scSplitter software for splitting 10x Genomics raw fastq files
  into individual cells: <https://github.com/zzhu33/scSplitter>

<!-- - For available covariance structures, see the help page; -->
<!-- ```{r, eval = FALSE} -->
<!-- ?Mclust_SEP_cpp -->
<!-- ``` -->
<!-- - As for initial assignment of cluster membership, each sample is assigned randomly to clusters. -->

## Citation

To cite this package, please use this bibtex format:

``` latex
@article{Lu:2021,
    title = "Overcoming Expressional Drop-outs in Lineage Reconstruction from Single-Cell RNA-Sequencing Data",
    journal = "Cell Reports",
    volume = "34",
    number = "1",
    pages = "108589",
    year = "2021",
    issn = "2211-1247",
    doi = "https://doi.org/10.1016/j.celrep.2020.108589",
    url = "http://www.sciencedirect.com/science/article/pii/S2211124720315783",
    author = "Tianshi Lu and Seongoh Park and James Zhu and Yunguan Wang and Xiaowei Zhan and Xinlei Wang and Li Wang and Hao Zhu and Tao Wang",
    keywords = "scRNA-seq, lineage tracing, drop-out, genetics"
}
```

## Issues

We are happy to troubleshoot any issue with the package;

- please contact to the maintainer by <seongohpark6@gmail.com>, or

- please open an issue in the github repository.

<!-- ## Error and warning messages you may get -->
<!-- ## References  -->

## License

GPL 3.0

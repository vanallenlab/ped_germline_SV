# PedSV helper R library  

This subdirectory contains the source and compiled code for a helper library of R functions that will be useful for various analyses throughout this project.  

## Quickstart  

The library comes pre-installed in one of the project's Docker images, `vanallenlab/pedsv-r`.

```
$ docker run --rm -it vanallenlab/pedsv-r
$ R
> library(PedSV)
```

## Installation  

In an R session:  
```
> install.packages("/path/to/PedSV_0.0.1.tar.gz", repos=NULL, type="source")
```

After installing the package, you can load it as any other R library with `library(PedSV)` or `require(PedSV)`.  

## Some useful functions  

The main purpose of this package is to facilitate analysis and make it as easy as possible to work with the study data in R.  

To that end, the list below includes just a few of helper functions that might come up most often. Note that [help text is available for each of these functions](https://github.com/vanallenlab/ped_germline_SV/tree/main/src#looking-for-more-help) as for any other R functions.  

- `load.sample.metadata()` loads and cleans a .tsv of sample metadata (e.g., `gatk_sv_pediatric_cancers_combined_cohort_metadata_3_31_23.txt`).  
- `load.sv.bed()` loads and cleans a .bed for a single SV callset (e.g., `PedSV.v1.trio_cohort.analysis_samples.wAFs.bed.gz`).  
- `query.ad.matrix()` queries an allele dosage matrix without loading the entire matrix into memory (_see Vignette below_).  
- `query.ad.from.sv.bed()` is a faster, more efficient method to query an allele dosage matrix based on the records present in a BED file loaded by `load.sv.bed()`.  
- `load.constants()` imports various contants, such as colors or labels, for the study (_see Vignette below_).  
- `prep.glm.matrix()` cleans and prepares sample metadata (loaded by `load.sample.metadata()`) as covariates association testing.  
- `get.eligible.samples()` returns lists of eligible case and control samples for association testing.  
- `pedsv.glm()` is a flexible function for conducting case-control association tests in our dataset (_see Vignette below_).  

---  

### Vignette: querying allele dosage matrixes  

The package has helper functions to load allele dosages for one or more SVs without loading the entire allele dosage matrix into memory.  

To accomplish this, use the `query.ad.matrix()` function, which has two optional arguments that will dramatically improve runtime and reduce memory requirements:  

1. If `query.regions` is provided as a list of vectors where each vector is `c(chrom, start, end)`, this function will use [Tabix](http://www.htslib.org/doc/tabix.html) to index into the intervals specified by `query.regions` and only load records overlapping those regions.  

2. If `query.ids` is provided as a vector of SV IDs, this function will only extract data for those SVs.  

For example, this command will only extract allele dosages for SVs `my_SV_ID_1` and `my_SV_ID_2`, provided that those SVs can be found in either `chr1:100-500` or `chr5:30000000-31000000`:  
```
ad.matrix.path <- "/path/to/my/ad.matrix.bed.gz"
query.regions <- list(c("chr1", 100, 500), c("chr5", 30000000, 31000000))
query.ids <- c("my_SV_ID_1", "my_SV_ID_2")

my.ad.results <- query.ad.matrix(path_to_ad_matrix, query.regions, query.ids)
```

If you are only interested in `my_SV_ID_1` and `my_SV_ID_2`, and you know those SVs are present in an SV BED file (like `PedSV.v1.trio_cohort.analysis_samples.wAFs.bed.gz`), you can use the `query.ad.from.sv.bed()` function to abstract away most of the messy nuts and bolts as follows:  

```
# Load SV BED file into memory
my.sv.bed <- load.sv.bed("/path/to/my/sv.bed.gz")

# Subset SV BED file to your variants of interest
my.SVs.of.interest <- c("my_SV_ID_1", "my_SV_ID_2")
subsetted.sv.bed <- my.sv.bed[my.SVs.of.interest, ]

# Query the AD matrix
my.ad.results <- query.ad.from.sv.bed("/path/to/my/ad.matrix.bed.gz", subsetted.sv.bed)
```

Both `query.ad.matrix()` and `query.ad.from.sv.bed()` have several options for returning the results of your query, which can be specified with the `action` argument:  

* `"verbose"` will return the full query matrix. This is the default behavior.  

* `"any"` will return a numeric indicator (0 or 1) if _any_ of the query rows are non-zero for each sample.  

* `"all"` will return a numeric indicator (0 or 1) if _all_ of the query rows are non-zero for each sample.  

* `"sum"` will return the sum of allele dosages for all query rows per sample.  

Finally, note that passing neither `query.regions` nor `query.ids` to `query.ad.matrix()` will load the whole AD matrix into memory and format it for subsequent queries. The output of this basic `query.ad.matrix()` call can be fed to `query.ad.matrix()` again as follows:  

```
## This approach is not recommended unless doing thousands of queries in serial

# Step 1. Load entire AD matrix into memory (note: this will be slow and consume a lot of memory; not recommended)
ad.matrix.path <- "/path/to/my/ad.matrix.bed.gz"
ad.df <- query.ad.matrix(ad.matrix.path)

# Step 2. Query the AD matrix that has already been loaded into memory
query.idx <- c("my_SV_ID_1", "my_SV_ID_2")
my.ad.results <- query.ad.matrix(ad.df, query.ids=query.ids)
```

---  

### Vignette: conducting case-control association tests  

A main purpose of this library is to enable rapid, standardized, and easy association testing between cases and controls.  

This is accomplished with the `pesv.glm()` function, which requires three arguments:  

1. `meta` : a `data.frame` of sample metadata as loaded by `load.sample.metadata()`;  

2. `X` : numeric or boolean/logical independent variable to test for association (_e.g._, SV carrier status); 

3. `Y` : numeric or boolean/logical dependent variable to test for association (_e.g._, case/control labels).  

By default, `pedsv.glm()` runs a basic linear regression model, but you can implement other models with the `family` parameter. See `?family` for more info.  

There are three options for specifying `X` and `Y`:  

1. As an **unnamed vector**. In this case, the values are assumed to be in the same order as the samples in ‘meta’.  

2. As a **named vector**. In this case, the vector names are assumed to be sample IDs, and any sample ID failing to match in ‘meta’ will be dropped.  

3. As a **data frame**. In this case, the row names are assumed to be sample IDs, and any sample ID failing to match in ‘meta’ will be dropped.  

As an example, let's say you want to test for association between Ewing sarcoma diagnosis and samples that carry either of the SVs `my_SV_ID_1` or `my_SV_ID_1`. You could run this test as follows:

```
# Load SV BED file into memory and subset to variants of interest
my.sv.bed <- load.sv.bed("/path/to/my/sv.bed.gz")
my.SVs.of.interest <- c("my_SV_ID_1", "my_SV_ID_2")
subsetted.sv.bed <- my.sv.bed[my.SVs.of.interest, ]

# Load sample metadata
meta <- load.sample.metadata("/path/to/sample/metadata.tsv.gz")

# The independent variable in this analysis is a 0/1 indicator 
# of whether each sample carries either "my_SV_ID_1" or "my_SV_ID_2"
X <- query.ad.from.sv.bed("/path/to/my/ad.matrix.bed.gz", subsetted.sv.bed, action="any")

# The dependent variable is a 0/1 indicator of case/control status
# Here, we can take advantage of the get.eligible.samples() function 
# to get the case and control sample IDs to include in our analysis.
# Note that "EWS" is the study-standardized abbreviation for Ewing sarcoma.
# We will also use the get.phenotype.vector() function to auto-format 
# our dependent variable for us
sample.ids <- get.eligible.samples(meta, "EWS")
Y <- get.phenotype.vector(sample.ids$cases, sample.ids$controls)

# We can now run a logistic regression association test as follows:
results <- pedsv.glm(meta, X, Y, family=binomial())
```

#### A note on covariate adjustment  

The `pedsv.glm` function calls `prep.glm.matrix()`, which transforms `meta` into a matrix of normalized covariates appropriate for association testing.  

By default, these covariates include:  
* One-hot indicator for whether the sample is female (defined as inferred sex label or chrX ploidy > 1.5)  
* Top `N` genetic principal components, where `N` defaults to 5 but can be specified by passing `keep.N.pcs` to `pedsv.glm()` or `prep.glm.matrix()`.   

Other covariates can be added to this default list; see `?prep.glm.matrix` for more details.  

---  

### Vignette: loading helpful project-wide constants  

The package has several built-in constants that are accessible with the `load.constants()` function. See `?load.constants` for full documentation of all options.  

For example, `load.constants("colors")` will load all project-associated colors and palettes, whereas `load.constants("scales")` will load default scale vectors commonly used in SV analyses (for example, log10-scaled genomic distances).  

Most users will probably just want to call `load.constants("all")` at the beginning of most analyses to get all project constants.  

---  

## Looking for more help?  

Can't seem to find what you need? Every function in `PedSV` comes with help text, which you can invoke as for any other R function by prefixing the function name with `?`.  

For example:  
```
> ?load.constants
```  

If the help text still isn't quite helpful enough, get in touch with [Ryan Collins](mailto:Ryan_Collins@dfci.harvard.edu).  

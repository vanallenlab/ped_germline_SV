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
- `load.constants()` imports various contants, such as colors or labels, for the study (_see Vignette below_).  
- `prep.glm.matrix()` cleans and prepares sample metadata as covariates association testing.  

### Vignette: querying allele dosage matrixes  

The package has helper functions to load allele dosages for one or more SVs without loading the entire allele dosage matrix into memory.  

To accomplish this, use the `query.ad.matrix()` function, which has two optional arguments that will dramatically improve runtime and reduce memory requirements:  

1. If `query.regions` is provided as a list of vectors where each vector is `c(chrom, start, end)`, this function will use [Tabix](http://www.htslib.org/doc/tabix.html) to index into the intervals specified by `query.regions` and only load records overlapping those regions.  

2. If `query.ids` is provided as a vector of SV IDs, this function will only extract data for those SVs.  

For example, this command will only extract allele dosages for SVs `my_SV_ID_1` and `my_SV_ID_2`, provided that those SVs can be found in either chr1:100-500 or chr5:30000000-31000000:  
```
ad.matrix.path <- "/path/to/my/ad.matrix.bed.gz"
query.regions <- list(c("chr1", 100, 500), c("chr5", 30000000, 31000000))
query.ids <- c("my_SV_ID_1", "my_SV_ID_2")

my.ad.results <- query.ad.matrix(path_to_ad_matrix, query.regions, query.ids)
```

Note that specifying neither `query.regions` nor `query.ids` will load the whole AD matrix into memory and format it for subsequent queries. The output of this basic `query.ad.matrix()` call can be fed to `query.ad.matrix()` again as follows:  

```
## This approach is not recommended unless doing thousands of queries in serial

# Step 1. Load entire AD matrix into memory (note: this will be slow and consume a lot of memory; not recommended)
ad.matrix.path <- "/path/to/my/ad.matrix.bed.gz"
ad.df <- query.ad.matrix(ad.matrix.path)

# Step 2. Query the AD matrix that has already been loaded into memory
query.idx <- c("my_SV_ID_1", "my_SV_ID_2")
my.ad.results <- query.ad.matrix(ad.df, query.ids=query.ids)
```

Finally, `query.ad.matrix()` has several options for returning the results of your query, which can be specified with the `action` argument:  

* `"verbose"` will return the full query matrix. This is the default behavior.  

* `"any"` will return a numeric indicator (0 or 1) if _any_ of the query rows are non-zero for each sample.  

* `"all"` will return a numeric indicator (0 or 1) if _all_ of the query rows are non-zero for each sample.  

* `"sum"` will return the sum of allele dosages for all query rows per sample.  

### Vignette: loading helpful project-wide constants  

The package has several built-in constants that are accessible with the `load.constants()` function. See `?load.constants` for full documentation of all options.  

For example, `load.constants("colors")` will load all project-associated colors and palettes, whereas `load.constants("scales")` will load default scale vectors commonly used in SV analyses (for example, log10-scaled genomic distances).  

Most users will probably just want to call `load.constants("all")` at the beginning of most analyses to get all project constants.  

## Looking for more help?  

Can't seem to find what you need? Every function in `PedSV` comes with help text, which you can invoke as for any other R function by prefixing the function name with `?`.  

For example:  
```
> ?load.constants
```  

If the help text still isn't quite helpful enough, get in touch with [Ryan Collins](mailto:Ryan_Collins@dfci.harvard.edu).  

# PedSV helper R library  

This subdirectory contains the source and compiled code for a helper library of R functions that will be useful for various analyses throughout this project.  

## Quick-start guide  

### Installation  

In an R session:  
```
> install.packages("PedSV_0.0.1.tar.gz", repos=NULL, type="source")
```

After installing the package, you can load it as with any other R library with `library(PedSV)` or `require(PedSV)`  

#### Example utility: loading helpful project-wide constants  

The package has several built-in constants that are accessible with the `load.constants()` function. See `?load.constants` for full documentation of all options.  

For example, `load.constants("colors")` will load all project-associated colors and palettes, whereas `load.constants("scales")` will load default scale vectors commonly used in SV analyses (for example, log10-scaled genomic distances).  

Most users will probably just want to call `load.constants("all")` at the beginning of most analyses to get all project constants.  

## Looking for more help?  

Can't seem to find what you need? Every function in `PedSV` comes with help text, which you can invoke as for any other R function by prefixing the function name with `?`.  

For example:  
```
> ?load.constants
```  

If the help text still isn't quite helpful enough, get in touch with (Ryan Collins)[mailto:Ryan_Collins@dfci.harvard.edu].  
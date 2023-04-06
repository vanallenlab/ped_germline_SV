# Pediatric Germline SVs
### Structural variant calling and analyses from germline whole-genome sequencing (WGS) in pediatric cancers and controls

Copyright (c) 2023-Present, [Riaz Gillani](RNGILLANI1@partners.org), [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu) and the Van Allen laboratory at Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

#### _Note: this repository is under active development. More documentation will be added as the project evolves._

---  

## Synopsis    

This repository contains the working code and scripts used to discover, genotype, filter, annotate, and analyze germline structural variants (SVs) from WGS across various pediatric cancer and cancer-free control cohorts  

---  

## Table of Contents  

| Directory | Description |  
| :--- | :--- |  
| [`config/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/config) | Environment & Docker configuration files |  
| [`docker/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/docker) | Docker build files |  
| [`gatksv_scripts/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/gatksv_scripts) | Stand-alone scripts for _post hoc_ filtering and manipulation of GATK-SV callsets |  
| [`src/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/src) | Source code for helper packages (see below) |  
| [`sv_filtering/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/gatksv_scripts) | Bash code for filtering the GATK-SV callsets |  
| [`wdl/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/wdl) | Stand-alone WDL workflows |  

---  

### Helper R package  

The `src/` directly contains the pre-compiled source for a helper library of R functions, `PedSV`.  

For more documentation, see the `README` within the `src/` subdirectory.  

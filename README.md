# Pediatric Germline SVs
### Structural variant calling and analyses from germline whole-genome sequencing (WGS) in pediatric cancers and controls

Copyright (c) 2023-Present, [Riaz Gillani](RNGILLANI1@partners.org), [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu), [Jett Crowdis](JettP_Crowdis@DFCI.HARVARD.EDU) and the Van Allen laboratory at Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

---  

## Synopsis    

This repository contains the working code and scripts used to discover, genotype, filter, annotate, and analyze germline structural variants (SVs) from WGS across various pediatric cancer and cancer-free control cohorts.  

For more information, please refer to Gillani*, Collins*, Crowdis, _et al._. Rare germline structural variants increase risk for pediatric solid tumors. [_bioRxiv_](https://www.biorxiv.org/content/10.1101/2024.04.27.591484v1) (2024).  

---  

## Table of Contents  

| Directory | Description |  
| :--- | :--- |  
| [`analysis/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/analysis) | Scripts and other code used for formal analysis |  
| [`config/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/config) | Environment & Docker configuration files |  
| [`docker/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/docker) | Docker build files |  
| [`gatksv_scripts/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/gatksv_scripts) | Stand-alone scripts for _post hoc_ filtering and manipulation of GATK-SV callsets |  
| [`src/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/src) | Source code for helper packages (see below) |  
| [`sv_filtering/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/gatksv_scripts) | Bash code for filtering the GATK-SV callsets |  
| [`version_changelogs/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/version_changelogs) | README for version changelogs |  
| [`wdl/`](https://github.com/vanallenlab/ped_germline_SV/tree/main/wdl) | Stand-alone WDL workflows |  

---  

### Helper R package  

The `src/` directly contains the pre-compiled source for a helper library of R functions, `PedSV`.  

For more documentation, see the `README` within the `src/` subdirectory.  

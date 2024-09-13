# jett_crowdis

This folder serves as a place to document Jett Crowdis's analysis in the germline pediatric SV project.

## Layout

`data`: This folder contains files needed for these analyses, as well as the results of the analyses. 
* In some instances, we cannot give access to underlying data (for example for GMKF RNA files) because their access is controlled by other entities. In these instances, how to obtain the files is detailed in `data/files.txt`
* Large files are not directly committed to the repo. These can be found either as supplementary files or online.
* Some files point to online google bucket paths (for example, the SV files). These files can be obtained according to the manuscript's data access section.
    
`figures`: Figures made by these analyses for use in the manuscript. `manuscript-figures` are figures compiled using illustrator.

`scripts`: Helper scripts for various notebooks.

## Running the notebooks

Each notebook can be run using the virtual environment defined below (with the exception of `cwas-tree-figures.ipynb`, see that notebook for details). Some notebooks must be run before others:

* `preprocess-neuroblastoma-rna` > `sv-impact on rna-generate` > `rna-analysis`
* `define-intergenic-gene-sv-distances` > `generate-cwas-gene-set-enrichment` > `gsea-analysis-figures`
* `generate-cwas-gene-set-enrichment` > `gtex-analysis`

## Python Virtual Environment

To create the virtual environment used to generate this code:

1. Clone this repository using `git`:
   ```sh
   git clone {giturl}
   ```
2. Create conda environment
   ```sh
   conda create --name pediatric-germline-svs-3.7.13 python=3.7.13
   ```
3. Activate the conda environment:
   ```sh
   conda activate pediatric-germline-svs-3.7.13
   ```
   
4. Install dependencies, set environment variables, then connect to the jupyter interface.
   ```sh
   python -m pip install -r jett_analysis/requirements.txt
   conda install ipykernel
   ipython kernel install --user --name=pediatric-germline-svs-3.7.13
   conda deactivate
   ```
   
   The kernel should then be selectable whenever a notebook is opened (see top right of the Jupyter interface)

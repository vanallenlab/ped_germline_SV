# jett_analysis

This folder serves as a place to document Jett's analysis in collaboration on Riaz and Ryan's pediatric SV project.

## Layout

* `lib`: top-level scripts used by various notebooks.
* `ref`: reference files used by various notebooks, with purposes documented in a `README.md`.
* `info`: files containing useful notes. Notably includes `methods.md` and `style.md`, which explain the methods and a style guide for code.

## Python Virtual Environment

I use a virtual environment for reproducibility. To create the virtual environment used to generate this code:

1. Clone this repository using `git`:
   ```sh
   git clone https://github.com/vanallenlab/ped_germline_SV.git
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
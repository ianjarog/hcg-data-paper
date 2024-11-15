# MeerKAT view of Hickson Compact Groups: I. Data description and release
=====

Introduction
------------

This repository contains a set of Python code used to generate the MeerKAT view of Hickson Compact Groups: I. Data description and release paper. 
We use [Snakemake](https://snakemake.readthedocs.io) to organize our scripts into a streamlined and automated workflow that efficiently manages dependencies, 
ensures reproducibility, and allows for scalable execution across various computing environments. Each step of the data analysis is defined in workflow/Snakefile, 
specifying the required inputs, expected outputs, and the set of Snakemake rules to execute.   

Requirements
------------
We have incorporated all necessary softwares into a singularity file that can be downloaded at [Link to be determined]. 

To produce the current paper, all the user has to do is type 'python workflow/run.py' and the pipeline executes each script sequencially and produce the final pdf paper. 
The required data should be put in the data/ folder and they can be retrieved from [Link to be determined].
The folder should be organized as follows

hcg-pipeline/
├── README.md          # Main documentation file for the pipeline
├── config/            # Configuration files and settings
├── data/              # Input data
├── workflow/          # Main directory for the workflow
│   ├── LICENCE        # License file for the workflow
│   ├── Snakefile      # Workflow definition file
│   ├── __init__.py    # Python module initializer
│   ├── __pycache__/   # Python bytecode cache
│   ├── aa.bst
│   ├── aa.cls         
│   ├── cache/         # Cache directory for intermediate data
│   ├── envs/          # Environment definitions for Snakemake (e.g., Conda environments)
│   ├── hcg_paper.tex  # LaTeX document file
│   ├── rules/         # Snakemake rule definitions
│   ├── run.py         # Python script to execute the workflow
│   ├── scripts/       # Custom scripts used in the workflow
│   │   ├── __init__.py
│   │   ├── __pycache__/  
│   │   ├── change_bitpix.py
│   │   ├── compress_pdf.csh
│   │   ├── cut_cube_hcg90_center.csh
│   │   ├── cut_hcg90_center1.txt
│   │   ├── download_decals.py
│   │   ├── download_decals_hcg30_sources.py
│   │   ├── download_decals_hcgs_center.py
│   │   ├── figure_properties.py
│   │   ├── freqtovel.sh
│   │   ├── make_nice_png_hcg31.py
│   │   ├── plot_chanmap.py
│   │   ├── plot_gbt_meerkat_spectra.py
│   │   ├── plot_gbt_vs_meerkat_flux.py
│   │   ├── plot_global_prof.py
│   │   ├── plot_moments.py
│   │   ├── plot_moments_hcg90.py
│   │   ├── plot_moments_hcg90_view_slice.py
│   │   ├── plot_moments_hcg91_delete.py
│   │   ├── plot_noise_specaxis.py
│   │   ├── plot_optical_mom0.py
│   │   ├── plot_segmented_pvdiag.py
│   │   ├── plot_segmented_pvdiag_hcg90.py
│   │   ├── plot_segmented_pvdiag_hcg90_tail.py
│   │   ├── plot_segmented_pvdiag_hcg90_tail_recenter.py
│   │   ├── plot_sofia_products.py
│   │   ├── plot_velocity_ra.py
│   │   ├── run_changebitpix.py
│   │   ├── run_sofia.py
│   │   ├── update_figures.sh
│   │   ├── utility_functions.py
│   │   ├── view_pdfs.sh


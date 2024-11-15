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

.<br />
├── README.md<br />
├── data # folder to store input data<br />
├── config # necessary settings<br />
│   ├── config.yaml<br />
│   ├── figure_properties.yaml<br />
│   ├── figures.yaml
│   ├── figures_backup.yaml
│   ├── parameters.yaml
│   ├── sofia
│   │   ├── hcg16_sofia.par
│   │   ├── hcg30_sofia.par
│   │   ├── hcg31_sofia.par
│   │   ├── hcg90_sofia.par
│   │   ├── hcg91_sofia.par
│   │   └── hcg97_sofia.par
│   └── sofia_pbc
│       ├── hcg16_sofia.par
│       ├── hcg30_sofia.par
│       ├── hcg31_sofia.par
│       ├── hcg90_sofia.par
│       ├── hcg90_sofia_cut_center.par
│       ├── hcg91_sofia.par
│       └── hcg97_sofia.par
└── workflow
    ├── $(document).rb
    ├── ::w
    ├── <h2 id="list-of-make-collocations">List .html
    ├── =7.3
    ├── LICENCE
    ├── Snakefile # a set of snakefile rules to execute
    ├── __init__.py
    ├── __pycache__
    │   └── __init__.cpython-310.pyc
    ├── aa.bst
    ├── aa.cls
    ├── cache
    ├── envs
    ├── hcg_paper.tex
    ├── rules
    │   ├── plot_chanmap.smk
    │   ├── plot_gbt_meerkat_spectra.smk
    │   ├── plot_gbt_vs_meerkat_flux.smk
    │   ├── plot_global_prof.smk
    │   ├── plot_moments.smk
    │   ├── plot_noise_specaxis.smk
    │   ├── plot_optical_mom0.smk
    │   ├── plot_segmented_pvdiag.smk
    │   ├── plot_sofia_products.smk
    │   ├── plot_velocity_ra.smk
    │   ├── run_sofia.smk
    │   └── run_sofia_latest.smk
    ├── run.py
    └── scripts
        ├── :w
        ├── Untitled-1.css
        ├── __init__.py
        ├── __pycache__
        │   ├── __init__.cpython-310.pyc
        │   ├── cubetokms.cpython-310.pyc
        │   ├── figure_properties.cpython-310.pyc
        │   ├── figure_properties.cpython-39.pyc
        │   ├── font_handler.cpython-310.pyc
        │   └── utility_functions.cpython-310.pyc
        ├── change_bitpix.py
        ├── compress_pdf.csh
        ├── cut_cube.csh
        ├── cut_cube_hcg90_center.csh
        ├── cut_hcg90_center1.txt
        ├── download_decals.py
        ├── download_decals_hcg30_sources.py
        ├── download_decals_hcgs_center.py
        ├── download_decals_hcgs_chanmap91.py
        ├── figure_properties.py
        ├── freqtovel.sh
        ├── hcg90_tails.py
        ├── make_nice_png_hcg31.py
        ├── miriad-imcomb.csh
        ├── miriad-imcomb2.csh
        ├── plot_chanmap.py
        ├── plot_gbt_meerkat_spectra.py
        ├── plot_gbt_vs_meerkat_flux.py
        ├── plot_global_prof.py
        ├── plot_moments.py
        ├── plot_moments_hcg90.py
        ├── plot_moments_hcg90_view_slice.py
        ├── plot_moments_hcg91_delete.py
        ├── plot_noise_specaxis.py
        ├── plot_optical_mom0.py
        ├── plot_segmented_pvdiag.py
        ├── plot_segmented_pvdiag_hcg90.py
        ├── plot_segmented_pvdiag_hcg90_tail.py
        ├── plot_segmented_pvdiag_hcg90_tail_recenter.py
        ├── plot_sofia_products.py
        ├── plot_velocity_ra.py
        ├── run_changebitpix.py
        ├── run_sofia.py
        ├── update_figures.sh
        ├── utility_functions.py
        └── view_pdfs.sh

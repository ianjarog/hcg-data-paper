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



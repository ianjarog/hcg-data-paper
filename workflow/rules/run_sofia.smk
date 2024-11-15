import subprocess

rule run_sofia:
    input:
        datacube=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['datacube'],
	parfile=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['sofia_par'],
        datacube_pbc=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['datacube_pbc'],
        weight_cube=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['weight_cube'],
	parfile_pbc=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['sofia_par_pbc']
    output:
        "outputs/sofia/{gal_idx}/{gal_idx}_sofia_updated_vopt_central.par",
        "outputs/sofia/{gal_idx}/{gal_idx}_line60_masked.image_vopt_cat.txt",
        "outputs/sofia/{gal_idx}/{gal_idx}_line60_masked.image_vopt_noise.fits",
        "outputs/sofia/{gal_idx}/{gal_idx}_line60_masked.image_vopt_cubelets/{gal_idx}_line60_masked.image_vopt_1_mom0.fits",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_sofia_updated_vopt_central.par",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_cat.txt",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_mask.fits",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_mom1.fits",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_mom0.fits",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_noise.fits",
        "outputs/data/{gal_idx}/{gal_idx}_line60_masked.image_vopt.fits",
        "outputs/data/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt.fits",
        "outputs/sofia_pbc/{gal_idx}/{gal_idx}_line60_masked.pb_corr_vopt_cubelets/{gal_idx}_line60_masked.pb_corr_vopt_1_mom0.fits"
    priority: 100
    log:
        "outputs/sofia/logs/{gal_idx}/{gal_idx}.log"
    params:
        scfind_threshold=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['scfind_threshold'],
        reliability_minpix=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['reliability_minpix'],
        reliability_threshold=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['reliability_threshold'],
        reliability_enable=lambda wildcards: config['sofia_params'][wildcards.gal_idx]['reliability_enable']
    run:
        # run SoFiA for non-primary beam corrected cube
        command1 = (
            f"python workflow/scripts/run_sofia.py "
            f"--parfile {input.parfile} "
            f"--datacube {input.datacube} "
            f"--weight_cube {input.weight_cube} " # this is not used
            f"--outname outputs/sofia/{wildcards.gal_idx} "
            f"--scfind_threshold {params.scfind_threshold} "
            f"--reliability_minpix {params.reliability_minpix} "
            f"--reliability_threshold {params.reliability_threshold} "
            f"--reliability_enable {params.reliability_enable} | tee {log}"
           )
        # run SoFiA for primary beam corrected cube
        command2 = (
            f"python workflow/scripts/run_sofia.py "
            f"--parfile {input.parfile_pbc} "
            f"--datacube {input.datacube_pbc} "
            f"--weight_cube {input.weight_cube} "
            f"--outname outputs/sofia_pbc/{wildcards.gal_idx} "
            f"--scfind_threshold {params.scfind_threshold} "
            f"--reliability_minpix {params.reliability_minpix} "
            f"--reliability_threshold {params.reliability_threshold} "
            f"--reliability_enable {params.reliability_enable} | tee {log}"
           )
        subprocess.run(command1, shell=True, check=True)
        subprocess.run(command2, shell=True, check=True)
       

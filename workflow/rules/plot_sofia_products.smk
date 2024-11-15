import yaml
import glob
# Load figure names from the YAML file
with open("config/figures.yaml", 'r') as f:
    figures = yaml.safe_load(f)

figure_files = ["outputs/publication/figures/" + fig for fig in figures['figure1']]

rule plot_sofia_products:
    input:
        ngc1622_mom0 = config['sofia_params']['hcg30']['ngc1622_mom0'],
        ngc1622_catalog = "outputs/sofia_pbc/hcg30/hcg30_line60_masked.pb_corr_vopt_cat.txt",
        cubeletes = (f"outputs/sofia_pbc/hcg30/hcg30_line60_masked.pb_corr_vopt_cubelets/hcg30_line60_masked.pb_corr_vopt_1_mom0.fits" 
                    for hcg30 in config['sofia_params'].keys()),
    output:
        "outputs/publication/figures/ngc1622_column_density.pdf",
        "outputs/publication/figures/ngc1622_snr.pdf",
        "outputs/publication/figures/ngc1622_mom1.pdf",
        "outputs/publication/figures/ngc1622_chan.pdf",
        "outputs/publication/figures/ngc1622_mom2.pdf",
        "outputs/publication/figures/ngc1622_spec.pdf",
        "outputs/publication/figures/ngc1622_pv.pdf",
        "outputs/publication/figures/ngc1622_pv_min.pdf"
    priority: 99
    params:
        source_id = lambda wildcards: config['sofia_params']['hcg30']['source_id']
    shell:
        "python workflow/scripts/plot_sofia_products.py"
        " --filetoplot {input.ngc1622_mom0}"
        " --source_id {params.source_id}"
        " --catalogue {input.ngc1622_catalog}"
        

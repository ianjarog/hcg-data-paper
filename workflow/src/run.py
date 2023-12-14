import subprocess

def run_snakemake():
    try:
        # Run the Snakemake workflow
        snakefile_path = 'workflow/src/Snakefile'
        subprocess.run(["snakemake", "--snakefile", snakefile_path, "--use-conda", "--conda-frontend", "mamba", "--default-resources", "--cores", "all"], check=True)
    except subprocess.CalledProcessError as e:
        print("Error running Snakemake:", e)
    else:
        print("Snakemake workflow executed successfully.")

if __name__ == "__main__":
    run_snakemake()


import sys
import os
import argparse
import subprocess
from shutil import which
from astropy.io import fits

def get_args():
    '''This function parses and returns arguments passed in'''
    description = 'Select dataset'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-p', '--parfile', dest='parfile', \
                        help='Link to Sofia-2 paramenters file')
    parser.add_argument('-o', '--outname', dest='outname', help='Name of output\
                        directory for Sofia products')
    parser.add_argument('-d', '--datacube', dest='datacube', help='Data cube to\
                        process.')
    parser.add_argument('-s', '--scfind_threshold', dest='scfind_threshold',
                        help='Define scfind.threshold parameter', default=5)
    parser.add_argument('-f', '--reliability_minpix', dest='reliability_minpix',
                        help='Define reliability.minPixels parameter', default=0)
    parser.add_argument('-t', '--reliability_threshold', dest='reliability_threshold',
                        help='Define reliability.threshold parameter', default=0.9)
    parser.add_argument('-e', '--reliability_enable', dest='reliability_enable',
                        help='Define reliability.enable parameter', default=False)
    parser.add_argument('-w', '--weight_cube', dest='weight_cube',
                        help='Weight cube for SoFiA', default=False)
    args = parser.parse_args()
    return args

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None

def squared_cube(cube, output_cube):
    with fits.open(cube) as hdul:
        data_cube = hdul[0].data
        hdr = hdul[0].header
    # Square the data
    squared_data_cube = data_cube * data_cube
    fits.writeto(output_cube, squared_data_cube, hdr, overwrite=True)
    return output_cube

def update_parfile(parfile, output_path, datacube,
              scfind_threshold, reliability_minpix,
              reliability_threshold, reliability_enable, weight_cube):
    '''Updates file with paramenters

    Parameters
    ----------
    parfile: str
        File contanining sofia parameters
    output_path: str
        Path of output file
    datacube: str
        Path to datacube
    scfind_threshold: float
        Sofia parameter scfind.threshold
    reliability_threshold: float
        Sofia parameter reliability.threshold
    reliability_minpix: int
        Sofia parameter reliability.minPixels
    reliability_enable: str
        Sofia parameter reliability.enable
    Returns
    -------
    updated_parfile: str
        Path of file with updated parameters
    '''
    datacube_path = datacube
    datacube_name = os.path.basename(datacube)
    updated_parfile = os.path.join(output_path, datacube_name.split('_')[0] + '_sofia_updated_vopt_central.par')
    print(f'Updated parfile: {updated_parfile}')
    with open(parfile, 'r') as filein, open(updated_parfile, 'w') as fileout:
        lines = filein.read().replace('output_path', output_path)
        lines = lines.replace('datacube', datacube_path)
        lines = lines.replace('outname', f'{datacube_name}')
        lines = lines.replace('scfind_threshold', scfind_threshold)
        lines = lines.replace('reliability_minpix', reliability_minpix)
        lines = lines.replace('reliability_threshold', reliability_threshold)
        lines = lines.replace('reliability_enable', reliability_enable)
        if "pbc" in updated_parfile:
            lines = lines.replace('weight_cube', weight_cube) 
        fileout.write(lines)
    return updated_parfile

def eliminate_time(cat):
    '''Eliminates timestamp from sofia catalog. Updates the file

    Parameters
    ----------
    cat: str
        Path to sofia catalog
    '''
    # Read in the file
    with open(cat, 'r') as infile :
        lines = infile.readlines()
    # Replace the target string
    for i, line in enumerate(lines):
        if line[:7] == '# Time:':
            lines[i] = '# Time:\n'
    # Write the file out again
    with open(cat, 'w') as infile:
        for i in lines:
            infile.write(i)

def miriad_velsw(cube, out_path = "./outputs/data"):
    """
    Convert frequency axis to velocity axis using MIRIAD task velsw 

    Parameters
    ----------
    cube: str
         Input data cube whose spectral axis in frequency needs to converted to velocity
    out_path: str
        Output directory path to save results to
    """
    cube_name = os.path.basename(cube)
    dir_path = out_path + "/" + cube_name[:5]
    cube_out = dir_path + "/" + cube_name[:-5] + "_vopt.fits"
    os.makedirs(dir_path, exist_ok=True)
    os.system("chmod a+x ./workflow/scripts/freqtovel.sh")
    command = ['./workflow/scripts/freqtovel.sh', 
                cube, cube[:-5], cube_out]
    subprocess.run(command, check=True)
    return cube_out

def execute_sofia(executable, output_catalog, updated_parfile):
    """
    Running sofia.
 
    Parameters
    ----------
    executable: str
        Name of executable, either sofia or sofia2, 
        depends on how SoFia-2 environment was set up
    output_catalog: str
        Name of sofia output catalogue
    updated_parfile: 
        updated sofia parameter file
    """
    print('Executing Sofia-2')
    subprocess.call([executable, updated_parfile])
    print("Catalogue ", output_catalog[:-4]+".xml")
    #subprocess.run(["sofia_image_pipeline", "-c", output_catalog[:-4]+".xml"])
    eliminate_time(output_catalog)

def run_sofia(parfile, outname, datacube, 
              scfind_threshold, reliability_minpix,
              reliability_threshold, reliability_enable, weight_cube):
    """Only runs Sofia if the output catalog  does not exist

    Parameters
    ----------
    parfile: str
        File contanining parameters
    outname: str
        Name of output file
    datacube: str
        Path to data cube
    scfind_threshold: float
        Sofia parameter scfind.threshold
    reliability_minpix: float
        Sofia parameter reliability.minPixels
    reliability_threshold: float
        Sofia parameter reliability.threshold
    reliability_enable: str
        Sofia parameter reliability.enable
    weight_cube: str
        Sofia parameter input.weights
    """
    #It makes sense to not run this when the results exist but maybe a check
    #on an existing catalog is better
    print('outname ', outname, datacube)
    output_path = outname
    # convert data cube in spectral axis to velocity axis using Miriad
    print("DATACUBE", datacube)
    datacube_or = datacube
    datacube = miriad_velsw(datacube)
    #datacube = datacube[:-5] + "_vopt.fits"
    weight_cube = miriad_velsw(weight_cube)
    weight_cube_squared = squared_cube(weight_cube, weight_cube[:-5]+"_squared.fits")
    datacube_name = os.path.basename(datacube.replace('.fits',''))
    output_catalog = os.path.join(output_path, f'{datacube_name}_cat.txt')
    if True: #not os.path.isfile(output_catalog):
        print(f'Parfile: {parfile}')
        updated_parfile = update_parfile(parfile, output_path, datacube,
              scfind_threshold, reliability_minpix,
              reliability_threshold, reliability_enable, weight_cube_squared)
        if is_tool('sofia'):
            execute_sofia('sofia', output_catalog, f"{updated_parfile}")
        elif is_tool('sofia2'):
            execute_sofia('sofia2', output_catalog, f"{updated_parfile}")
        else:
            print('SoFia-2 not available. Please install it.')
            sys.exit(1)
    else:
        print(f"We have already found the catalogue {output_catalog}. \
                Sofia will not be executed" )

def main():
    ''' Run Sofia if the output catalog does not exist'''
    args = get_args()
    if not os.path.isdir(args.outname):
        os.makedirs(args.outname)
    run_sofia(parfile=args.parfile, 
              outname=args.outname,
              datacube=args.datacube, 
              scfind_threshold=args.scfind_threshold,
              reliability_minpix=args.reliability_minpix,
              reliability_threshold=args.reliability_threshold,
              reliability_enable=args.reliability_enable,
              weight_cube=args.weight_cube)

if __name__ == '__main__':
    main()

#python workflow/scripts/run_sofia.py --parfile config/sofia_pbc/hcg90_sofia.par --datacube data/hcg90/hcg90_line60_masked.pb_corr.fits --weight_cube data/hcg90/hcg90_line60_masked.pb.fits --outname results/sofia_pbc/hcg90/ --scfind_threshold 4 --reliability_minpix 0 --reliability_threshold 0.95 --reliability_enable true
#python workflow/scripts/run_sofia.py --parfile config/sofia_pbc/hcg97_sofia.par --datacube data/hcg97/hcg97_line60_masked.pb_corr.fits --weight_cube data/hcg97/hcg97_line60_masked.pb.fits --outname results/sofia_pbc/hcg97/ --scfind_threshold 4 --reliability_minpix 0 --reliability_threshold 0.95 --reliability_enable true

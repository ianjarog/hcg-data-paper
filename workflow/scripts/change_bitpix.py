import numpy as np
from astropy.io import fits
import argparse

#results/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_cubelets
def changebitpix(input, output):
    hdu = fits.open(input)
    hdr = hdu[0].header
    img = hdu[0].data.astype(np.int32)
    img = hdu[0].data.astype(np.int32)
    hdrc = hdr.copy()
    hdrc["BIT32"] = 32
    fits.writeto(output, img, hdrc, overwrite=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FITS file to 32-bit floating point format.")
    parser.add_argument("--input", help="Input FITS file", required=True)
    parser.add_argument("--output", help="Output FITS file", required=True)
    
    args = parser.parse_args()
    
    changebitpix(args.input, args.output)
#python change_bitpix.py --input ../../results/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_cubelets/hcg90_line60_masked.pb_corr_vopt_9_cube.fits --output ../../results/sofia_pbc/hcg90/hcg90_line60_masked.pb_corr_vopt_cubelets/hcg90_line60_masked.pb_corr_vopt_9_cube_bitpix.fits

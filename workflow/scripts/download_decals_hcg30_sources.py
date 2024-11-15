import os
import legacystamps
from analysis_tools.functions import radec2deg

hcg_position = {
        "hcg30": {"ra": "4:36:37", "dec": "-3:11:05"},
}

band = ['g', 'z', 'i','r']
radec = radec2deg(hcg_position["hcg30"]["ra"], hcg_position["hcg30"]["dec"])
ra = round(radec["ra"],6)
dec= round(radec["dec"],6)
print(ra, dec, " RADEC")

for gal in hcg_position.keys():
    for filter in band:
        legacystamps.download(ra=ra, dec=dec, mode='fits', bands=filter, 
                size=0.1833, layer='ls-dr10', pixscale=0.5)
        filename = ("legacystamps_" + str(ra) + '_' 
                + str(dec) + "_ls-dr10.fits")
        os.system('mv ' + filename + " data/" + gal + "/optical/" + "ngc1622" + "-" + filter + "filter.fits")

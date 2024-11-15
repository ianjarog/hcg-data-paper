import os
import legacystamps
hcg_position = {
    "hcg16": {"ra": 32.38881, "dec": -10.16298},
    "hcg30": {"ra": 69.11918, "dec": -2.83239},
    "hcg31": {"ra": 75.40372, "dec": -4.25671},
    "hcg90": {"ra": 330.52343, "dec": -31.96680},
    "hcg91": {"ra": 332.30172, "dec": -27.77593},
    "hcg97": {"ra": 356.86224, "dec": -2.30542}
}

band = ['g', 'z', 'i', 'r']
band = ['r']
pixscale = 2.4
size = 1.8
for gal in ['hcg97']:#hcg_position.keys():
    for filter in band:
        legacystamps.download(ra=hcg_position[gal]["ra"], dec=hcg_position[gal]["dec"], 
             mode='fits', bands=filter, size=size, layer='ls-dr10', pixscale=pixscale)
        if gal != "hcg90":
            filename = ("legacystamps_" + str(hcg_position[gal]["ra"]) + '0_' 
                + str(hcg_position[gal]["dec"]) + "0_ls-dr10.fits")
        else:
            filename = ("legacystamps_" + str(hcg_position[gal]["ra"]) 
                        + '0_' + str(hcg_position[gal]["dec"]) + "00_ls-dr10.fits")
        os.system('mv ' + filename + " data/" + gal + "/optical/" + filename[:-5] + "-" + filter + "filter.fits")

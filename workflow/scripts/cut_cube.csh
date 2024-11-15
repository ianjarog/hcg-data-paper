fits in=$1.fits op=xyin out=$1
imsub in=$1 out=$1_imsub2 region='boxes('$2,$3,$4,$5')('$6', '$7')'
fits in=$1_imsub2 op=xyout out=$1_imsub2.fits
rm -r $1_imsub2
rm -r $1

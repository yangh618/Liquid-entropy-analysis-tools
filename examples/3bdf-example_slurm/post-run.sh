#!/bin/bash

#  Post-analysis
# ============================================

if [ -f 3bdf_bins.dat-sum ];
then
    rm 3bdf_bins.dat-sum
fi

3bdfxdat_sum 3bdf_bins.dat-sum 3bdf_bins.dat_*
rm 3bdf_bins.dat_*		# check your file names to avoid mis-deleting!!!
3bdf2s3 3bdf_bins.dat-sum 2
3bdf2s3_interp.py > tot.s3

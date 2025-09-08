#!/bin/bash

# For this run, POSCAR, POTCAR, XDATCAR are required.

#  Pre-run
# ============================================
cat XDATCAR | pdfxdat_cubic 1> /tmp/Nul
pdf2s2_v2 -x XDATCAR > tot.s1

#  Execute [MODIFY COMPLETELY TO YOUR NEEDS]
# ============================================

# -----------------------------------------------
# define your histogram
rmin=0.7			# maximum distance of non-zero pdf based on pdf.XX
rmax=2.3			# recommend to be 1/4 of the box size
dr=0.02				# bin size, recommend to have 100~500 bins
# -----------------------------------------------

# -----------------------------------------------
# divided "nsample" XDATCAR into "njobs" jobs each has "nxdat" files
ninit=0			# remove the first ninit number of  inital structures
nsamples=100		# total number of structures to calculate
njobs=4 		# devided nsample in to njobs
nxdat=25		# each job deal with nxdat structures
# -----------------------------------------------

date > RunAt
hostname >> RunAt
echo running in `pwd` >> RunAt

init_bins.sh $rmin $rmax $dr
cat XDATCAR | split_xdat.awk	# adjust the name "XDATCAR" to your XDATCAR files

status=3bdfxdat.status		# log file of the run
rm $status;
touch $status;

for i in `seq 0 1 $njobs`;	# a queueing system for running jobs in parallel
do
    if [ $i -eq $njobs ];
    then
	break;
    fi
    
    nstart=$(( $i*$nxdat+1 +$ninit ))
    nend=$(( ($i+1)*$nxdat +$ninit ))
    echo $nstart $nend >> 3bdfxdat_bundle.sh.out
    3bdfxdat_bundle.sh $nstart $nend $nstart >> $status &
done

while true;			# waiting all jobs to finish
do
    sleep 3;
    ans=`grep DONE $status | wc -l`;
    if [ $ans -eq $njobs ];
    then
	break;
    fi
done

echo Done `date` >> RunAt
sleep 10			# necessary to save all results into files when these files are large in memory

rm run_dir* -r			# check your file names to avoid mis-deleting!!!


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

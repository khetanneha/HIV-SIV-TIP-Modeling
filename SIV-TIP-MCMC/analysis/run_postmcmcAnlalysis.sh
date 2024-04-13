
samNo=$1
nRuns=$1 #00

opath='./outA1M1/'
ipathPD=opath
opathPDp='./outA1M1/PDplots/'
opathPDf='./outA1M1/PDfiles/'
dattype='A1M1'


for iS in $( eval echo "{1..$samNo}" )
	do
	
		ipath="/Volumes/NKhd/V24/Jul10.23/dat_070423/A1M1/s${iS}/output/plots/"



		echo  $iS $nRuns  $ipath  $opath
		python poolanalysisnRuns_V10.py $iS $nRuns  $ipath  $opath 17500 20

		# computes meadian and percentiles / Mean and SD
		python plot_PDparam_FitCurve.py $opath $opathPDp $opathPDf $iS $dattype 


		# Plot for KS_test on LogProb scores
		python mcmc_statsLogProbScores.py $opath $opathPDp $iS
		# Test MCMC convergence 
		python mcmc_testconvergence.py  $ipath  $opathPDp  $opathPDf $nRuns $iS
		# > "${opathPDf}/converstatsA${iS}.csv"



		# PLOTS MEDIAN AND 95 % percentile from posterior distribution
		#python plotMedianCI_curve.py 2 1 ./ ./ 1

done

"""Neha Khetan, Nov 2022-2023
   Post-MCMC runs to pool all the independent runs with "N" walkers in each
   Ensemble statistis and plotting for the ensemble population

   Compute Gelman Rubin Score 

   USAGE: 
        mcmc_testconvergence.py Ipath Opathplots Opathfiles NumIterations samNum """
      

import emcee
import corner
import sys
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from emcee.autocorr import integrated_time
plt.rcParams.update({'font.size': 18 })
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rc('legend', fontsize = 8 )
plt.tight_layout()


ipath	= sys.argv[1] 
opath   = sys.argv[2] 
opathFile   = sys.argv[3] 
numRuns = int( sys.argv[4] )
nSam    =int(sys.argv[5])
AllWalkers = None
AllLogProb = None
ndim = 8

# Gelman-Rubin statistic
def gelman_rubin(samples):
    nchains, nsamples  = samples.shape
    means = samples.mean(axis=1)
    variances = np.var(samples, axis=1, ddof=1)
    mean_of_means = means.mean(axis=0)
    B = nsamples * np.var(means, axis=0, ddof=1)
    W = variances.mean(axis=0)
    var_hat = ((nsamples - 1) / nsamples) * W + (1 / nsamples) * B
    R = np.sqrt(var_hat / W)
    return R

AllWalkers = None
AllLogProb = None
FinalVal     = []
sortedThetas = []
sortedScores = []
SuccessfulEvents = 0
numsteps  = 17000#17500
burnin    =  3000#2500 + 500 


SampSize = []
parmCol = 0
#print( "No. of steps:", numsteps )
for nfiles in range(  1 , numRuns  ):		
	filename     = ipath + "back_" + str( nfiles ) + ".hdf5"
	sampler      = emcee.backends.HDFBackend( filename  )

	# FILE EXISTS?
	try:
		tmp_samples     = sampler.get_chain(  discard=burnin , thin = 1 ,  flat=False )
		tmp_samples     = tmp_samples[0:numsteps,0:1 ,parmCol]
		#print( nfiles , tmp_samples.shape[0])
		SampSize.append( tmp_samples.shape[0])
	except:
		#print("Skipping:", nfiles )
		continue

	if tmp_samples.shape[0] < numsteps :
		continue

	# If desired steps!
	if SuccessfulEvents==0:
		AllWalkers     		   = sampler.get_chain( discard=burnin , thin= 1, flat =False )
		AllWalkers             = AllWalkers[0:numsteps,0:1,parmCol]
	
	else:
				
		AllWalkers             = np.append( AllWalkers , tmp_samples , axis = 1  )
	
	SuccessfulEvents       = SuccessfulEvents + 1
	if SuccessfulEvents > 100:
		#break
		pass

AllWalkers = np.transpose( AllWalkers )
GRstat     = np.array( [ 1 , round( gelman_rubin( AllWalkers ),2) ])
'''
dfGR			= pd.DataFrame(  {'GR':    GRstat } )
dfGR.to_csv( opathFile +  "GRstat_A" + str(nSam ) + ".csv"  )
'''
np.savetxt(opathFile +  "GRstat_A" + str(nSam ) + ".csv" , GRstat )
plt.hist( SampSize )
plt.savefig( opath + 'Dist_A' + str(nfiles) + '.png')
plt.close()

'''
Last Modified: 19 Nv. 2022
Last Nodified: 6 Nov. 2022

Sep - Oct, 2022

1. Pool different independent runs 
2. Pick the parameters from "top-best i.e. TopN " samples

Input: filename.py SampleNo NumRuns


USAGE: file.py SamNo NoRuns Ipath Opath NumIterations TopN 

'''



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import solve_ivp
import emcee
from scipy.optimize import minimize

from emcee.autocorr import integrated_time

import corner
import sys
import seaborn as sns
import pandas as pd 
plt.rcParams.update({'font.size': 18 })
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rc('legend', fontsize = 8 )
plt.tight_layout()

# Automated windowing procedure following Sokal (1989)
def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


# Following the suggestion from Goodman & Weare (2010)
def autocorr_gw2010(y, c=5.0):
    f = autocorr_func_1d(np.mean(y, axis=0))
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]

def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf


def autocorr_ml(y, thin=1, c=5.0):
    # Compute the initial estimate of tau using the standard method
    init = autocorr_new(y, c=c)
    z = y[:, ::thin]
    N = z.shape[1]

    # Build the GP model
    tau = max(1.0, init / thin)
    kernel = terms.RealTerm(
        np.log(0.9 * np.var(z)), -np.log(tau), bounds=[(-5.0, 5.0), (-np.log(N), 0.0)]
    )
    kernel += terms.RealTerm(
        np.log(0.1 * np.var(z)),
        -np.log(0.5 * tau),
        bounds=[(-5.0, 5.0), (-np.log(N), 0.0)],
    )
    gp = celerite.GP(kernel, mean=np.mean(z))
    gp.compute(np.arange(z.shape[1]))

    # Define the objective
    def nll(p):
        # Update the GP model
        gp.set_parameter_vector(p)

        # Loop over the chains and compute likelihoods
        v, g = zip(*(gp.grad_log_likelihood(z0, quiet=True) for z0 in z))

        # Combine the datasets
        return -np.sum(v), -np.sum(g, axis=0)

    # Optimize the model
    p0 = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    soln = minimize(nll, p0, jac=True, bounds=bounds)
    gp.set_parameter_vector(soln.x)

    # Compute the maximum likelihood tau
    a, c = kernel.coefficients[:2]
    tau = thin * 2 * np.sum(a / c) / np.sum(a)
    return tau


'''
======================================================================
...
...
...
...
======================================================================
'''
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
numsteps  = 1500#17500
burnin    =  500#2500


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


#-----------
#print( AllWalkers.shape)
AllWalkers = np.transpose( AllWalkers )
#print( AllWalkers.shape)
#print( "GR-statistic ", round( gelman_rubin( AllWalkers ),2) )
GRstat = np.array( [ 1 , round( gelman_rubin( AllWalkers ),2) ])
'''
dfGR			= pd.DataFrame(  {'GR':    GRstat } )
dfGR.to_csv( opathFile +  "GRstat_A" + str(nSam ) + ".csv"  )
'''
np.savetxt(opathFile +  "GRstat_A" + str(nSam ) + ".csv" , GRstat )


plt.hist( SampSize )
plt.savefig( opath + 'Dist_A' + str(nfiles) + '.png')
plt.close()
#print("Allwakker shape", AllWalkers.shape )

chain = AllWalkers
N = np.exp(np.linspace(np.log(10), np.log(chain.shape[1]), 3 )).astype(int)
gw2010 = np.empty(len(N))
new = np.empty(len(N))



for i, n in enumerate(N):
    gw2010[i] = autocorr_gw2010(chain[:, :n])
    new[i]    = autocorr_new(chain[:, :n])

# Plot the comparisons
fig = plt.figure()
plt.loglog(N, gw2010, "o-", label="G&W 2010")
plt.loglog(N, new, "o-", label="new")
ylim = plt.gca().get_ylim()
plt.plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
plt.ylim(ylim)
plt.xlabel("number of samples, $N$")
plt.ylabel(r"$\tau$ estimates")
plt.legend(fontsize=14);
plt.savefig( opath + 'hist_autocorr_A' + str(nfiles) + '.png')
plt.close()


#print( "Total Files:", SuccessfulEvents )
#print( AllWalkers.shape )
#Auto-correlation
#print(  "Integrated-time ", integrated_time(AllWalkers))
#print(   np.mean( integrated_time(AllWalkers)))
plt.hist( integrated_time(AllWalkers))
plt.savefig( opath + 'Autocorrelation.png')
plt.close()





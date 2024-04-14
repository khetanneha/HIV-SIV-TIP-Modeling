"""Neha Khetan, Nov 2022- 2023
   Post-MCMC runs to pool all the independent runs with "N" walkers in each
   Ensemble statistis and plotting for the ensemble population

   Test Homogenity across chains, calculate KS-scores on different chains"""
 

import sys
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

ipath 		  = sys.argv[1]
opath 		  = sys.argv[2]
nR    		  =	int( sys.argv[3] )
nSamplePlot   = 2000 #sys.argv[4]


for nfiles in range(  1, nR+1  ):	

	filename            = ipath + "AllScore_A" + str( nfiles ) + ".csv"
	#filename            = ipath + "AllScore_A5"  + ".csv"
	df 		 	        = pd.read_csv(  filename   )	
	df2   			    = df[df.columns[1:]]

	randDat = df2.sample( frac =0.50 , random_state = 100 )
	rc , nc = randDat.shape
	indx  = int( rc*0.5 )

	randA =  randDat.iloc[:indx, :]
	randB =  randDat.iloc[ indx:, :]
	#print( randA.shape , randB.shape )

	randA = randA.to_numpy().flatten()
	randB = randB.to_numpy().flatten()
	
	KSres = stats.ks_2samp( randA , randB )
	#print( KSres  )
	newdf = pd.DataFrame(  {'Sample A': randA , 'Sample B': randB })
	#ax= sns.histplot( newdf ,  bins=50 , stat="density", fill=0  ,element = "step" )
	ax= sns.histplot( newdf.iloc[1:nSamplePlot,: ] ,  bins=50 , stat="density", fill=0  ,element = "step" )
	#plt.xlim([ -50 , 0])
	plt.xlabel('Score' ) 
	plt.ylabel('Density' ) 
	plt.tight_layout()	
	ax.text(  -30 , 0.08 , np.round(  KSres.pvalue , 3 ) )
	plt.savefig( opath +'ScoreKS_' + str( nfiles )  + '.pdf' , dpi = 300 )
	plt.close() 
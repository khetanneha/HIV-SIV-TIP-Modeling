'''
Neha Khetan, Oct 2022 -  Dec,2023
Model parameterization of Viral kinetics from
rhesus macaques: SHIV infected [ Control Group ]
and SHIV infected + TIP-treated [ Treated Group ]

This is the main file that needs to be executed
Usage: 
	python run_mcmc_infanttips.py AnimalNo NumIndependentRuns OutputFolder_plots OutputFolder_files
'''



import emcee
import corner
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy import integrate, optimize
plt.rcParams.update({'font.size': 18 })
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rc('legend', fontsize = 8)
plt.tight_layout()




import extractExptData as ExpData
import compute_R0 as calR0
from GlobalparamInfantPVL import SysParam
PARAM = SysParam()
global ViralT0
from tipmodel1_InfantPVL_A1M1 import tipmodel1
'''
TIP 1:  no super-infection
TIP 2:  clearance rates
TIP 3:  superinfection m = 3
TIP 4 : super-infection m =2
TIP 5 : heterozygous | Assuming they do not infect new cells
TIP 6 : Cell division and growth
TIP 7:  immature virions that do not infect: with and w/o feedback on infectivity 
'''


def main( p0 ,nwalkers, niter,  ndim,  lnprob,  data , numRun ):
	"""Main mcmc module 
	"""

	OutFname = OPATH_files + 'back_' + str( numRun ) + '.hdf5'
	backEnd = emcee.backends.HDFBackend(  OutFname , name='mcmc', read_only=False, dtype=None, compression=None, compression_opts=None)
	backEnd.reset( nwalkers, ndim )

	sampler = emcee.EnsembleSampler( nwalkers, ndim, lnprob  , backend = backEnd , args=data  )
	print("Running burn-in...")
	p0, _, _ = sampler.run_mcmc( p0, numBurn , progress=True)
	print("Running production...")
	pos, prob, state = sampler.run_mcmc( p0 , niter, progress=True )	
	print("Mean acceptance fraction: {0:.3f}".format( np.mean(sampler.acceptance_fraction)))
	return sampler, pos, prob, state 


def plot_thetaHist( sam , numRun , varNo ):
	"""plots parameter explored in successive steps
	"""
	ss = sam.get_chain( flat= True )
	f10 , ax10 = plt.subplots(  varNo , 1 , figsize=( 12 , 10 ) )
	for fh in range( 0 , varNo ):
		#ax10[fh].hist( ss[:,fh] , color ="k" , histtype = "step")
		ax10[fh].plot( sam.chain[:,:,fh].T , color ="gray" , alpha = 0.5 )
		ax10[fh].plot( np.mean(sam.chain[:,:,fh].T , axis = 1 ) , color ="k" , alpha = 0.3 )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'ThetaTrajectories_' + str( numRun ) + '.pdf'  , dpi= 300  )
	plt.close()
	

def execute_mmc( ModelObj , AllDat , TIPDat , varNo , numRun , datatype ):
	"""Prepares and execute  MCMC: 
		1. Prepares experiment data 
		2. Initial conditions based on #1
		3. Initialize parameter-guess sets for walkers
		4. Invoke actual-mcmc-run"""

	nR , nC =  AllDat.shape 
	ylab    = AllDat.columns.values.tolist()[1:]
	print( "Sample #:", samNo , "Run #:" , numRun )	
	Fit_vals = []
	tmp 	 = AllDat.iloc[  : , [ 0, samNo ] ].dropna()  
	backGrnd = AllDat.iloc[  0 , [ 0, samNo ] ] 
	xx  	 = np.array( tmp.iloc[ StartRow:,0 ] )*cFWk2Day
	yy  	 = np.array( tmp.iloc[ StartRow:,1 ] )
				
	if datatype =='control':
		indx = yy[:] < LOD_siv
		if (indx.any() ):
			yy[ indx ] = LOD_siv
		yy[0] = 2*ViralT0 #LOD_siv # vRNA

	else:
		#print( backGrnd[0] , backGrnd[1] )
		# For. start norm
		if ( backGrnd[1] > LOD_siv ):
			yy[:] = yy[:] - backGrnd[1] 
		
		indx  = yy[:] < LOD_siv
		if (indx.any() ):
			yy[ indx ] = LOD_siv
		yy[0] = 2*ViralT0 #LOD_siv # vRNA	
	eXY = np.column_stack( (xx , yy ))  
	data = ( xx , yy , yerr )  


	eXY_tip = []
	if If_TIPDat:
		# tip data
		tmptip      = TIPDat.iloc[ : , [ 0, samNo ] ].dropna()  
		tipBackGrnd = TIPDat.iloc[ 0 , [ 0, samNo ] ]
		xtip  		= np.array( tmptip.iloc[ StartRow:,0 ] )*cFWk2Day
		ytip  		= np.array( tmptip.iloc[ StartRow:,1 ] )
			
		if ( tipBackGrnd[1] > 1.1 ):
			#ytip[:] = ytip[:] - tipBackGrnd[1]
			indx2   = ytip[:] < LOD_tip
			if (indx2.any() ):
				ytip[ indx2 ] = LOD_tip
		eXY_tip = np.column_stack( (xtip , ytip ))

	# ========================================
	initial_theta   = np.array( ModelObj.params )
	ndim            = len( initial_theta )
	iv0             = ModelObj.get_AllIV( )
   
	intgdTIP        = ExpData.get_dataInfantIntegratedTIP( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
	intgdTIPSam     = intgdTIP.iloc[  0 , samNo  ]	
	iv              = ModelObj.set_AllIV(  V0=ViralT0 , VT0=0 , IT10=intgdTIPSam ,  TotalT0=PARAM.TotalT0 ) # 

	# ==================================================
	p0              = np.empty(( nwalkers, ndim ))
	p0[:,0] 		= initial_theta[0] + 1e-3*np.random.randn( nwalkers )
	p0[:,1] 		= initial_theta[1] + 1e-3*np.random.randn( nwalkers )
	p0[:,2] 		= initial_theta[2] + 1e-5*np.random.randn( nwalkers )
	p0[:,3] 		= initial_theta[3] + 1e-3*np.random.randn( nwalkers )
	p0[:,4] 		= initial_theta[4] + 1e-5*np.random.randn( nwalkers )
	p0[:,5] 		= initial_theta[5] + 1e-3*np.random.randn( nwalkers )
	p0[:,6] 		= initial_theta[6] + 1e-3*np.random.randn( nwalkers )
	p0[:,7] 		= initial_theta[7] + 1e-3*np.random.randn( nwalkers )
	new_sampler, newpos, newprob, newstate = main( p0 , nwalkers , niter , ndim , ModelObj.lnprob , data , numRun )		
	new_samples    = new_sampler.flatchain
	sampler_lnprob = new_sampler.flatlnprobability
	new_theta_max  = new_samples[np.argmax(sampler_lnprob)]
	Fit_vals.append( new_theta_max )
	whoamI = ylab[ samNo -1 ]
	print( whoamI )
	#====
	'''
	# Calculate Reproductve ratio
	# estimate R0 from exponetial growth
	maxVLexp = np.max( yy[0:5] )
	maxInd   = yy[0:5].argmax()
	# Fitting the exponential growth model to the data
	params, covariance = curve_fit( exponential_growth, xx[0:maxInd] , yy[0:maxInd] )
	V0_fit, r_fit = params
	estGrowthR0 = 1 + ( r_fit/initial_theta[4] )
	print( "params-fit" , xx[0:maxInd] , yy[0:maxInd] , params , "R0:" ,r_fit  )	
	tmpval = np.append( new_theta_max , r_fit  )
	Fit_vals.append( tmpval )
	sys.exit()
	'''

	# plot the parameter trace
	plot_thetaHist( new_sampler  , numRun  , varNo )        
	plot_traj( ModelObj ,   eXY , whoamI , new_sampler  , numRun , eXY_tip , datatype )
	plt.plot( new_sampler.get_log_prob( flat=False )	 )
	plt.savefig( OPATH_plots + 'LogProbScore_' + str(numRun) + '.pdf'  , dpi = 300 )

	df2 = pd.DataFrame( Fit_vals , columns= ModelObj.labels )
	df2.to_csv( OPATH_files + 'Treated_Fit_' + str( numRun ) + '.csv' )
	#print( df2.mean())




def plot_traj(  ModelObj2 ,  expXY , whoamI  , new_sampler  , numRun  , eXYtip , datatype ):	

	"""All the plotting functions invoked: Evolution of each population in ODE model, Corner plots for posterior distribution, parameter traces, histograms of
	   of posterior distributions"""
	new_samples    		=  new_sampler.flatchain
	sampler_lnprob 		=  new_sampler.flatlnprobability
	new_theta_max  		=  new_samples[np.argmax(sampler_lnprob)]
	iv    		   		=  ModelObj2.get_AllIV()
	fineC 		   		=  np.linspace( np.min( expXY[:,0] ) , np.max( expXY[:,0] ) , tot_timPt  )
	new_best_fit_model  =  ModelObj2.model( new_theta_max, iv, fineC )

	# ============================================================
	# Plots T, and I  from the best fit
	# ============================================================
	f2 , ax2 = plt.subplots()#( figsize=( 6,4 ))
	CF = PARAM.CFml2ul
	ax22 = ax2.twinx()

	if TIPModel == 1:
		HealthyCells   = ( new_best_fit_model.y[0,:] + new_best_fit_model.y[3,:] )*CF
		UnHealthyCells = ( new_best_fit_model.y[1,:] + new_best_fit_model.y[4,:] )*CF  
		TIPCells       = ( new_best_fit_model.y[3,:]   )*CF

	elif TIPModel == 5:
		HealthyCells   = ( new_best_fit_model.y[0,:] + new_best_fit_model.y[3,:] )*CF
		UnHealthyCells = ( new_best_fit_model.y[1,:] + new_best_fit_model.y[4,:] )*CF  
		TIPCells       = ( new_best_fit_model.y[3,:]   )*CF

	elif TIPModel == 3:
		HealthyCells   = ( new_best_fit_model.y[0,:] + new_best_fit_model.y[3,:] + new_best_fit_model.y[4,:] + new_best_fit_model.y[5,:] )*CF
		UnHealthyCells = ( new_best_fit_model.y[1,:] + new_best_fit_model.y[6,:] + new_best_fit_model.y[7,:] + new_best_fit_model.y[8,:] )*CF  
		TIPCells       = ( new_best_fit_model.y[3,:] + new_best_fit_model.y[4,:] + new_best_fit_model.y[5,:]  )*CF

	elif TIPModel == 4:
		HealthyCells   = ( new_best_fit_model.y[0,:] + new_best_fit_model.y[3,:] + new_best_fit_model.y[4,:]  )*CF
		UnHealthyCells = ( new_best_fit_model.y[1,:] + new_best_fit_model.y[5,:] + new_best_fit_model.y[6,:]  )*CF  
		TIPCells       = ( new_best_fit_model.y[3,:] + new_best_fit_model.y[4,:]   )*CF

	ax2.plot( new_best_fit_model.t , HealthyCells   , '-', c='k'  ,   markersize = 0.5  , label='T + Tt')
	ax2.plot( new_best_fit_model.t , UnHealthyCells , '-',  c ='gray'    ,   markersize = 0.5  , label='I + Id')  
	ax22.plot( new_best_fit_model.t , (TIPCells*100)/( HealthyCells + UnHealthyCells ) , '-',  c ='cornflowerblue'    ,    markersize = 0.5  , label='I + Id')   
	if isCD4data:
		cdData = ExpData.get_cd4dataInfant( datatype , ExpData.Animal2ExcludeInfant( datatype ))
		ax2.scatter(  np.array( cdData.iloc[:,0] )*cFWk2Day , np.array( cdData.iloc[:,samNo] ) ,   c='gray' , linewidth = 1.2 )  
	ax2.set_xlabel('Days' ) 
	ax2.set_ylabel('CD4+ cells/uL' )
	ax2.legend()
	plt.tight_layout()
	plt.savefig( OPATH_plots +  'BestFit_CD4_TIP_' + str( numRun ) + '.pdf'  , dpi = 300 )
	plt.close()   




	# ============================================================
	# VIRAL populations : V and VT
	# ============================================================
	f3   =  plt.figure( figsize=( 6 ,4 ) ) 
	plt.yscale('log')
	plt.plot( expXY[:,0] , expXY[:,1] , 'ko-.' ,  label= whoamI ,alpha = 0.8 , linewidth = 1.2 , markersize = 4 )  
	if TIPModel == 1:
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[5,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  

	elif TIPModel == 5:
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[5,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[6,:],'c-' , alpha = 0.3 , markersize = 1 , label='HH' )  

	elif TIPModel == 3:
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[9,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  
		
	elif TIPModel == 4:
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( new_best_fit_model.t , v2RNA*new_best_fit_model.y[7,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  

	plt.axhline( y= ModelObj2.LOD_siv, color='k', linestyle='-.')
	TimeExpt = expXY[:,0]
	newFitC  = ModelObj2.model( new_theta_max, iv, TimeExpt )
	if TIPModel == 1:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[5,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
	elif TIPModel == 3:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[9,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    

	elif TIPModel == 4:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[7,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    

	elif TIPModel == 5:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[5,:] + newFitC.y[6,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT+Vhh' )    



	if If_TIPDat:
		plt.plot( eXYtip[:,0] , eXYtip[:,1] , 'o-.' , color='green' , alpha = 0.5 ,   label= 'TIP-GFP' , linewidth = 1.2  , markersize = 3 )  	
	plt.xlabel('Days' ) 
	plt.ylabel('Gag RNA (copies/ml) ') 
	plt.legend()
	plt.tight_layout()
	plt.savefig(  OPATH_plots + 'BestFitVLs_Treated_' + str( numRun )  + '.pdf'  , dpi = 300 )
	plt.close()
	


	f30   =  plt.figure( figsize=( 6 ,4 ) ) 
	plt.yscale('log')
	plt.plot( expXY[:,0] , expXY[:,1] , 'ko-.' ,  label= whoamI ,alpha = 0.8 , linewidth = 1.5 , markersize = 4 )  
	plt.axhline( y=ModelObj2.LOD_siv, color='k', linestyle='-.')
	TimeExpt = expXY[:,0]
	newFitC = ModelObj2.model( new_theta_max, iv, TimeExpt )
	if TIPModel == 1:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[5,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
	elif TIPModel == 3:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[9,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
	elif TIPModel == 5:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[5,:] + newFitC.y[6,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT + V$_{H}$' )    
	elif TIPModel == 4:
		plt.plot( newFitC.t , v2RNA*( newFitC.y[2,:]+ newFitC.y[7,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VH' )    
	plt.xlabel('Days' ) 
	plt.ylabel('Gag RNA (copies/ml) ') 
	plt.legend()
	plt.tight_layout()
	plt.savefig(  OPATH_plots + 'BestFitGAG_Treated_' + str( numRun )  + '.pdf'  , dpi = 300 )
	plt.close()


	# ============================================================
	# CORNER PLOTS:
	# ============================================================
	f10    = plt.figure( figsize =(6,4) )
	labels = ModelObj2.labels
	f10    = corner.corner( new_samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles= [0.25, 0.5, 0.75]) 
	plt.tight_layout()
	plt.savefig( OPATH_plots + 'cornerTreated_' + str( numRun ) + '.pdf'  , dpi = 300 )
	plt.close()


	f300   =  plt.figure( figsize=( 6 ,4 ) ) 
	plt.yscale('log')
	plt.plot( expXY[:,0] , expXY[:,1] , 'ko' ,  label= whoamI ,alpha = 0.6  , markersize = 2 )  
	plt.axhline( y=ModelObj2.LOD_siv, color='k', linestyle='-.')

	TimePredict = np.linspace( 0, 365*3 , int( (365*3) ) )
	newFitC2 = ModelObj2.model( new_theta_max, iv, TimePredict )
	if TIPModel == 1:
		#plt.plot( newFitC2.t , v2RNA*( newFitC2.y[2,:]+ newFitC2.y[5,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[5,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  
	if TIPModel == 5:
		#plt.plot( newFitC2.t , v2RNA*( newFitC2.y[2,:]+ newFitC2.y[5,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[5,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[6,:],'b-' , alpha = 0.3 , markersize = 1 , label='HH' )  


	elif TIPModel == 3:
		#plt.plot( newFitC2.t , v2RNA*( newFitC2.y[2,:]+ newFitC2.y[9,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[9,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  

	elif TIPModel == 4:
		#plt.plot( newFitC2.t , v2RNA*( newFitC2.y[2,:]+ newFitC2.y[9,:] ), 'o-', color= 'gray' , alpha = 0.8 , markersize = 4 , label='V+VT' )    
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[2,:],'r-' , alpha = 0.3 , markersize = 1 , label='SIV' )  
		plt.plot( newFitC2.t , v2RNA*newFitC2.y[7,:],'b-' , alpha = 0.3 , markersize = 1 , label='TIP' )  
	plt.xlabel('Days' ) 
	plt.ylabel('Gag RNA (copies/ml) ') 
	plt.legend()
	plt.tight_layout()
	plt.savefig(  OPATH_plots + 'PredictBestFitGAG_Treated_' + str( numRun )  + '.pdf'  , dpi = 300 )
	plt.close()






if __name__ == '__main__':
	"""Reads in all the parameters-associated with simulation-system from the GlobalparaminfantPVL.py. 
	   Extracts experimental data for PVL and TIP 
	   Initiates MCMC module 
	"""

	datatype            = PARAM.datatype
	varNo               = PARAM.varNo
	TIPModel  	   		= PARAM.TIPModel      
	numBurn     		= PARAM.numBurn
	niter       		= PARAM.niter   
	yerr        		= PARAM.MLsigma
	nwalkers    		= PARAM.nwalkers
	predicttime 		= PARAM.predicttime
	tot_timPt   		= PARAM.tot_timPt
	isCD4data           = PARAM.isCD4data
	If_TIPDat   		= PARAM.If_TIPDat
	cFWk2Day    		= PARAM.cFWk2Day
	StartRow       		= PARAM.StartRow
	LOD_tip        		= PARAM.LOD_tip
	LOD_siv		   		= PARAM.LOD_siv
	v2RNA 		   		= PARAM.v2RNA
	TotalT0		   		= PARAM.TotalT0
	ViralT0             = PARAM.ViralT0
	OPATH_plots 		= sys.argv[3] 
	OPATH_files 		= sys.argv[3]    
	samNo  				= int( sys.argv[1] )
	numRun 				= int( sys.argv[2] ) 
	np.random.seed( numRun )
	print("Performing MCMC fits for Infant PVL.........")
	# Reading data
	AllDat = ExpData.get_dataInfantPVL(    datatype , ExpData.Animal2ExcludeInfant( datatype ) )
	# TIP- data
	TIPDat = []
	if If_TIPDat:
		TIPDat = ExpData.get_TIPdataInfant(    datatype , ExpData.Animal2ExcludeInfant( datatype ) )	
	# Run mcmc
	ModelObj        = tipmodel1()
	execute_mmc( ModelObj  , AllDat  ,  TIPDat ,  varNo , numRun , datatype )
"""Neha Khetan, 2023
   Post-MCMC runs to pool all the independent runs with "N" walkers in each
   Ensemble statistis and plotting for the ensemble population

   USAGE: 
   	poolanalysisnRuns.py SamNo NoRuns Ipath Opath NumIterations TopN """


import numpy as np
import pandas as pd
import emcee
import corner
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import solve_ivp
plt.rcParams.update({'font.size': 18 })
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rc('legend', fontsize = 8 )
plt.tight_layout()

import extractExptData as ExpData
from GlobalparamInfantPVL import SysParam
PARAM = SysParam()
from tipmodel1_InfantPVL_A1M1 import tipmodel1
#from tipmodel11_InfantPVL import tipmodel1
import compute_R0 as calR0


def analysis( M0, AllDat , TIPDat ,  nR , numsteps ,  burnin  , TopN , datatype ):
	"""Read in-all the MCMC chains across all the runs, sorts based on scores can select "TopN" or choose all
	   solve ODE model
	   Select the parameter set that satify "realistic biological criterion" and Ro-TIP > 1 and 
	   and calculate statistics from population and use for plotting and o/p into files
	"""
	RTIPChoiceVal = 1  # Cutoff-value; select only parameter sets with RTIP > 1
	# EXPT DATA
	ylab    = AllDat.columns.values.tolist()[1:]
	tmp 	 = AllDat.iloc[  : , [ 0, samNo ] ].dropna()  
	backGrnd = AllDat.iloc[  0 , [ 0, samNo ] ] 
	xx  	 = np.array( tmp.iloc[ StartRow:,0 ] )*cFWk2Day
	yy  	 = np.array( tmp.iloc[ StartRow:,1 ] )	
	whoamI   = ylab[ samNo - 1]
	print(" Data-ID" , whoamI )
	# Account for LOD: VIRAL LOAD
	if datatype =='control':
		indx = yy[:] < LOD_siv
		if (indx.any() ):
			yy[ indx ] = LOD_siv
		yy[0] = 2*ViralV0
	else:
		#print( backGrnd[0] , backGrnd[1] )
		# For. start norm
		if ( backGrnd[1] > LOD_siv ):
			yy[:] = yy[:] - backGrnd[1] 
			indx  = yy[:] < LOD_siv
			if (indx.any() ):
				yy[ indx ] = LOD_siv

			yy[0] = 2*ViralV0
	eXY = np.column_stack( (xx , yy ))  
	data = ( xx , yy ) 

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
	# ==================================================================================================
	AllWalkers = None
	AllLogProb = None
	AllScores  = []
	FinalVal     = []
	sortedThetas = []
	sortedScores = []
	SuccessfulEvents = 0
	print( "No. of steps:", numsteps )

	for nfiles in range(  1 , nR +1   ):		
		filename     = ipath + "back_" + str( nfiles ) + ".hdf5"
		sampler      = emcee.backends.HDFBackend( filename  )	
		# FILE EXISTS
		try:
			tmp_samples     = sampler.get_chain(  discard=burnin , flat=False )
			tmp_samples     = tmp_samples[0:numsteps,:,:]

			tmp_logprob     = sampler.get_log_prob( discard=burnin, flat=False) 
			tmp_logprob     = tmp_logprob[0:numsteps,:]
			#print( nfiles , tmp_samples.shape[0])
		except:
			print("Skipping:", nfiles )
			continue

		if tmp_samples.shape[0] < numsteps:
			continue

		# If desired steps!
		if SuccessfulEvents==0:
			AllWalkers     		   = sampler.get_chain( discard=burnin , flat =False )
			AllLogProb     		   = sampler.get_log_prob(discard=burnin, flat=False )
			AllWalkers             = AllWalkers[0:numsteps,:,:]
			AllLogProb             = AllLogProb[0:numsteps,:]
			tmpVal 				   = select_topScores( AllWalkers , AllLogProb , TopN  )
	
			for tmpv in tmpVal:
				FinalVal.append( tmpv )		
			AllScores.append(  AllLogProb  )
			#print( AllLogProb )	
		else:
					
			AllWalkers             = np.append( AllWalkers , tmp_samples , axis = 1  )
			# Choose number of fits to pick
			tmpVal                 = select_topScores( tmp_samples , tmp_logprob , TopN  )		
			for tmpv in tmpVal:
				FinalVal.append( tmpv )

			tmp_score = sampler.get_log_prob(discard=burnin, flat=False )
			tmp_score             = tmp_score[0:numsteps,:]
			AllScores.append(  tmp_score  )
			#print( tmp_score )
			
		SuccessfulEvents       = SuccessfulEvents + 1
		if SuccessfulEvents > 500:
			pass#break
	'''
	===================================================================================================================
	SCORES- Log-likelihood Scores
	'''
	AllScores = np.array( AllScores )	
	newarr = np.empty( ( numsteps , 0 ) )
	for ii in range( 0 , AllScores.shape[0] ):
		tmp = AllScores[ ii , : , :]
		newarr = np.hstack( ( newarr, np.array(tmp ) ) )
	df_LogScore   =pd.DataFrame( newarr )
	df_LogScore.to_csv( OPATH_files + 'AllScore_A' + str( samNo ) + '.csv' , float_format='%.2f' )
	
	# Sorting by scores - max to min
	maxL 		   = sorted( FinalVal , key=lambda x: x[0] , reverse = True )
	FinalVal 	   = maxL
	fineTime       = np.linspace( np.min( eXY[:,0] ) , np.max( eXY[:,0] ) , tot_timPt  )
	#fineTime      = eXY[:,0] #np.linspace( np.min( eXY[:,0] ) , np.max( eXY[:,0] ) , tot_timPt  )
	
	#------------------------------------
	AllR0       = np.zeros( len(FinalVal)  )
	AllRT       = np.zeros( len(FinalVal)  )		
	AllCD4Cells = []
	AllVL       = []
	vTotGag     = []
	vTotTIP     = []
	AllTIP      = [] 
	PerInfected = np.zeros( len(FinalVal)  )
	perTipInfectedCells = np.zeros( len(FinalVal)  )
	EndTotalCD4 = np.zeros( len(FinalVal)  )
	OutData    	= [] 
	SelectData 	= []
	ALLT_TipCells =[]
	ALLT_VirCells =[]
	All_VirOnly   = []

	f1 = plt.figure( 1 ,  figsize=( 6,4 ))
	f2 = plt.figure( 2 ,  figsize=( 6,4 ))
	f3 = plt.figure( 3 ,  figsize=( 6,4 ))
	f4 = plt.figure( 4 ,  figsize=( 6,4 ))
	f5 = plt.figure( 5 ,  figsize=( 6,4 ))
	cmap  = plt.get_cmap("copper")
	ccmap = cmap( np.linspace( 0 , 1,  len(FinalVal) ) )
	CF 	  = 10**(-3) # from mL to uL
	cc 	  = 0
	
	v2RNA = V2RNA
	print("============================================================ \n") 
	print(" In poolAnalysis_V10 \n" )
	tmpScore = 0
	for ik in range( 0 , len( FinalVal ) ) :
		THETA           = FinalVal[ik][1]
		iv0             = M0.get_AllIV( )
		# ==================================================
		intgdTIP        = ExpData.get_dataInfantIntegratedTIP( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
		intgdTIPSam     = intgdTIP.iloc[  0 , samNo  ]		

		if PARAM.TIPModel == 11:
			iv                    = ModelObj.set_AllIV(  V0=ViralT0 , VIM0=0 ,   TotalT0=PARAM.TotalT0 ) # 
		else:
			iv                    = M0.set_AllIV(  V0=ViralV0 , VT0=0 , IT10=intgdTIPSam ,  TotalT0=PARAM.TotalT0 ) # 		


		new_mdl =  M0.model( THETA , iv  ,  fineTime )
		# Standard TIP-model		
		if PARAM.TIPModel == 1:
			HealthyCells   = ( new_mdl.y[0,:] + new_mdl.y[3,:] )*CF
			UnHealthyCells = ( new_mdl.y[1,:] + new_mdl.y[4,:] )*CF  
			TIPCells       =  new_mdl.y[3,:]*CF	
			vTotGag            =  V2RNA*( new_mdl.y[2,:] + new_mdl.y[5,:] )
			vTotTIP            =  V2RNA*new_mdl.y[5,:]
			vTotHIV            =  V2RNA*new_mdl.y[2,:]
		# super-infection with m =3
		elif PARAM.TIPModel == 3:
			HealthyCells   = ( new_mdl.y[0,:] + new_mdl.y[3,:] + new_mdl.y[4,:] + new_mdl.y[5,:] )*CF
			UnHealthyCells = ( new_mdl.y[1,:] + new_mdl.y[6,:] + new_mdl.y[7,:] + new_mdl.y[8,:] )*CF  
			TIPCells       =  ( new_mdl.y[3,:] + new_mdl.y[4,:] + new_mdl.y[5,:] ) *CF
			vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[9,:] )
			vTotTIP            =  v2RNA*new_mdl.y[9,:] 			
			vTotHIV            =  v2RNA*new_mdl.y[2,:]
		# Heterozygous model
		elif PARAM.TIPModel == 5:
			HealthyCells   = ( new_mdl.y[0,:] +  new_mdl.y[3,:]  )*CF
			UnHealthyCells = ( new_mdl.y[1,:] +  new_mdl.y[4,:] )*CF  
			TIPCells       =  ( new_mdl.y[3,:]  ) *CF
			vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[5,:] + new_mdl.y[6,:]  )
			vTotTIP            =  v2RNA*new_mdl.y[5,:] 			
			vTotHIV            =  v2RNA*new_mdl.y[2,:]
		# super-infection with m =2 
		elif PARAM.TIPModel == 4:
			HealthyCells   = ( new_mdl.y[0,:] + new_mdl.y[3,:] + new_mdl.y[4,:]  )*CF
			UnHealthyCells = ( new_mdl.y[1,:] + new_mdl.y[5,:] + new_mdl.y[6,:]  )*CF  
			TIPCells       = ( new_mdl.y[3,:] + new_mdl.y[4,:]   )*CF
			vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[7,:] )
			vTotTIP            =  v2RNA*new_mdl.y[7,:] 			
			vTotHIV            =  v2RNA*new_mdl.y[2,:]
		# Immature virions in Basic model ; can be due to AART model or "seeding-hypothesis" model - 4 equation
		elif ( PARAM.TIPModel == 10 or PARAM.TIPModel == 11 ):
			HealthyCells   = ( new_mdl.y[0,:] )*CF
			UnHealthyCells = ( new_mdl.y[1,:]  )*CF  
			TIPCells       = 0 
			vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[3,:] )
			vTotTIP            =  v2RNA*new_mdl.y[3,:] 			
			vTotHIV            =  v2RNA*new_mdl.y[2,:]

		TotalCD4       	    = ( HealthyCells + UnHealthyCells ) 
		percentInfected     = ( UnHealthyCells/TotalCD4 )*100
		perTipInfectedCells = ( TIPCells/TotalCD4 )*100
		PerInfected[ik]     =  percentInfected[ -1 ]
		EndTotalCD4[ik]     =  TotalCD4[ -1 ]



		if ( (np.max( TotalCD4 ) < CutTotalCD4 )  and ( percentInfected[ -1 ] < CutOffPerInf ) and ( calR0.estimate_RTIP( THETA )> RTIPChoiceVal ) ):
			"""Select the parameter set that satify "realistic biological criterion"""
			if FinalVal[ik][0] < tmpScore:  
				"""Choose unique parameter-set """
				AllCD4Cells.append( TotalCD4 )
				AllVL.append(  vTotGag )
				AllTIP.append( vTotTIP )		
				ALLT_TipCells.append( perTipInfectedCells )
				ALLT_VirCells.append( percentInfected )
				All_VirOnly.append( vTotHIV  )
				# ========================================================================================
				AllR0[ik]         =  calR0.estimate_RR( THETA )   
				AllRT[ik]         =  calR0.estimate_RTIP( THETA )   
				# ========================================================================================

				t1 = np.array( FinalVal[ik][1] ).tolist()
				OutData.append( [ FinalVal[ik][0] , t1[0] , t1[1] ,t1[2] , t1[3] , t1[4] , t1[5] ,  t1[6] , t1[7] ,    AllR0[ik] , PerInfected[ik] , EndTotalCD4[ik] ,  AllRT[ik] ])

				plt.figure( f1 )
				plt.plot( new_mdl.t , vTotGag , 'o-',   c=ccmap[cc] ,   linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )  

				plt.figure( f3 )		
				plt.plot( new_mdl.t ,  TotalCD4 , 'o-',   c=ccmap[cc] ,   linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )

				plt.figure( f4 )	
				plt.plot( new_mdl.t ,  percentInfected  , 'o-',  c=ccmap[cc] ,  linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )

				plt.figure( f5 )
				plt.plot( new_mdl.t ,  vTotTIP , 'o-',  c=ccmap[cc] ,  linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )  

				SelectData.append( [ FinalVal[ik][0] , t1[0] , t1[1] ,t1[2] , t1[3] , t1[4] , t1[5] , t1[6] , t1[7] ,   AllR0[ik] , PerInfected[ik] , EndTotalCD4[ik] ,  AllRT[ik] ])
				tmpScore = np.Inf#FinalVal[ik][0]
			cc = cc + 1

	## ========================================================================================================
	AnimalMarker = [ 'o', '^' , 'v' , 's' , '*' ,'d']
	# Plotting Confidence Interval
	datMean , sdV = calMeanCILog( new_mdl.t ,  AllVL    )
	lowBV  = datMean - sdV
	upBV   = datMean + sdV
	lowBV[ lowBV < LOD_siv ] = LOD_siv

	f1000         = plt.figure( 1000 ,  figsize=( 6,4 ))		
	plt.fill_between( new_mdl.t, lowBV   , upBV , color = 'cornflowerblue' , alpha = 0.3  , label = whoamI )
	plt.plot( new_mdl.t ,  datMean   ,   color = 'cornflowerblue' , linewidth=2 , label = 'fit')
	plt.plot( eXY[:,0] , ( eXY[:,1] ) , marker = AnimalMarker[samNo] ,  color='cornflowerblue',   linewidth = 0 ,  markersize=8 ,markeredgecolor ='k' ,)		
	#peakVL = round( np.max( ( datMean)  ), 3 )
	#spVL   = round( ( datMean[-1] ), 3 )
	#print( "PVL", peakVL )
	#print( "SPVL", spVL )
	#print(  "ratio", round( peakVL/spVL , 3 ) )
	dfGag  = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': datMean , 'Column3': sdV })
	dfGag.to_csv( OPATH_files + 'Total_GagA' + str( samNo ) + '.csv' )
		
	datMean , sdV = calMeanCILog( new_mdl.t ,  All_VirOnly   )
	lowBV  = datMean - sdV
	upBV   = datMean + sdV
	lowBV[ lowBV < LOD_siv ] = LOD_siv
	f1000         = plt.figure( 1000 ,  figsize=( 6,4 ))		
	plt.fill_between( new_mdl.t, lowBV   , upBV , color = 'salmon' , alpha = 0.3  , label = whoamI )
	plt.plot( new_mdl.t ,  datMean   ,   color = 'salmon' , linewidth=2 , label = 'fit')
	dfVirOnly = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': datMean , 'Column3': sdV })
	dfVirOnly.to_csv( OPATH_files + 'OnlyVIR_GagA' + str( samNo ) + '.csv' )
	

	# Plotting Confidence Interval
	TIPdatMean ,  sdV= calMeanCItipLog( new_mdl.t , AllTIP    )
	lowBV  = TIPdatMean - sdV
	upBV   = TIPdatMean + sdV
	lowBV[ lowBV < LOD_tip ] = LOD_tip
	plt.fill_between( new_mdl.t, lowBV , upBV  , color = 'darkseagreen' , alpha = 0.3   )
	plt.plot( new_mdl.t , TIPdatMean  , color = 'darkseagreen' , linewidth=2 )

	if If_TIPDat:
		plt.plot( eXY_tip[:,0] , ( eXY_tip[:,1] ) , marker = AnimalMarker[samNo]  , color=  'darkseagreen',   linewidth=0 , markersize=8 , alpha = 0.7 , markeredgecolor = 'k')
	dfTGag  = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': TIPdatMean , 'Column3': sdV })
	dfTGag.to_csv( OPATH_files + 'Only_TIPA' + str( samNo ) + '.csv' )	
	plt.ylim( [ 10**0 , 10**9] )
	plt.xlim( [  -2 , 220 ] )
	plt.axhline( y=  5  , color = 'cornflowerblue' ,  alpha = 0.5 ,  linestyle='-.')
	plt.axhline( y= 500 , color = 'darkseagreen' ,  alpha = 0.5 , linestyle='-.')
	plt.yscale('log' )
	plt.xlabel('Days' ) 
	plt.ylabel('Log$_{10}$ RNA copies/ml' ) 
	plt.tight_layout()	
	plt.legend()
	plt.savefig( OPATH_plots +'sdVL_both' + str( samNo )  + '.pdf' , dpi = 300 )
	plt.close() 
	
	# Plotting Confidence Interval
	datMean ,  sdV = calMeanCI( new_mdl.t ,  AllCD4Cells  )
	f1000 = plt.figure( 1000 ,  figsize=( 6,4 ))
	plt.fill_between( new_mdl.t, (datMean - sdV) , ( datMean + sdV ) , color = 'b' , alpha = 0.1   , label = whoamI )
	plt.plot( new_mdl.t , datMean , color = 'k' , alpha =0.1, linewidth=2 )
	df_allCD4  = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': datMean , 'Column3': sdV })
	df_allCD4.to_csv( OPATH_files + 'Only_CD4_A' + str( samNo ) + '.csv' )	

	if isCD4data:
		cdData = ExpData.get_cd4dataInfant( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
		cdxy   = ( cdData.iloc[:,[0 , samNo ]].dropna()  ) 
		plt.plot(  np.array( cdxy.iloc[:,0] )*cFWk2Day , np.array( cdxy.iloc[:,1] ) , 'ko'  ,alpha = 0.6 )  
	plt.xlabel('Days' ) 
	plt.ylabel('CD4+ cells/uL' ) 
	plt.ylim( [ 0 , 4000 ] )
	plt.xlim( [  -2 , 220 ] )
	plt.tight_layout()	
	#plt.legend()
	plt.savefig( OPATH_plots +'sdCD4_S' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 
	

	# Plotting Confidence Interval
	datMean , sdV = calMeanCI( new_mdl.t ,  ALLT_TipCells  )
	f2003 = plt.figure( 2003 ,  figsize=( 6,4 ))
	plt.plot( new_mdl.t   , datMean , color = 'coral' , alpha =0.1, linewidth=2 )
	plt.fill_between( new_mdl.t, (datMean - sdV) , ( datMean + sdV ) , color = 'coral' , alpha = 0.1   , label = whoamI )
	#plt.ylim( 0 , 5 )	
	#plt.xlim( [  0 , 220 ] )
	#plt.xlabel('Days' ) 
	#plt.ylabel('% TIP integrated cells' ) 
	#plt.tight_layout()	
	#plt.legend()
	df_allTIPCD4  = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': datMean , 'Column3': sdV })
	df_allTIPCD4.to_csv( OPATH_files + 'Only_CD4withTIP_A' + str( samNo ) + '.csv' )	
	#plt.savefig( OPATH_plots +'perIntegratedTIP_' + str( samNo ) + '.pdf' , dpi = 300 )
	#plt.close() 

	# Plotting Confidence Interval
	datMean , sdV = calMeanCI( new_mdl.t ,  ALLT_VirCells  )
	#f2031 = plt.figure( 2031 ,  figsize=( 6,4 ))
	plt.plot( new_mdl.t   , datMean , color = 'darkseagreen' , alpha =0.1, linewidth=2 )
	plt.fill_between( new_mdl.t, (datMean - sdV) , ( datMean + sdV ) , color = 'darkseagreen' , alpha = 0.1   , label = whoamI )
	plt.ylim( 0 , 100 )	
	plt.xlim( [  -2 , 220 ] )
	plt.xlabel('Days' ) 
	plt.ylabel('% Vir infected cells' ) 
	plt.tight_layout()	
	#plt.legend()
	df_allTIPCD4  = pd.DataFrame(  {'Column1': new_mdl.t , 'Column2': datMean , 'Column3': sdV })
	df_allTIPCD4.to_csv( OPATH_files + 'Only_CD4withVir_A' + str( samNo ) + '.csv' )	
	plt.savefig( OPATH_plots +'perIntegrated_TIP_nVir_' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 
	
	# =======================
	plt.figure(1)
	plt.plot( eXY[:,0] , eXY[:,1] , 'bo' , label=ylab[ samNo - 1] )
	plt.yscale('log')
	plt.xlabel('Days' ) 
	plt.ylabel('Gag-RNA (copies/ml)' ) 
	plt.tight_layout()	
	plt.legend()
	plt.savefig( OPATH_plots +'Fits_VL_S' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 


	plt.figure(3)
	if isCD4data:
		cdData = ExpData.get_cd4dataInfant( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
		cdxy   = ( cdData.iloc[:,[0 , samNo ]].dropna()  ) 
		plt.plot(  np.array( cdxy.iloc[:,0] )*cFWk2Day , np.array( cdxy.iloc[:,1] ) , 'bo'   , label=ylab[ samNo - 1]  )  
	plt.xlabel('Days' ) #,fontsize = 16)
	plt.ylabel('CD4+ cells/uL' )
	plt.legend()
	plt.tight_layout()
	plt.savefig( OPATH_plots +'Fits_TotalCD4_S' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 

	plt.figure(4)
	plt.xlabel('Days' ) #,fontsize = 16)
	plt.ylabel('% Infected cells' )
	#plt.legend()
	plt.tight_layout()
	plt.savefig( OPATH_plots +'PercentInfected_S' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 

	plt.figure(5)
	if If_TIPDat:
		plt.plot( eXY_tip[:,0] , eXY_tip[:,1] , 'bo-.')
	plt.yscale('log')
	plt.xlabel('Days' ) 
	plt.ylabel('Gag-RNA (copies/ml)' ) 
	plt.tight_layout()	
	plt.legend()
	plt.savefig( OPATH_plots +'Fits_TIP_' + str( samNo ) + '.pdf' , dpi = 300 )
	plt.close() 
	print( "Total files:......." , SuccessfulEvents )
	print( "Saving the values into file....")
	df2 			= pd.DataFrame(  OutData   , columns=[ 'Score' , 'Lam' , 'bet' , 'd' , 'k' , 'a' , 'c'  , 'D' , 'P' ,'R0' ,'PerInf' , 'TotCD4'  ,'RT'] )
	df2.to_csv( OPATH_files + 'TIP_FitN_S' + str( samNo ) + '.csv' )
	dfSelect        =	pd.DataFrame(  SelectData   , columns=[ 'Score' , 'Lam' , 'bet' , 'd' , 'k' , 'a' , 'c' , 'D' , 'P' , 'R0' ,'PerInf' ,  'TotCD4' ,'RT'] )
	dfSelect.to_csv( OPATH_files + 'FilterTIP_FitN_S' + str( samNo ) + '.csv' )
	print( "Plotting the stats....")
	plotHistFitStats( dfSelect  , varNo , TopN  , fineTime , eXY , datatype ,  'Filter' )
	plot_thetaHist( AllWalkers  , varNo , TopN , 'All' ) 
	#plot_MedCIpd( dfSelect )	
	return

def select_topScores( allwalkers , alllogprob , topVals ):

		flatAllWalkerChains = allwalkers.reshape(-1, allwalkers.shape[-1] )
		flatAllLogProb      = alllogprob.flatten()		
		tmpsortedThetas 	= flatAllWalkerChains[ np.argsort( flatAllLogProb )[::-1]] 
		tmpsortedScores     = np.sort( flatAllLogProb )[::-1]
		tmpVal = []
		for i in range( 0 , topVals ):
			tmpVal.append( [ tmpsortedScores[i] , tmpsortedThetas[i] ])
		return tmpVal


def plotHistFitStats( sam , varNo  , topVn , fineTime ,  eXY , datatype , Sel_type ):
	"""Plot histogram of all fit-statistics"""
	F_lam  = 	sam.loc[:,"Lam"]
	F_beta =    sam.loc[:,"bet"]
	F_k    =    sam.loc[:,"k"]
	F_d    =    sam.loc[:,"d"]
	F_a    =    sam.loc[:,"a"]
	F_c    =    sam.loc[:,"c"]
	F_rho  =    sam.loc[:,"P"]
	F_psi  =    sam.loc[:,"D"]
	F_R0   =    sam.loc[:,"R0"]
	F_RT   =    sam.loc[:,"RT"]
	f100 , ax100 = plt.subplots(  varNo+2 , 1 , figsize=( 12 , 10 ) )	
	ax100[0].hist( F_lam  , bins =100,  density = True )
	ax100[1].hist( F_beta , bins =100,  density = True )
	ax100[2].hist( F_k    , bins =100,  density = True )
	ax100[3].hist( F_d    , bins =100,  density = True )
	ax100[4].hist( F_a    , bins =100,  density = True )
	ax100[5].hist( F_c    , bins =100,  density = True )
	ax100[6].hist( F_rho  , bins =100,  density = True )
	ax100[7].hist( F_psi  , bins =100,  density = True )
	ax100[8].hist( F_R0   , bins =100,  density = True )
	ax100[9].hist( F_RT   ,bins =100,   density = True )
	#plt.tight_layout() 
	#plt.savefig( OPATH_plots + 'FlatHist_Param_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	#plt.close()
	'''
	f10 = plt.figure( 10,  figsize=( 6,4 ))
	sns.violinplot(data=sam[['Lam' , 'bet' ,  'k' ]], orient="V")
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'vioLog_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()
	f11 = plt.figure( 11, figsize=( 6,4 ))
	sns.violinplot(data=sam[[ 'd' , 'a' , 'c']], orient="V")
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'vioRates'  + str( topVn ) + '_' + str( type ) +  '.pdf'  , dpi= 300  )
	plt.close()
	'''
	f12 = plt.figure(  12, figsize=( 6,4 ))
	plt.hist( sam.loc[:,"Score"]  ,  density = True )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'endScore_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()

	cmap  = plt.get_cmap("copper" )
	cMAP  = cmap( np.linspace( 0 , 1,  len(sam.loc[:,"R0"])  ) )
	f13 = plt.figure( 13 , figsize=( 6,4 ))
	#f14 = plt.figure( 14 , figsize=( 6,4 ))

	cc = 0
	for ik in range( 0 , len(sam.loc[:,"R0"]) ):
		plt.figure(13)
		plt.plot( sam.loc[ik,"R0"]  , sam.loc[ik,"PerInf"] , '*' , c=cMAP[cc] , markersize = 8 )
		plt.figure(14)
		plt.plot( sam.loc[ik,"R0"]  , sam.loc[ik,"RT"] , '*' , c=cMAP[cc] , markersize = 8 )
		cc = cc + 1
	plt.figure(13)
	plt.xlabel('R0' ) 
	plt.ylabel('% Infected cells' )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'R0-Inf_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()


	plt.figure(14 , figsize=( 6,4 ) )
	plt.xlabel('R0' ) 
	plt.ylabel('RT' )
	#plt.xlim( [ 0 , 5 ])
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'R0-RT'  + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()



def plot_thetaHist( sam ,  varNo , topVals  , thetatype):
	flatWalkers = sam.reshape(-1, sam.shape[-1])
	f1 , ax1 = plt.subplots(  varNo , 1 , figsize=( 10 , 8 ) )

	for fh in range( 0 , varNo ):
		ax1[fh].plot( sam[:,:,fh] , color ="gray" , alpha = 0.1 )
		ax1[fh].plot( np.mean(sam[:,:,fh], axis = 1 ) , color ="k" , alpha = 0.3 )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'ThetaTraj_S' + str( samNo ) + '_' + str( thetatype ) + '.png'  , dpi= 150  )
	plt.close()	
	# ============================================================
	# CORNER PLOTS:
	# ============================================================
	f10    = plt.figure( figsize =(6,4) )
	#labels = [ 'lam' , 'beta' , 'd', 'k' , 'a' , 'c' , 'D' , 'P' ]
	labels = [ r"$\log \lambda$" , r"$\log$k"  , 'd' , r"$\log$n" , r"$\delta$"  , 'c' , 'D' , r"$\log P$"  ]
	print( flatWalkers.shape )
	flatWalkers = flatWalkers[ flatWalkers[:,6] < 0.2 ] 
	print( flatWalkers.shape )
	f10    = corner.corner( flatWalkers ,show_titles=True,labels=labels,  max_n_ticks= 3 , color="cornflowerblue",
	                                     plot_datapoints=False, fill_contours=True , plot_contours = 1, 
	                                     quantiles= [0.16, 0.5, 0.84],
	                                     title_fmt = '.2f' , labelpad = 0.2 ) 
	plt.tight_layout()
	plt.savefig( OPATH_plots + 'cornerTreated_S' + str( samNo ) + '_' + str( thetatype ) +  '.pdf'  , dpi = 300 )
	plt.close()


def compute_statistics(column):
    median = np.nanmedian(column)
    lower_bound = np.nanpercentile(column, 2.5)
    upper_bound = np.nanpercentile(column, 97.5)
    return pd.Series([median, lower_bound, upper_bound], index=['Median', 'LB', 'UB'])


def calMeanCILog(  tim , yVals ):
	#print( "In cal Mean CI ")
	yVals = np.array( yVals)
	yVals[ yVals < LOD_siv ] = LOD_siv
	sdNum   = 2
	yVals   = yVals#np.nan_to_num( np.log10( yVals ) , neginf=0  )	
	yMean   = np.mean( yVals , axis = 0 )
	sd      = np.std(  yVals , axis = 0 )
	sem     = sd/np.sqrt( len( yVals ) )
	ci      = 1.96* sem
	sdL     = yMean - ( sdNum*sd )
	sdU     = yMean + ( sdNum*sd )
	#print( 'ci:' , ci , 'sd:', sd )
	return [ yMean,  sd ] #(yMean - ci ) , (yMean + ci ) , sdL , sdU ]


def calMeanCItipLog(  tim , yVals ):
	#print( "In cal Mean CI ")
	yVals = np.array( yVals)
	yVals[ yVals < LOD_tip] = LOD_tip
	sdNum   = 2
	yVals   = yVals#np.nan_to_num( np.log10( yVals ) , neginf=0  )	
	yMean   = np.mean( yVals , axis = 0 )
	sd      = np.std(  yVals , axis = 0 )
	sem     = sd/np.sqrt( len( yVals ) )
	ci      = 1.96* sem	
	return [ yMean, sd ] #sdNum*sd] # (yMean - ci ) , (yMean + ci ) , sdL , sdU ]

def calMeanCI(  tim , yVals ):
	#print( "In cal Mean CI ")
	sdNum   = 3
	yMean   = np.nanmean( yVals , axis = 0 )
	sd      = np.nanstd(  yVals , axis = 0 )
	sem     = sd/np.sqrt( len( yVals ) )
	ci      = sem #1.96* sem
	sdL     = yMean - ( sdNum*sd )
	sdU     = yMean + ( sdNum*sd )
	#sdL[ sdL < 0] = LOD_siv
	return [ yMean, sdNum*sd] # (yMean - ci ) , (yMean + ci ) , sdL , sdU ]

def calMedianCI(  tim , yVals ):
	#print( "In cal Mean CI ")
	sdNum   = 3
	yMean   = np.nanmedian( yVals , axis = 0 )
	sdL     = np.nanpercentile( yVals , 2.5)
	sdU     = np.nanpercentile( yVals , 97.5)
	#sdL[ sdL < 0] = LOD_siv
	return [ yMean, (yMean-sdL) , (yMean+sdU) ] # (yMean - ci ) , (yMean + ci ) , sdL , sdU ]


if __name__ == '__main__':
	global isCD4data
	global datatype
	global nwalkers , niter , yerr 
	global predicttime
	global tot_timPt
	global numBurn
	global cFWk2Day
	global StartRow    #   0 week should correspond to start
	global LOD_siv
	global LOD_tip
	global OPATH_plots
	global OPATH_files    
	global samNo
	global numRun
	global ipath
	global CutOffPerInf
	global CutTotalCD4
	global V2RNA
	V2RNA 				= 2
	If_TIPDat 			= 0
	datatype            = 'treated'
	isCD4data           = 1
	varNo               = 8	
	if datatype == 'control':
		StartRow    = 0
	else:
		StartRow    = 1
	if datatype == 'control':
		LOD_siv    = PARAM.LOD_siv
	else:
		LOD_siv    = PARAM.LOD_siv   
	LOD_tip        = LOD_siv #500    


	cFWk2Day    	 	= 7
	tot_timPt   		= 420 
	burnin 				= 500
	CutOffPerInf        = 100
	CutTotalCD4         = 10000  # cells/uL


	ipath  			= sys.argv[3]
	OPATH_plots 	= sys.argv[4] #'./output/plots/'
	OPATH_files 	= sys.argv[4] #'./output/files/'
	samNo     		= int( sys.argv[1] )
	totalRun  		= int( sys.argv[2] ) 
	numsteps  		= int( sys.argv[5] )
	TopFitNum       = int( sys.argv[6] )
	ViralV0         = PARAM.ViralV0
	ViralT0         = PARAM.ViralT0

	## Reading data
	AllDat = ExpData.get_dataInfantPVL(    datatype , ExpData.Animal2ExcludeInfant( datatype ) ) 
	# TIP- data
	If_TIPDat = 1
	TIPDat = []
	if If_TIPDat:
		TIPDat = ExpData.get_TIPdataInfant(    datatype , ExpData.Animal2ExcludeInfant( datatype ) )
		ylab  = TIPDat.columns.values.tolist()[1:]

	#Run mcmc
	#execute_mmc( AllDat  ,  varNo , numRun  )
	ModelObj        = tipmodel1()
	analysis( ModelObj , AllDat ,   TIPDat , totalRun , numsteps , burnin  , TopFitNum  , datatype )
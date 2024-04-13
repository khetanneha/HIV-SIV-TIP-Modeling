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
import corner
import sys
import seaborn as sns
plt.rcParams.update({'font.size': 18 })
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rc('legend', fontsize = 8 )
plt.tight_layout()




import extractExptData as ExpData
import compute_R0 as calR0
from GlobalparamInfantPVL import SysParam
PARAM = SysParam()



from tipmodel1_InfantPVL_A1 import tipmodel1





def analysis( MO, AllDat , TIPDat ,  nR , numsteps ,  burnin  , TopN , datatype ):
	colch = [ 'cornflowerblue', 'darkseagreen' ,'gray' , 'lightcoral']
	for ii in range( 0, samNo ):

		#  EXPT DATA
		ylab    = AllDat.columns.values.tolist()[1:]
		tmp 	 = AllDat.iloc[  : , [ 0, ii+1 ] ].dropna()  
		backGrnd = AllDat.iloc[  0 , [ 0, ii+1 ] ] 
		xx  	 = np.array( tmp.iloc[ StartRow:,0 ] )*cFWk2Day
		yy  	 = np.array( tmp.iloc[ StartRow:,1 ] )	
		whoamI   = ylab[ samNo - 1]
		print(" Data-ID" , whoamI )
		# ===================== 
		# Account for LOD:
		
		# VIRAL LOAD
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

			#plt.semilogy( xx , yy )
			#plt.show()
		eXY = np.column_stack( (xx , yy ))  
		data = ( xx , yy )  

		eXY_tip = []
		if If_TIPDat:
			# tip data
			tmptip      = TIPDat.iloc[ : , [ 0, ii+1 ] ].dropna()  
			tipBackGrnd = TIPDat.iloc[ 0 , [ 0, ii+1 ] ]
			xtip  		= np.array( tmptip.iloc[ StartRow:,0 ] )*cFWk2Day
			ytip  		= np.array( tmptip.iloc[ StartRow:,1 ] )
				
			if ( tipBackGrnd[1] > 1.1 ):
				ytip[:] = ytip[:] - tipBackGrnd[1]
				indx2   = ytip[:] < LOD_tip
				if (indx2.any() ):
					ytip[ indx2 ] = LOD_tip
			#plt.semilogy( xtip , ytip )
			#plt.show()
			eXY_tip = np.column_stack( (xtip , ytip ))
		# ==================================================================================================
		AllWalkers = None
		AllLogProb = None

		FinalVal     = []
		sortedThetas = []
		sortedScores = []
		SuccessfulEvents = 0

		ipath = ipath0 + './s' + str(ii +1) + '/output/plots/'
		print( ipath )

		print( "No. of steps:", numsteps )
		for nfiles in range(  1 , nR +1   ):		
			filename     = ipath + "back_" + str( nfiles ) + ".hdf5"
			sampler      = emcee.backends.HDFBackend( filename  )
		

			# FILE EXISTS?
			try:
				tmp_samples     = sampler.get_chain(  discard=burnin , flat=False )
				tmp_samples     = tmp_samples[0:numsteps,:,:]

				tmp_logprob     = sampler.get_log_prob(discard=burnin, flat=False) 
				tmp_logprob     = tmp_logprob[0:numsteps,:]
				print( nfiles , tmp_samples.shape[0])
			
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
					
			else:
						
				AllWalkers             = np.append( AllWalkers , tmp_samples , axis = 1  )
				# Choose number of fits to pick
				tmpVal                 = select_topScores( tmp_samples , tmp_logprob , TopN  )		
				for tmpv in tmpVal:
					FinalVal.append( tmpv )
					
			SuccessfulEvents       = SuccessfulEvents + 1
			if SuccessfulEvents > 100:
				break

				
		# =============================================================	
		#print( type(FinalVal ))
		# Sorting by scores - max to min
		maxL 		   = sorted( FinalVal , key=lambda x: x[0] , reverse = True )
		FinalVal 	   = maxL

		#fineTime       = np.linspace( 0 , 365*2 , int(365*2/7.0)  )
		fineTime        = eXY[:,0] #np.linspace( np.min( eXY[:,0] ) , np.max( eXY[:,0] ) , tot_timPt  )
			


		#------------------------------------
		AllR0       = np.zeros( len(FinalVal)  )
		AllRT       = np.zeros( len(FinalVal)  )
			

		AllCD4Cells = []
		AllVL       = []
		vTotGag     = []
		vTotTIP     = []
		AllTIP      = [] 
		AllHIV      = []
		ATTipCells =[]
		HHvir = []
		TTvir = []
		HTvir = []

		PerInfected = np.zeros( len(FinalVal)  )
		EndTotalCD4 = np.zeros( len(FinalVal)  )

		OutData    	= [] 
		SelectData 	= []

		f1 = plt.figure( 1 ,  figsize=( 6,4 ))
		f2 = plt.figure( 2 ,  figsize=( 6,4 ))
		f3 = plt.figure( 3 ,  figsize=( 6,4 ))
		f4 = plt.figure( 4 ,  figsize=( 6,4 ))
		f5 = plt.figure( 5 ,  figsize=( 6,4 ))
		cmap  = plt.get_cmap("copper")
		ccmap = cmap( np.linspace( 0 , 1,  len(FinalVal) ) )
		CF 	  = 10**(-3) # from mL to uL
		cc 	  = 0
		

		tmpScore = 0
		for ik in range( 0 , len( FinalVal ) ) :
			THETA   = FinalVal[ik][1]
			IV0     = MO.get_AllIV( )
			
		    # ==================================================
			intgdTIP        = ExpData.get_dataInfantIntegratedTIP( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
			intgdTIPSam     = intgdTIP.iloc[  0 , samNo  ]		
			# IC - 1 
			#iv                   = MO.set_AllIV(  V0=ViralV0 , VT0=0 , IT10=2.5*10**4 ,  TotalT0=PARAM.TotalT0 ) # 
			#iv                   = M0.set_AllIV(  V0=ViralV0 , VT0=0 , IT10=0.5*eXY_tip[0,1] ,  TotalT0=TotalT0 ) # 
			iv                    = M0.set_AllIV(  V0=ViralV0 , VT0=0 , IT10=intgdTIPSam ,  TotalT0=TotalT0 ) # 		

			
			
			new_mdl =  MO.model( THETA ,  IV0  ,  fineTime )

		
			if TIPModel == 1:
				HealthyCells   = ( new_mdl.y[0,:] + new_mdl.y[3,:] )*CF
				UnHealthyCells = ( new_mdl.y[1,:] + new_mdl.y[4,:] )*CF  
				TIPCells       =  new_mdl.y[3,:]*CF

			elif TIPModel == 3:
				HealthyCells   = ( new_mdl.y[0,:] + new_mdl.y[3,:] + new_mdl.y[4,:] + new_mdl.y[5,:] )*CF
				UnHealthyCells = ( new_mdl.y[1,:] + new_mdl.y[6,:] + new_mdl.y[7,:] + new_mdl.y[8,:] )*CF  
				TIPCells       =  ( new_mdl.y[3,:] + new_mdl.y[4,:] + new_mdl.y[5,:] ) *CF
			


			TotalCD4       	   = ( HealthyCells + UnHealthyCells ) 
			percentInfected    = ( UnHealthyCells/TotalCD4 )*100
		
			PerInfected[ik]    =  percentInfected[ -1 ]
			EndTotalCD4[ik]    =  TotalCD4[ -1 ]


			if TIPModel == 1:
				vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[5,:] )
				vTotTIP            =  v2RNA*new_mdl.y[5,:]
				vTotHIV            =  v2RNA*new_mdl.y[2,:]
			elif TIPModel ==3:
				vTotGag            =  v2RNA*( new_mdl.y[2,:] + new_mdl.y[9,:] )
				vTotTIP            =  v2RNA*new_mdl.y[9,:] 			
				vTotHIV            =  v2RNA*new_mdl.y[2,:]

	
			if ( (np.max( TotalCD4 ) < CutTotalCD4 )  and ( percentInfected[ -1 ] < CutOffPerInf ) and ( calR0.estimate_RTIP( THETA )>1) ):
				#print( np.max( TotalCD4 ))
				AllCD4Cells.append( TotalCD4 )
				ATTipCells.append( TIPCells/TotalCD4 )
				AllVL.append(  vTotGag )
				AllTIP.append( vTotTIP )
				AllHIV.append( vTotHIV )
				HHvir.append( new_mdl.y[2,:] )
				TTvir.append( new_mdl.y[5,:] )
				       

				# ========================================================================================
				AllR0[ik]         =  calR0.estimate_RR( THETA )   
				AllRT[ik]         =  calR0.estimate_RTIP( THETA )   
				# ========================================================================================

				t1 = np.array( FinalVal[ik][1] ).tolist()
				OutData.append( [ FinalVal[ik][0] , t1[0] , t1[1] ,t1[2] , t1[3] , t1[4] , t1[5] ,  t1[6] , t1[7] ,    AllR0[ik] , PerInfected[ik] , EndTotalCD4[ik] ,  AllRT[ik] ])




				if FinalVal[ik][0] < tmpScore:  

					#print( "Score:", tmpScore )
					plt.figure( f1 )
					plt.plot( new_mdl.t , vTotGag , 'o-',   c=ccmap[cc] ,   linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )  


					plt.figure( f3 )		
					plt.plot( new_mdl.t ,  TotalCD4 , 'o-',   c=ccmap[cc] ,   linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )

					plt.figure( f4 )	
					plt.plot( new_mdl.t ,  percentInfected  , 'o-',  c=ccmap[cc] ,  linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )

					plt.figure( f5 )
					plt.plot( new_mdl.t ,  vTotTIP , 'o-',  c=ccmap[cc] ,  linewidth = 0.1 , alpha = 0.8 ,  markersize = 4 )  

					SelectData.append( [ FinalVal[ik][0] , t1[0] , t1[1] ,t1[2] , t1[3] , t1[4] , t1[5] , t1[6] , t1[7] ,   AllR0[ik] , PerInfected[ik] , EndTotalCD4[ik] ,  AllRT[ik] ])
					tmpScore = FinalVal[ik][0]


			cc = cc + 1
		



		# Plotting Confidence Interval
		datMean , datL , datU , sdLow, sdUp = calMeanCILog( new_mdl.t ,  AllVL    )
		f1000 = plt.figure( 1000 ,  figsize=( 6,4 ))	
		plt.plot( new_mdl.t , datMean , color = colch[ii] , linewidth=2 , label = 'fit')
		plt.fill_between( new_mdl.t, sdLow , sdUp , color = colch[ii] , alpha = 0.2  , label = str(ii+1) )
		#plt.fill_between( new_mdl.t, datL , datU , color = 'red' , alpha = 0.2  )
		plt.plot( eXY[:,0] , np.log10( eXY[:,1] ) , 'o',  color = colch[ii] , linewidth = 0.8, label = 'expt', markersize=5 )
		#plt.yscale('log' )
		plt.xlabel('Days' ) 
		plt.ylim( 10**-1 , 10**9 )
		plt.ylabel('Log$_{10}$Gag-RNA copies/ml' ) 
		plt.tight_layout()	
		#plt.legend()

		
		# Plotting Confidence Interval
		datMean , datL , datU , sdLow, sdUp= calMeanCILog( new_mdl.t , AllTIP    )
		f1001 = plt.figure( 1001,  figsize=( 6,4 ))	
		plt.plot( new_mdl.t , datMean , color = colch[ii] , linewidth=2 )
		plt.fill_between( new_mdl.t, sdLow , sdUp ,  color = colch[ii] , alpha = 0.2 , label = str(ii+1) )
		#plt.fill_between( new_mdl.t, datL , datU , color = 'red' , alpha = 0.2  )
		if If_TIPDat:
			plt.plot( eXY_tip[:,0] , np.log10( eXY_tip[:,1] ) , 'o',  color = colch[ii] , linewidth=2 , markersize=5)
		#plt.yscale('log' )
		plt.xlabel('Days' ) 
		plt.ylim( 10**-1 , 10**9 )
		plt.ylabel('Log$_{10}$ TIP-RNA copies/ml') 
		plt.tight_layout()	
		#plt.legend()

		#%$$$$$$$$$$$$$
		f10001 = plt.figure( 10001,  figsize=( 6,4 ))	
		datMean , datL , datU , sdLow, sdUp= calMeanCILog( new_mdl.t , HHvir    )
		plt.plot( new_mdl.t , datMean , color = 'r' , linewidth=2 )

		datMean , datL , datU , sdLow, sdUp= calMeanCILog( new_mdl.t , TTvir    )
		plt.plot( new_mdl.t , datMean , color = 'b' , linewidth=2 )

		datMean , datL , datU , sdLow, sdUp= calMeanCILog( new_mdl.t , HTvir    )
		plt.plot( new_mdl.t , datMean , color = 'k' , linewidth=2 )
		'''
		plt.plot( eXY[:,0] , np.log10( eXY[:,1] ) , 'o',  color = 'r' , linewidth = 0.8, label = 'expt', markersize=5 )
		if If_TIPDat:
			plt.plot( eXY_tip[:,0] , np.log10( eXY_tip[:,1] ) , 'bo'  , linewidth=2 , markersize=5)
		'''
		#plt.yscale('log' )
		plt.ylim( 10**-1 , 10**9 )
		plt.xlabel('Days' ) 
		plt.ylabel('Log$_{10}$ virions/ml') 
		plt.tight_layout()	
		#plt.legend()		
		# Plotting Confidence Interval
		datMean , datL , datU , sdLow , sdUp = calMeanCI( new_mdl.t ,  AllCD4Cells  )
		f1002 = plt.figure( 1002 ,  figsize=( 6,4 ))
		plt.plot( new_mdl.t   , datMean , color = colch[ii] )	
		plt.fill_between( new_mdl.t, sdLow  , sdUp ,  color = colch[ii] , alpha = 0.2   )
		if isCD4data:
			cdData = ExpData.get_cd4dataInfant( datatype , ExpData.Animal2ExcludeInfant( datatype ) )
			cdxy   = ( cdData.iloc[:,[0 , ii+1 ]].dropna()  ) 
			plt.plot(  np.array( cdxy.iloc[:,0] )*cFWk2Day , np.array( cdxy.iloc[:,1] ) , 'o', color = colch[ii] , markersize=5  )  

		#plt.fill_between( new_mdl.t, datL , datU , color = 'red' , alpha = 0.2 )
		plt.xlabel('Days' ) 
		plt.ylabel('CD4+ cells/uL' ) 
		plt.tight_layout()	
		#plt.legend()
		#plt.savefig( OPATH_plots +'sdCD4_S' + str( samNo ) + '.pdf' , dpi = 300 )
		#plt.close() 
		

		# Plotting Confidence Interval
		datMean , datL , datU , sdLow , sdUp = calMeanCI( new_mdl.t ,  ATTipCells  )
		f1003 = plt.figure( 1003 ,  figsize=( 6,4 ))
		plt.plot( new_mdl.t   , datMean , color = colch[ii] )	
		plt.fill_between( new_mdl.t, sdLow  , sdUp ,  color = colch[ii] , alpha = 0.2   )
		plt.ylim( 10**-1 , 10**9 )	
		#plt.fill_between( new_mdl.t, datL , datU , color = 'red' , alpha = 0.2 )
		plt.xlabel('Days' ) 
		plt.ylabel('% TIP integrated cells' ) 
		plt.tight_layout()	
		#plt.legend()
		#plt.savefig( OPATH_plots +'sdCD4_S' + str( samNo ) + '.pdf' , dpi = 300 )
		#plt.close() 
	

	plt.figure(f1000)
	plt.savefig( OPATH_plots +'sdVL' + str( samNo )  + '.pdf' , dpi = 300 )
	#plt.close() 

	plt.figure(f1001)
	plt.savefig( OPATH_plots +'sdTIPVL' + str( samNo )  + '.pdf' , dpi = 300 )
	#plt.close() 
		
	plt.figure(f1002)
	plt.savefig( OPATH_plots +'sdcd4' + str( samNo )  + '.pdf' , dpi = 300 )
	#plt.close() 

	plt.figure(f1003)
	plt.savefig( OPATH_plots +'integratedTip' + str( samNo )  + '.pdf' , dpi = 300 )
	#plt.close() 
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
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'FlatHist_Param_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()


	f10 = plt.figure( 10,  figsize=( 6,4 ))
	sns.violinplot(data=sam[['Lam' , 'bet' ,  'k' ]], orient="V")
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'vioLog_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()

	'''
	f11 = plt.figure( 11, figsize=( 6,4 ))
	sns.violinplot(data=sam[[ 'd' , 'a' , 'c']], orient="V")
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'vioRates'  + str( topVn ) + '_' + str( type ) +  '.pdf'  , dpi= 300  )
	plt.close()
	'''


	f12 = plt.figure(  12, figsize=( 6,4 ))
	plt.hist( sam.loc[:,"Score"]  ,  density = True )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'Score_S' + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()


	cmap  = plt.get_cmap("copper")
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


	plt.figure(14)
	plt.xlabel('R0' ) 
	plt.ylabel('RT' )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'R0-RT'  + str( samNo ) + '_' + str( Sel_type ) +  '.pdf'  , dpi= 300  )
	plt.close()



def plot_thetaHist( sam ,  varNo , topVals  , thetatype):
	flatWalkers = sam.reshape(-1, sam.shape[-1])

	f1 , ax1 = plt.subplots(  varNo , 1 , figsize=( 12 , 10 ) )
	for fh in range( 0 , varNo ):

		ax1[fh].plot( sam[:,:,fh] , color ="gray" , alpha = 0.1 )
		ax1[fh].plot( np.mean(sam[:,:,fh], axis = 1 ) , color ="k" , alpha = 0.3 )
	plt.tight_layout() 
	plt.savefig( OPATH_plots + 'ThetaTraj_S' + str( samNo ) + '_' + str( thetatype ) + '.pdf'  , dpi= 300  )
	plt.close()




	# ============================================================
	# CORNER PLOTS:
	# ============================================================
	f10    = plt.figure( figsize =(6,4) )
	labels = ['lam' , 'beta' , 'd', 'k' , 'a' , 'c' , 'D' , 'P' ]
	f10    = corner.corner( flatWalkers ,show_titles=True,labels=labels,plot_datapoints=True,quantiles= [0.25, 0.5, 0.75]) 
	plt.tight_layout()
	plt.savefig( OPATH_plots + 'cornerTreated_S' + str( samNo ) + '_' + str( thetatype ) +  '.pdf'  , dpi = 300 )
	plt.close()


def calMeanCILog(  tim , yVals ):
	#print( "In cal Mean CI ")
	sdNum   = 3
	yVals   = np.nan_to_num( np.log10( yVals ) , neginf=0  )	
	yMean   = np.mean( yVals , axis = 0 )


	sd      = np.std(  yVals , axis = 0 )
	sem     = sd/np.sqrt( len( yVals ) )
	ci      = 1.96* sem

	sdL     = yMean - ( sdNum*sd )
	sdU     = yMean + ( sdNum*sd )
	#sdL[ sdL < 0] = LOD_siv
	return [ yMean,  (yMean - ci ) , (yMean + ci ) , sdL , sdU ]

def calMeanCI(  tim , yVals ):
	#print( "In cal Mean CI ")
	sdNum   = 3
	yMean   = np.mean( yVals , axis = 0 )
	sd      = np.std(  yVals , axis = 0 )
	sem     = sd/np.sqrt( len( yVals ) )
	ci      = 1.96* sem

	sdL     = yMean - ( sdNum*sd )
	sdU     = yMean + ( sdNum*sd )
	#sdL[ sdL < 0] = LOD_siv
	return [ yMean,  (yMean - ci ) , (yMean + ci ) , sdL , sdU ]


if __name__ == '__main__':

	datatype            = PARAM.datatype
	varNo               = PARAM.varNo
	TIPModel  	   		= PARAM.TIPModel      #1     # 1: 2003, 3: 2011
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
	OPATH_plots 		= sys.argv[3] #'./output/plots/'
	OPATH_files 		= sys.argv[3] #'./output/files/'
	#======================      
	samNo  				= int( sys.argv[1] )
	numRun 				= int( sys.argv[2] ) 
	np.random.seed( numRun )
	print("In Infant TIP PVL.....")





	ViralV0             = PARAM.ViralT0
	burnin 				= 0
	CutOffPerInf        = 100
	CutTotalCD4         = 10000  # cells/uL



	samNo     		= int( sys.argv[1] )
	totalRun  		= int( sys.argv[2] ) 

	ipath0  			= sys.argv[3]
	OPATH_plots 	= sys.argv[4] #'./output/plots/'
	OPATH_files 	= sys.argv[4] #'./output/files/'

	numsteps  		= int( sys.argv[5] )
	TopFitNum       = int( sys.argv[6] )






	## Reading data
	AllDat = ExpData.get_dataInfantPVL(    datatype , ExpData.Animal2ExcludeInfant( datatype ) ) 
	


	# TIP- data
	TIPDat = []
	if If_TIPDat:
		TIPDat = ExpData.get_TIPdataInfant(    datatype , ExpData.Animal2ExcludeInfant( datatype ) )
		ylab  = TIPDat.columns.values.tolist()[1:]

	Mobj = tipmodel1()
	analysis( Mobj , AllDat ,   TIPDat , totalRun , numsteps , burnin  , TopFitNum  , datatype )


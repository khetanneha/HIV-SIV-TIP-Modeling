'''
Neha Khetan,
Dec 2022 - Dec 2023,
TIP-model	
'''

from scipy.integrate import solve_ivp
from scipy import integrate, optimize
import numpy as np
from GlobalparamInfantPVL import SysParam
PARAM = SysParam()




class tipmodel1:
	"""Class for TIP-model defines all the variables, parameters and model
	   Sets up the mcmc-model """

	def __init__( self ):
		"""Define all the model-parameters and variables based on the TIP-model of choice .."""

		self.LOD_siv = PARAM.LOD_siv
		self.v2RNA   = PARAM.v2RNA
		self.labels  = ['lam' , 'k' , 'd', 'n' , 'delta' , 'c' , 'psi' , 'rho']
		# units are in mL, days
		self.lam     	= np.log10( 10**PARAM.Lam )               # birth rate of T
		self.k    	    = np.log10( 10**PARAM.Infectivity  )      # infection rate: T -> I
		self.d       	= PARAM.HealthyDEATH                      # death rate of T
		self.n       	= np.log10( 10**PARAM.Burst )             # bust size
		self.delta      = PARAM.delta                             # death rate of I
		self.c       	= PARAM.c                                 # death rate of V       
		self.psi       	= PARAM.psi
		self.rho       	= np.log10( 10**PARAM.rho  )   
		self.params     = [ self.lam , self.k , self.d , self.n , self.delta , self.c , self.psi , self.rho ]
		# Initial Values
		self.cf    	=  PARAM.CFul2ml 
		self.VT0   	=  0          					# Assuming 1 V has 2 RNA so Gag[0]/2
		self.V0     =  0     		
		self.It0    =  0
		self.I0     =  0  
		self.Id0    =  0
		self.T0     =  ( PARAM.TotalT0*self.cf ) - self.I0 - self.Id0 -self.It0
		self.AllIV  =   [ self.T0 , self.I0 , self.V0 , self.It0 , self.Id0 ,  self.VT0 ]		
		return

	def get_TIPmodel_param( self ):
		return self.params

	def get_TIPmodel_iv( self ):
		self.AllIV =   [ self.T0 , self.I0 , self.V0 , self.It0 ,  self.Id0 ,  self.VT0 ]
		return [ self.T0 , self.I0 , self.V0 , self.It0 , self.Id0 ,   self.VT0 ]

	def get_AllIV( self ):
		return self.AllIV

	def set_AllIV( self , V0 ,  VT0 , IT10 , TotalT0  ):
		# All conc. are in cell/mL hence, conversion
		self.cf     =  PARAM.CFul2ml
		self.VT0    =  VT0          # Assuming 1 V has 2 RNA so Gag[0]/2
		self.V0     =  V0     		
		self.I0     =  0 
		self.It0    =  IT10          # GFP-RNA copies/2
		self.Id0    =  0
		self.T0     =  ( (TotalT0 )* self.cf ) - self.I0 - self.Id0  - self.It0 
		self.AllIV  =   [ self.T0 , self.I0 , self.V0 , self.It0 , self.Id0 , self.VT0 ]
		return self.AllIV

	def TIPmodel_ode( self, t , iv , par ):
		"""ODE model for TIP: with populations:
			T : Healthy T-cells
			I : Infected T-cells
			V : Viral population
			It: TIP-infected population
			Id : Dually infected population
			VT : TIP-population
			Model parameters: P and D 
			Here: model is as in Weingberger et.al., JVI 2003
		"""
		T , I , V , It ,   Id ,  VT = iv
		lam1 ,  k1 , d , n1 , delta , c ,  D , rho1        = par
		# lam, beta , d , k , a      , c , D    , P
		lam   	   = 10**lam1
		k     	   = 10**k1
		n     	   = 10**n1
		P   	   = 10**rho1

		T_dt       = lam - T*( d + (k * ( V + VT ) ) )
		I_dt       = ( k * T * V ) - ( delta * I )
		V_dt       = ( n * delta * I) + ( D*n*(D * delta )*Id ) - ( c * V )
		It_dt      = ( k * T * VT )   - ( It * d ) - ( It * k * V )
		Id_dt      = ( It * k * V )   - ( (D *  delta ) * Id )
		VT_dt      =  ( (P**2) * D * n *(D *  delta ) *Id ) - ( c * VT )                    
		return [ T_dt , I_dt , V_dt , It_dt , Id_dt , VT_dt ]

	def model( self , theta , iv , x  ):
		"""Solve ODE"""
		ss = integrate.solve_ivp( self.TIPmodel_ode ,  ( x[0] , x[-1] )  , iv , args= (theta,)  , t_eval = x , method = "BDF" )
		while( ss.success == 0 ):
			ss = integrate.solve_ivp( self.TIPmodel_ode ,  ( x[0] , x[-1] )  , iv , args= (theta,)  , t_eval = x , method = "BDF" )
		return ss

	def lnlike( self , theta, x, y, yerr ):
		"""Calculate Log-likelihood"""
		yM          = self.model(  theta, self.get_TIPmodel_iv(), x )
		# Viral population; RNA copies will be twice
		VL          = self.v2RNA*( yM.y[2,:] + yM.y[5,:] )	
		indx = ( VL[1:] < self.LOD_siv )
		if ( indx.any()  ):
			VL[:] = self.LOD_siv
		## Log-scale
		yMLog 		= np.nan_to_num( np.log10( VL ) , neginf=0  )
		yELog       = np.log10( y )		
		ssd   		= ( yMLog  - yELog )**2  
		sum_measure = np.sum( ssd)/PARAM.MLsigma**2 
		LnLike  	= -0.5*sum_measure -( 0.5*np.log( 2*np.pi*PARAM.MLsigma**2 )*len( yELog ) )
		return LnLike

	def lnprior( self , theta):
		"""
		This module: corresponds to the TIP-model, as in Weinberger et al. 2003	
		Priors are based on fits to the untreated animals
		In TIP-fits: P and D: rest are uninformed priors; rest as
		Informed priors : estimated from the posterior distributions obtained fits to the un-treated animals, 
		with decay rates are from literature : informed priors in untreated animals with rest as
		uninformed priors
		Bounds are estimated from  posterior distributions of MCMC-fits on un-treated animals
		"""
		#lam , beta , d, k , a , c , D , P = theta
		lam ,  k , d , n , delta , c ,  psi , rho = theta
		mu1, sigma1 =   PARAM.Lam   		 , PARAM.LamSD  
		mu2, sigma2 =	PARAM.Infectivity    , PARAM.InfectivitySD 
		mu3, sigma3 =   PARAM.HealthyDEATH   , PARAM.HealthyDEATHSD  
		mu4, sigma4 =   PARAM.Burst   		 , PARAM.BurstSD  
		mu5, sigma5 =   PARAM.delta    	     , PARAM.deltaSD  
		mu6, sigma6 =   PARAM.c    		     , PARAM.cSD             
		'''
		bounds on parameters from  posterior distributions of MCMC-fits on un-treated animals
		'''
		c1 = ( 4 < lam  < 6 )
		c2 = ( -9 < k < -6 ) 
		c3 = ( delta >= d > 0 ) 
		c4 = ( 1 < n <  6 )
		c5 = ( delta >= d > 0 ) 
		c6 = ( 36 >= c >= 9 )   
		c7 = ( 0.01 < psi <= 1 )
		c8 = ( 0 < rho <= 2 ) 
		c9 = delta >= psi*delta > d
		tot_c = (   np.log( 1.0/(np.sqrt(2*np.pi)*sigma1) ) -0.5*(lam-mu1)**2/sigma1**2 +
					np.log( 1.0/(np.sqrt(2*np.pi)*sigma2) ) -0.5*(k -mu2)**2/sigma2**2 +
					np.log( 1.0/(np.sqrt(2*np.pi)*sigma3) ) -0.5*(d-mu3)**2/sigma3**2 +
					np.log( 1.0/(np.sqrt(2*np.pi)*sigma4) ) -0.5*(n-mu4)**2/sigma4**2 +
					np.log( 1.0/(np.sqrt(2*np.pi)*sigma5) ) -0.5*(delta-mu5)**2/sigma5**2 +
					np.log( 1.0/(np.sqrt(2*np.pi)*sigma6) ) -0.5*(c-mu6)**2/sigma6**2  )	 
		if not ( c7 and c8 and c1 and c2 and c3 and c4 and c5 and c6 and c9 and tot_c) :
			return -np.inf
		return tot_c

	def lnprob( self , theta , x , y, yerr ):
		"""computes log probablity for sampling"""
		lp = self.lnprior(  theta )
		if not np.isfinite(lp): 
			return -np.inf
		scoreVal = self.lnlike(  theta, x, y, yerr ) 
		return lp + scoreVal

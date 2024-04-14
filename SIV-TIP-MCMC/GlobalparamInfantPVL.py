'''
Neha Khetan, 2022-2023, 
TIP 1:  no super-infection
TIP 2:  clearance rates
TIP 3:  superinfection m = 3
TIP 4 : super-infection m =2
TIP 5 : heterozygous | Assuming they do not infect new cells
TIP 6 : Cell division and growth
TIP 7:  immature virions that do not infect: with and w/o feedback on infectivity 
'''


class SysParam:
	"""Class for defining system parameters"""

	def __init__( self ):
		"""Defines variables associated with experimenta-data type 
				(i.e. SHIV infected only or SHIV infected + TIP-treated), system, method, and
		   simulation parameters for MCMC runs, priors  etc., and,
		   simulation parameters for ODE TIP-models"""
		self.numBurn      =  50
		self.niter        =  50		
		self.TIPModel  	  =  1     
		self.AnimalType   = 'INFANT'
		self.datatype     = 'treated'
		self.method       = 'PVL'	
		# INITIAL PRIORS
		self.psi           = 0.1
		self.rho           = 1
		self.Lam           = 5
		self.Infectivity   = (-7.5)
		self.Burst         = 3.5
		self.HealthyDEATH  = 0.02 
		self.delta         = 0.88  
		self.c             = 23  
		self.LamSD         = 0.5
		self.InfectivitySD = 0.5 
		self.BurstSD       = 0.5
		self.HealthyDEATHSD = 0.005 
		self.cSD            = 4
		self.deltaSD        = 0.4     
		self.gamma          = self.delta
		self.h0            = 0.03
		self.hmax          = 3.3*10**6
		self.ViralV0      = 0.1 #0.5 #0.05
		self.ViralT0      = 0

		if self.datatype == 'control':
			self.StartRow   = 0
			self.varNo      = 6 
		else:
			if self.TIPModel == 2:
				self.StartRow   = 1
				self.varNo      = 9
			else:
				self.StartRow   = 1
				self.varNo      = 8
		self.MLsigma      = 1
		self.nwalkers 	  = 20
		self.CFul2ml      = 10**3
		self.CFml2ul      = 10**-3
		self.cFWk2Day     = 7
		self.predicttime  = 600 
		self.tot_timPt    = 300
		self.v2RNA 		    = 2
		if self.AnimalType == "INFANT":
			self.isCD4data    = 1
			self.If_TIPDat    = 1
			self.LOD_tip      = 500  
			self.TotalT0	  = 3000
			if self.method  == "PVL":
				self.LOD_siv   =  5 #0.01
			else:
				if self.datatype =='control':
					self.LOD_siv = 1
				else:
					self.LOD_siv   = 100
		else:
			self.isCD4data    = 0
			self.If_TIPDat    = 0
			self.LOD_tip      = 500 #10**3   
			self.TotalT0	  = 1500 
			self.LOD_siv      = 20
		self.epsilonIM1  = 0.5  # immature in TIP
		self.epsilonIM2  = 0.5  # immature in TIP
		self.epsilon     = 0.1  # ART and Inf= f( v )
		self.epsilon2    = 0.1  # ART
		'''
		if PARAM.DATASET == "Adult": 
			from tipmodel3_Adult import tipmodel1
		elif PARAM.DATASET == "Inf-PVL":
			from tipmodelOne_InfantPVL import tipmodel1
		elif PARAM.DATASET == "Inf-ddPCR":
			from tipmodelOne_InfantddPCR import tipmodel1
		'''

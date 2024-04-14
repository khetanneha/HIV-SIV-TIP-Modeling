'''
Neha Khetan, 2022-2023
computes R0 of wt-virus and TIP
'''
import numpy as np

def estimate_RR( new_theta_max ):
	"""Estimate R0 of WT-virus based on basic model of Viral dynamics"""
	num = (10**new_theta_max[0]) * (10**new_theta_max[1]) * (10**new_theta_max[3])
	den = ( new_theta_max[2] ) * ( new_theta_max[4] ) * ( new_theta_max[5] )
	return np.array( num/den )


def estimate_RTIP( new_theta_max ):
	"""Estimate R0 of TIP-virus based on expanded basic model of Viral dynamics i.e. TIP-model"""
	rWT    = estimate_RR( new_theta_max ) #2
	prefac = ( 10**new_theta_max[7]) * (10**new_theta_max[7] ) * new_theta_max[6]
	fac    = 1 - (  1/rWT  )
	return np.array( prefac * fac )

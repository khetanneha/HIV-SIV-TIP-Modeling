'''
Sep, 2022, Neha Khetan, 
Read in data from csv files: untreated and TIP_treated animals
'''

import numpy as np
import pandas as pd

"""Reads data from csv files for All animals: PVL ( ddPCT and qPCR ) , TIP-GFP (ddPCR) and Integrated-TIP, CD4+
"""

def get_dataAdult( dataset , anim2ignore ):
	# Reading in the data , ignore if any data
	if dataset == 'control':
		Ctr = pd.read_csv( './data/AdultCTR_G1_20wk_filter.csv' , sep = ',' ,  skiprows= None, header = 0 )     
		print("Dataset in use: control")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	else:
		Ctr = pd.read_csv( './data/AdultTIP_G1_20wk_filter.csv' , sep = ',' ,  skiprows= None, header = 0 )
		print("Dataset in use: Treated")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	return Ctr


def Animal2ExcludeAdult( dataset ):
	if dataset == 'control':
		animalSkip    =  [] #[ '40444']
	elif dataset == 'treated':
		animalSkip    = [] #[ '40407' , '40408']
	else:
		print("Choose right dataset.....")
	return animalSkip

def get_TIPdataInfant( dataset , anim2ignore ):
	# Reading in the data , ignore if any data
	if dataset == 'control':
		tipInfo = pd.read_csv( './data/TIPGfpRNAControl.csv' , sep = ',' , skiprows= None, header = 0 )
		print("Dataset in use: control")
		if anim2ignore:
			tipInfo  = tipInfo.drop( anim2ignore , axis =1)
	else:
		#tipInfo = pd.read_csv( './data/TIPGfpRNATreated-Mod.csv' , sep = ',' , skiprows= None, header = 0 )
		#tipInfo = pd.read_csv( './data/SIV_GFP_TIP.csv' , sep = ',' , skiprows= None, header = 0 )
		tipInfo = pd.read_csv( './data/newtipdata_v2.csv' , sep = ',' , skiprows= None, header = 0 )

		print("Dataset in use: Treated")
		if anim2ignore:
			tipInfo   = tipInfo.drop( anim2ignore , axis =1)
	return tipInfo

def get_cd4dataInfant( datatype ,  anim2ignore ):
	if datatype == 'control':
		Cd4 = pd.read_csv( './data/cd4_ctrl.csv' , sep = ',' , skiprows= None, header = 0 )
		print("CD4 data from control")
		if anim2ignore:
			Cd4   = Cd4.drop( anim2ignore , axis =1)
	else:
		Cd4 = pd.read_csv( './data/cd4_treated.csv' , sep = ',' , skiprows= None, header = 0 )		
		print("CD4 data from treated")
		if anim2ignore:
			Cd4   = Cd4.drop( anim2ignore , axis =1)
	return Cd4

## IF INFANT DATA
def get_dataInfantPVL( dataset , anim2ignore ):
	# Reading in the data , ignore if any data
	if dataset == 'control':
		#Ctr = pd.read_csv( './data/SIVGagRNAControl.csv' , sep = ',' , skiprows= None, header = 0 )
		Ctr = pd.read_csv( './data/PVLInfant_CTR_10 TCID50.csv' , sep = ',' , skiprows= None, header = 0 )		
		print("Dataset in use: Infant PVL control")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	else:
		#Ctr = pd.read_csv( './data/SIVGagTreated_LOD_500skip.csv' , sep = ',' , skiprows= None, header = 0 )
		Ctr = pd.read_csv( './data/PVLInfant_TIP_10_TCID50.csv' , sep = ',' , skiprows= None, header = 0 )		
		print("Dataset in use: Infant PVL Treated")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	return Ctr

def get_dataInfantddPCR( dataset , anim2ignore ):
	# Reading in the data , ignore if any data
	if dataset == 'control':
		Ctr = pd.read_csv( './data/SIVGagRNAControl.csv' , sep = ',' , skiprows= None, header = 0 )
		print("Dataset in use: Infant ddPCR control")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	else:
		Ctr = pd.read_csv( './data/SIVGagTreated_LOD_500skip.csv' , sep = ',' , skiprows= None, header = 0 )
		print("Dataset in use: Infant ddPCR Treated")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	return Ctr

def Animal2ExcludeInfant( dataset ):
	if dataset == 'control':
		animalSkip    =  [ '40444']
	elif dataset == 'treated':
		animalSkip    = [] #[ '40407' , '40408']
	else:
		print("Choose right dataset.....")
	return animalSkip

## IF INFANT DATA
def get_dataInfantIntegratedTIP( dataset , anim2ignore ):
	# Reading in the data , ignore if any data
	if dataset == 'control':
		#Ctr = pd.read_csv( './data/SIVGagRNAControl.csv' , sep = ',' , skiprows= None, header = 0 )
		Ctr  = pd.read_csv( './data/PVLInfant_CTR_TIP.csv' , sep = ',' , skiprows= None, header = 0 )		
		print("Dataset in use: Infant PVL control")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	else:
		#Ctr = pd.read_csv( './data/SIVGagTreated_LOD_500skip.csv' , sep = ',' , skiprows= None, header = 0 )
		Ctr = pd.read_csv( './data/virTIP_Integrated.csv' , sep = ',' , skiprows= None, header = 0 )		
		print("Dataset in use: Infant PVL Treated for ")
		if anim2ignore:
			Ctr   = Ctr.drop( anim2ignore , axis =1)
	return Ctr
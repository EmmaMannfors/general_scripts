#################################################################################
#			Perform 2D FellWalker
#################################################################################
# 	Each instance of os.system() makes a new shell, and so Starlink-
#			commands must be run in each excecution
#################################################################################
#********************************************************************************
#	BEFORE RUNNING THIS SCRIPT: 
#	input the following commands into the terminal
#	( Change to your own FW path)

#	export STARLINK_DIR=/home/local/emannfor/star-2021A
#	source $STARLINK_DIR/etc/profile
#	convert
#	cupid
#********************************************************************************
#################################################################################

import os
from astropy.io import fits
import numpy as np
import aplpy
import matplotlib.pyplot as plt

#DIR	= '/home/local/emannfor/star-2021A/'
DIR	= '/PATH/TO/STARLINK/'

import sys
sys.path.append(DIR+'lib/')
os.environ['PATH'] += os.path.pathsep + DIR+'lib/'
os.environ['PATH'] += os.pathsep + DIR+'lib/'

#############################
# Init
#############################


cmd1 = "export STARLINK_DIR="+DIR+"star-2021A;"	# CHANGE star-2021A to whatever your version is
cmd2 = "source $STARLINK_DIR/etc/profile;"

convert = DIR+'bin/convert/convert.sh'
cupid = DIR+'bin/cupid/cupid.sh'


# Individual routines
fits2ndf = DIR+'bin/convert/fits2ndf'
ndf2fits = DIR+'bin/convert/ndf2fits'

smurf      =  DIR+"bin/smurf/smurf.csh"
kappa      =  DIR+"bin/kappa/kappa.csh"
makemap       = DIR+"bin/smurf/makemap"
makesnr       = DIR+"bin/kappa/makesnr"
findclumps    =  DIR+"bin/cupid/findclumps"
extractclumps =  DIR+"bin/cupid/extractclumps"
thresh        =  DIR+"bin/kappa/thresh"
findback      =  DIR+"bin/cupid/findback"
stats         =  DIR+"bin/kappa/stats"
ndfcopy       =  DIR+"bin/kappa/ndfcopy"
cmd3 = convert
cmd4 = cupid
init = cmd1+cmd2+cmd3+';'+cmd4+';'




#############################
# Used in get_rms()
#############################
def PointsInCirclesVecOrig(x,y,X,Y,R):	# Version from Mika!
    '''
    Check which points of a grid (x,y) are
    within a circle of radius R centered on point X,Y
    If Circles > 1, check several circles
    Assuming that the point is more likely outside the circle
    x,y = the point to check
    X,Y = the centre of the circle
    R   = the radius of the circe in pixels
    Returns a mask of where points within the circle are 1
    '''
    if np.ndim(X) == 0:
        X = np.asarray([X])
        Y = np.asarray([Y])
        R = np.asarray([R])
    GRID = np.zeros((len(x),len(y)))
    CC = len(R)
    
    if len(X) == len(Y) and len(X) == len(R):
        for N in range(CC):
            for i in range(len(x)):
                dx = abs(x[i]-X[N])
		#print dx, R[N]
                if dx < R[N]: # cut-off by x coordinate
                    for j in range(len(y)):
                        dy = abs(y[j]-Y[N])
                        if dy < R[N]:# cut-off by y coordinate
                            # We now have a square left
                            # cut it to a circle
                            if dx**2 + dy**2 <= R[N]**2:
                                GRID[i,j] = 1
        return GRID
    else:
        print('Mismatch in array lengths of circle coordinates')


#############################
# Calculate RMS
#############################
def get_rms(IN): 
	# Calculate rms within a "empty" reference region
	# rao0deg,dec0deg, R_deg are the RA,DEC,and radius (degrees) of the region
	
	ra0Deg,dec0Deg,R_deg	= 83.888555,-5.0379067,0.02
	
	f		= fits.open(IN)
	fData		= f[0].data
	
	pix_size	= abs(f[0].header['CDELT1'])
	R		= R_deg/pix_size

	fig		= aplpy.FITSFigure(IN)
	I,J		= np.indices(fData.shape)
	ra0,dec0	= fig.world2pixel(ra0Deg,dec0Deg)
	
	GRID		= PointsInCirclesVecOrig(I[:,0],J[0],ra0,dec0,R)	# radius needs to be in pix!

	m		= np.nonzero((GRID==1)&(np.isfinite(fData)))#&(fData>0.0)))			
	if fData[m].shape == (0,):
		STD = 'nan'
	else: 
		STD = np.nanstd(fData[m])

	return str(STD)




#############################
# Convert fits to sdf file
#############################
def makeSDF(IN):
	OUT	= IN.split('.fits')[0]+'.sdf'
	os.system(cmd1+cmd2+cmd3+';'+fits2ndf+' '+IN+' '+OUT)
	print('########## Written to:',OUT)
	return OUT

#############################
# Convert sdf to fits
#############################
def makeFits(FW_NAME,FITS_OUT):
	os.system(cmd1+cmd2+cmd3+';'+ndf2fits+' '+FW_NAME+' '+FITS_OUT)
	print('################# Written to ',OUT)


#############################
# FellWalker
#############################

def FW(SDFIN,rms,NAME): 
	CONF_FILE	= 'config_file.txt'
	fIN		= SDFIN				# SDF input file
	cat		= NAME+'_cat'			# cat file
	FW_OUT		= NAME+'_clumps'		# Output file
		
	command = init+findclumps+' '+fIN+' outcat= '+cat+' out='+FW_OUT+' method=FellWalker'+' RMS='+rms+' deconv=false config=^'+CONF_FILE

	os.system(command)
	print('################# Written to ',FW_OUT)
	return FW_OUT


#############################
# Extract clumps
#############################
def extract(FW_IN,FIELD,NAME): 
	data	= NAME+'.sdf'
	fileOUT	= NAME+'_out.sdf'
	cat	= NAME+'_cat'
	maskIN	= FW_IN
	command = init+extractclumps+' '+maskIN+' '+data+' '+fileOUT+' outcat='+cat+' deconv=true shape=ellipse2'
	os.system(command)
	print('################# Written to ',fileOUT)




#############################
# Run the FellWalker script
#############################
def RUN_FW(NAME,IN):
	rms	= 	get_rms(IN)
	SDFIN	= makeSDF(IN)
	print('Made SDF')
	FW_OUT	= FW(SDFIN,rms,NAME)
	print('Did FW')
	makeFits(FW_IN,FW_OUT)
	print('Made FITS')
	extract(FIELD,'FW/'+NAME)




	















""" This code is approximating the total cross section for the physics process
	of "Elastic neutrino-electron scattering" to its leading order"""
	
##########################################################################################
# Madhuranga Thilakasiri                                                                 #
# Final project for the course, Experimental Methods of High Energy Physics              #
# Instructor: Dr. Alexander Khanov                                                       #
# 12/11/2017	                                                                         #
##########################################################################################

print"\n\n###############################################################################################"
print "\n\v\tWelcome to the Monte-Carlo integration approach in estimating the total cross section"
print "\n\t\t\tfor the physics process enu -> Z -> enu \n"
print"\n#################################################################################################"

#------------------ Importing modules --------------------

import sys
import random
import numpy as np
import time
import matplotlib.pyplot as plt

#------------- Set seed for reproducibility ---------------

seed = 54327
np.random.seed(seed)


#---------------- Parameter definition --------------------

# all parameters are adopted from: 
# http://pdg.lbl.gov/2017/reviews/rpp2016-rev-phys-constants.pdf 
# and equations are from 
# Introduction to Elementary Particles by David Griffiths

pb_conv = 3.893793656E8     # GeV^-2 -> pb
alpha = 1./132.035999139    # QED alpha
sw2 = 0.2314                # sin^2(Weinberg angle) s.t. Weinberg angle = 28.75 degrees 
G = 1.1663787E-5            # Fermi coupling constant
pi = np.pi
#for electrons/muons/taus, C_V and C_A values are as follows
CV = -0.5 + 2*sw2
CA = -0.5

# lets consider the center of mass energy for the scattering is  = 80 GeV

CME = 80                    # in GeV
CME2 = CME**2

print "\nThe center of mass energy is = ", CME, "GeV\n"

#------------- End of parameter definition -----------------

start_time = time.time()

#------------ Definition of the function that approximate the total cross section ----------------

def dsigma(costheta):
		
	# refer eq 9.99 in page 333 in Introduction to Elementary Particles by David Griffiths 
	# and eq 11 in the report
	
	return (1./2*(pi)) * (G**2) * CME2 * ((CV+CA)**2 + (CV-CA)**2 * ((1 + costheta)**2/4))
	                   

#--------------------- End of function definition -----------------------	


#--------------------- Monte-Carlo Integration ---------------------------

# integration range is when costheta goes from -1 to +1 : (x_2 - x_1) in (x_2 - x_1)f(x)

rng = 2

# number of integral points

N = 1000000

#-------- loop over the phase space for sigma and variance calculations ----------------------

# define wsum and wsum_sq for sigma and variance

wsum = 0               # weight as defined in the report
wsum_sq = 0

# define maximum values for dsigma and costheta

wmax = 0
costheta_max = -2

trials = []
for k in range(0, N):
	trials.append(k)
val = []

print "integrating for the total cross-section approximation"
print "....."

for i in range(0,N):
	# printing progress
	sys.stdout.write("progress: %d%%   \r" % (float(i)*100./(N)) )
	sys.stdout.flush()
	
	# generating random values for costheta within -1 to +1
	costheta_i = -1 + random.random() * rng
	
	# calculating the point at the phase space
	w_i = dsigma(costheta_i) * rng
	
	# adding to the sum
	wsum = wsum + w_i
	wsum_sq = wsum_sq + wsum**2
	sigma = (wsum/N)*pb_conv
	val.append(sigma)
	
	# check if w_ii > wmax 
	if w_i > wmax:
		wmax = w_i
		costheta_max = costheta_i

#--------- cross section calculation --------------------

sigma = wsum/N


#--------- error calculation -------------------

variance = (wsum_sq/N) - (wsum/N)**2
error = np.sqrt(variance/N)

print "\nintegration completed...\n"
print "dsigma maximum: ", wmax, "found at costheta: ", costheta_max
print "\nTotal cross section approx: ", sigma*pb_conv, "+-", error, "pb"
tot = sigma*pb_conv
#---------------------------- end of the loop-----------------------------------------------


#------------------------- Analytical cross section calculation ----------------------------

# Refer to equation 9.100 in page 333 in Introduction to Elementary Particles by David Griffiths 
# and eq 13 in the report

an_sigma = (4. / 3*pi) * G**2 * CME2 * (CV**2 + CA**2 + CV*CA)
tot_an = an_sigma*pb_conv
print "\nAnalytical value for the cross section: ", an_sigma*pb_conv , "pb\n"
print "\nDiffernce (aprox_sigma / an_sigma): ", np.abs(sigma / an_sigma)
print "\n", N, "integral points were integrated in %s seconds" % (time.time() - start_time)


plt.axhline(tot_an, label="Analytical = %1.4f pb" % tot_an, color='g', linewidth='4')
plt.plot(trials, val, label="Estimation = %1.4f pb" % tot, linewidth='3', color='r')
plt.grid()
plt.title("Total cross section ($\sigma$) approx. for %d \nintegral points" % N)
plt.xlabel('# of Trials')
plt.ylabel('Estimated value for $\sigma$ pb')
plt.legend(loc='best')
plt.show()
print"\n+++++++++++++++++++++++++++++++++++++++++++ End of the code ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"


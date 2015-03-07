import numpy as np
import pylab as pl
import scipy.stats as stats


def velocityPlots():
	# Plotting of the solutions to the first 
	atoms = np.genfromtxt('../data/velocities.xyz',skip_header = 2)#, usecols = (0,1,2,3,4,5,6))

	nAtoms = atoms.shape[0]

	velocityx = sorted(atoms[:,1])
	velocityy = sorted(atoms[:,2])
	velocityz = sorted(atoms[:,3])

	# print velocityy

	fitx = stats.norm.pdf(velocityx, np.mean(velocityx), np.std(velocityx))  #this is a fitting indeed
	fity = stats.norm.pdf(velocityy, np.mean(velocityy), np.std(velocityy))  #this is a fitting indeed
	fitz = stats.norm.pdf(velocityz, np.mean(velocityz), np.std(velocityz))  #this is a fitting indeed

	print "std_x = " + str(np.std(velocityx)) 
	print "std_y = " + str(np.std(velocityy))
	print "std_z = " + str(np.std(velocityz))


	# h = sorted(velocityz)  #sorted
	velocityDistributions = pl.figure()


	x = pl.subplot(311)
	pl.title('Starting velocity distributions')

	pl.plot(velocityx,fitx,'-')
	pl.hist(velocityx,normed=True, bins=100)      #use this to draw histogram of your data

	y = pl.subplot(312)

	pl.plot(velocityy,fity,'-')
	pl.hist(velocityy,normed=True, bins=100)      #use this to draw histogram of your data

	z = pl.subplot(313)

	pl.plot(velocityz,fitz,'-')
	pl.hist(velocityz,normed=True, bins=100)      #use this to draw histogram of your data


	pl.savefig('velocityDistributions.png')

	pl.show()                   #use may also need add this 


def forceplot():
	sigma = 1.0
	epsilon = 1.0

	def lennardJonesForce(x):
		return 24/sigma*epsilon*(2*(sigma/x)**13 - (sigma/x)**7)

	xx = np.linspace(0.01,3,1000)


	forceFig = pl.figure()
	pl.plot( xx, lennardJonesForce(xx))
	pl.ylim( [-3,10] )
	pl.grid(True)
	pl.title('The force between to particles. sigma = 1, epsilon = 1')
	pl.xlabel('Distance')
	pl.ylabel('Force')

	pl.savefig('forcePlot.png')

def energyplots():
	# The coloumns are time, kinetic energy, potential energy, temperature in the  order
	stats =  np.genfromtxt("../statistics/stats", skip_header = 1)

	totalEnergy = stats[:,1] + stats[:,2]
	# print max(totalEnergy)
	# print min(totalEnergy)
	# print totalEnergy

	energyplots = pl.figure()
	

	total = pl.subplot(311)
	pl.title('Energy in the system')

	pl.plot(stats[:,0], totalEnergy)
	pl.ylabel("Total Energy")

	kinetic = pl.subplot(312)
	pl.plot(stats[:,0], stats[:,1])
	pl.ylabel("Kinetic Energy")


	potential = pl.subplot(313)
	pl.plot(stats[:,0], stats[:,2])
	pl.ylabel("Potential Energy")

	print max(stats[:,1]) - min(stats[:,1])
	print max(stats[:,2]) - min(stats[:,2])
	print max(stats[:,1] + stats[:,2]) - min(stats[:,1] + stats[:,2])
	
	pl.show()
	# print stats[:,0]
	

##
## What do I want to do this time...
##


velocityPlots()
# forceplot()
# energyplots()



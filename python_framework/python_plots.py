import numpy as np
import pylab as py
import scipy.stats as stats

def plot_function():

	# Plotting of the solutions to the first 
	statisticalValues = np.genfromtxt('../data/statisticalResults/statisticalValues.csv', delimiter='\t')

	timeToPicoSeconds = 1.00224e-13 / 10e-12
	temperatureToKelvin = 119.735
	energyToJoule = 1.65313e-21
	pressureToPa = 1.65313e+09

	time = timeToPicoSeconds * statisticalValues[1: , 0]
	temperature = temperatureToKelvin * statisticalValues[1: , 1]
	kineticEnergy = energyToJoule * statisticalValues[1: , 2]
	potentialEnergy =energyToJoule * statisticalValues[1: , 3]
	pressure = statisticalValues[1: , 5]

	totalEnergy = kineticEnergy + potentialEnergy

	relativeVariation = np.abs(( np.max(totalEnergy) - np.min(totalEnergy) ) /np.average(totalEnergy))

	print "The total energy has a relative maxvariation of = " , relativeVariation

	tempFigs = py.figure()
	py.plot(time, temperature)
	py.ylabel('Temperature in Kelvin')
	py.xlabel('Time in picoseconds')
	py.title('Temperature of argon gas')
	py.savefig('../data/Plots/temperature.png')

	# energyFig = py.figure()

	# py.subplot(211)
	# py.plot(time ,kineticEnergy)
	# py.ylabel('Energy in Joule')
	# py.xlabel('Time in picoseconds')
	# py.title('Kinetic energy in argon gas')

	# py.subplot(212)
	# py.plot(time , potentialEnergy)
	# py.ylabel('Energy in Joule')
	# py.xlabel('Time in picoseconds')
	# py.title('Potential energy in argon gas')
	# py.savefig('../data/Plots/energyFig.png')

	# totalEnergyFig = py.figure()
	# py.plot(time , totalEnergy)
	# py.ylabel('Energy in Joule')
	# py.xlabel('Time in picoseconds')
	# py.title('Total energy')
	# py.savefig('../data/Plots/totalEnergyFig.png')

	# pressurePlot = py.figure()
	# py.plot(time , pressure)
	# py.ylabel('Pressure in Pa')
	# py.xlabel('Time in picoseconds')
	# py.title('Pressure in the system')
	# py.savefig('../data/Plots/pressureFig.png')
	py.show()

def plotForce():
	## Plotting the force so we can estimate when it approaches zero.

	sigma = 1.0
	epsilon = 1.0

	def lennardJonesForce(x):
		return 24/sigma*epsilon*(2*(sigma/x)**13 - (sigma/x)**7)

	xx = np.linspace(0.01,3,1000)


	forceFig = py.figure()
	py.plot( xx, lennardJonesForce(xx))
	py.ylim( [-3,10] )
	py.grid(True)
	py.title('The force between to particles. sigma = 1, epsilon = 1')
	py.xlabel('Distance')
	py.ylabel('Force')

	py.savefig('../data/Plots/forcePlot.png')

	return

def plotVelocities(time):


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
	velocityDistributions = py.figure()
	xaxis = [-600, 600]
	yaxis = [0,0.003]

	# x = py.subplot(311)
	py.title('Velocity at time %i' %time)

	py.plot(velocityx,fitx,'-')
	py.hist(velocityx,normed=True, bins=100)      #use this to draw histogram of your data
	py.xlim(xaxis)
	py.ylim(yaxis)	


	# y = py.subplot(312)

	# py.plot(velocityy,fity,'-')
	# py.hist(velocityy,normed=True, bins=100)      #use this to draw histogram of your data
	# py.xlim(xaxis)
	# py.ylim(yaxis)	

	# z = py.subplot(313)

	# py.plot(velocityz,fitz,'-')
	# py.hist(velocityz,normed=True, bins=100)      #use this to draw histogram of your data
	# py.xlim(xaxis)
	# py.ylim(yaxis)
	
	py.savefig('../data/Plots/velocityDistribution_' + str(time) + '.eps')
	py.close()
	# py.show()                   #use may also need add this 


	return

def relEnergyError():
	#Just some conversion constants
	timeToPicoSeconds = 1.00224e-13 / 10e-12
	temperatureToKelvin = 119.735
	energyToJoule = 1.65313e-21
	pressureToPa = 1.65313e+09


	statisticalValues = np.genfromtxt('../data/statisticalResults/statisticalValues.csv', delimiter='\t')
	time = timeToPicoSeconds * statisticalValues[1: , 0]
	kineticEnergy = energyToJoule * statisticalValues[1: , 2]
	potentialEnergy =energyToJoule * statisticalValues[1: , 3]

	totalEnergy = kineticEnergy + potentialEnergy

	relativeVariation = np.abs(( np.max(totalEnergy) - np.min(totalEnergy) ) /np.average(totalEnergy))

	print "The total energy has a relative maxvariation of = " , relativeVariation


	return relativeVariation

def plotEnergyVsTimestep(data):

	timestep = data[:,0]
	relativeError = data[:,1]

	energyError = py.figure()
	py.plot(timestep , relativeError)
	py.ylabel('Relative Variation')
	py.xlabel('Timestep in picoseconds')
	py.title('Total energy variation')
	py.savefig('../data/Plots/totalEnergyError.png')

	py.show()


plot_function()
# plotForce()
# plotVelocities()



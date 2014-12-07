import numpy as np
import pylab as py

def plot_function():

	# Plotting of the solutions to the first 
	statisticalValues = np.genfromtxt('statisticalResults/statisticalValues.csv', delimiter='\t')

	timeToPicoSeconds = 1.00224e-13 / 10e-12
	temperatureToKelvin = 119.735
	energyToJoule = 1.65313e-21
	pressureToPa = 1.65313e+09

	time = timeToPicoSeconds * statisticalValues[1: , 0]
	temperature = temperatureToKelvin * statisticalValues[1: , 1]
	kineticEnergy = energyToJoule * statisticalValues[1: , 2]
	potentialEnergy = energyToJoule * statisticalValues[1: , 3]
	pressure = statisticalValues[1: , 5]

	totalEnergy = kineticEnergy + potentialEnergy

	relativeVariation = np.abs(( np.max(totalEnergy) - np.min(totalEnergy) ) /np.average(totalEnergy))

	print "The total energy has a relative maxvariation of = " , relativeVariation

	tempFigs = py.figure()
	py.plot(time, temperature)
	py.ylabel('Temperature in Kelvin')
	py.xlabel('Time in picoseconds')
	py.title('Temperature of argon gas')
	py.savefig('./Plots/temperature.png')

	energyFig = py.figure()

	py.subplot(211)
	py.plot(time ,kineticEnergy)
	py.ylabel('Energy in Joule')
	py.xlabel('Time in picoseconds')
	py.title('Kinetic energy in argon gas')

	py.subplot(212)
	py.plot(time , potentialEnergy)
	py.ylabel('Energy in Joule')
	py.xlabel('Time in picoseconds')
	py.title('Potential energy in argon gas')
	py.savefig('./Plots/energyFig.png')

	totalEnergyFig = py.figure()
	py.plot(time , totalEnergy)
	py.ylabel('Energy in Joule')
	py.xlabel('Time in picoseconds')
	py.title('Total energy')
	py.savefig('./Plots/totalEnergyFig.png')

	pressurePlot = py.figure()
	py.plot(time , pressure)
	py.ylabel('Pressure in Pa')
	py.xlabel('Time in picoseconds')
	py.title('Pressure in the system')
	py.savefig('./Plots/pressureFig.png')
	# py.show()

## Plotting the force so we can estiamte when it approaches zero.

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

py.savefig('./Plots/forcePlot.png')



# legend(loc='upper right')

# N = len(Solutions)

# legend(loc='upper right')
# xlabel('x')
# ylabel('y')
# #ylim(0, 1)
# title('Plot with ' + str(N) + ' steps')
# grid(True)
# savefig("Plot_N_" + str(N) + ".png")
# show()

# # Plotting the relative error
# ErrorTable = genfromtxt('ErrorTable.csv', delimiter=',')

# FigRelErr = figure()
# plot(log10(1/(ErrorTable[:,0]+1)), ErrorTable[:,1], label = 'Relative Error')

# title('Plot of the max relative error dependant on the stepsize')
# legend(loc='upper right')
# xlabel('log10(1/(N+1))')
# ylabel('Max Relative Error')
# grid(True)
# savefig("RelErrorPlot.png")

# show()

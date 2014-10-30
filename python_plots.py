import numpy as np
import pylab as py

# Plotting of the solutions to the first 
statisticalValues = np.genfromtxt('statisticalResults/statisticalValues.tsv', delimiter='\t')

time = statisticalValues[1: , 0]
temperature = statisticalValues[1: , 1]
kineticEnergy = statisticalValues[1: , 2]
potentialEnergy = statisticalValues[1: , 3]

totalEnergy = kineticEnergy + potentialEnergy

relativeVariation = np.abs(( np.max(totalEnergy) - np.min(totalEnergy) ) /np.average(totalEnergy))
relativePotentialVariation = np.abs(( np.max(potentialEnergy) - np.min(potentialEnergy) ) /np.average(totalEnergy))
relativeKineticVariation = np.abs(( np.max(kineticEnergy) - np.min(kineticEnergy) ) /np.average(totalEnergy))

print "The following is the max deviation of the energies relative to the total average energy"
print "Kinetic energy = ", relativeKineticVariation
print "Potential energy = " , relativePotentialVariation
print "Total energy = " , relativeVariation
print "The total energy deviates compared to the single energies deviation is in the order " , \
			 relativeVariation / np.max([relativeKineticVariation, relativePotentialVariation])

tempFigs = py.figure()
py.plot(time, temperature)
py.title('Temperature of argon gas')
py.savefig('./Report/temperature.png')

energyFig = py.figure()

py.subplot(211)
py.plot(time ,kineticEnergy)
py.title('Kinetic energy in argon gas')


py.subplot(212)
py.plot(time , potentialEnergy)
py.title('Potential energy in argon gas')
py.savefig('./Report/energyFig.png')

totalEnergyFig = py.figure()
py.plot(time , totalEnergy)
py.title('Total energy')
py.savefig('./Report/totalEnergyFig.png')

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

py.savefig('./Report/forcePlot.png')



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

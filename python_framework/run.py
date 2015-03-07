from mdconfig import *
from python_plots import *

program = MD(dt=10, name= "../build-molecular-dynamics-Desktop-Debug/molecular-dynamics" )#"molecular-dynamics-fys3150") # dt in femtoseconds, name is the name of the executable you have

###########################################################################
#Let us first create a new system and apply the thermostat, let it thermalise
#and then save it a ready state we can do statistics on. Let us do it for a variety of latticeconstants
#and temperatures
###########################################################################



def testRun():
	# # # Create new FCC lattice
	program.runNewSystem(numberOfUnitCells = 12, FCCLatticeConstant = 5.26) # lattice constant in angstroms

	# # Enable the thermostat for a while to reach a desired temperature. You need to figure out how long it needs to run.
	# program.runThermostat(temperature = 200, timesteps = 1000)
	# program.saveState(path="states/statesTermostat-" + str(latticeConstant))
	# plot_function()


	# # program.loadSavedState(path="states/statesTermostat-" + str(latticeConstant))
	# # # Thermalize for some time (Note: thermalizing means without the thermostat, equilibrate). Here as well you need to figure out how long it needs to thermalize before we start sampling
	program.runThermalize(timesteps = 1000)
	plot_function()
	# # program.saveState(path="../data/states/statesThermalized-lattice" + str(latticeConstant) + "-temperature-" + str(temperature))

	return

def runThermalizeStates():
	for temperature in range (100, 1001,100):

		for latticeConstant in range(3,10):

			# Create new FCC lattice
			program.runNewSystem(numberOfUnitCells = 10, FCCLatticeConstant = latticeConstant) # lattice constant in angstroms

			# Enable the thermostat for a while to reach a desired temperature. You need to figure out how long it needs to run.
			program.runThermostat(temperature = temperature, timesteps = 1000)
			# program.saveState(path="states/statesTermostat-" + str(latticeConstant))

			# program.loadSavedState(path="states/statesTermostat-" + str(latticeConstant))
			# # Thermalize for some time (Note: thermalizing means without the thermostat, equilibrate). Here as well you need to figure out how long it needs to thermalize before we start sampling
			program.runThermalize(timesteps = 1000)
			plot_function()
			# program.saveState(path="../data/states/statesThermalized-lattice" + str(latticeConstant) + "-temperature-" + str(temperature))


	return

def uniformVelocityDistrution():

	#Sets up a new system with uniform velocity distribution
	program.runNewSystem(numberOfUnitCells = 15, gaussianVelocities = False, FCCLatticeConstant = 5.260)
	# program.runThermostat(temperature = 300, timesteps = 1000)

	plotVelocities(0)
	timestep = 10

	for i in range(0,20):
		time = (i+1)*timestep
		program.runThermalize(timesteps = timestep)
		plotVelocities(time)


	return

def energyError():

	errors = np.zeros((20,2))

	#Set up a new system in and let it thermalizem then save the state
	program.runNewSystem(numberOfUnitCells = 10, gaussianVelocities = True, FCCLatticeConstant = 5.260)
	program.runThermostat(temperature = 300, timesteps = 1000)
	program.runThermalize(timesteps = 5000)

	program.saveState(path="../data/states/energyErrorState")

	dt = 0
	stepsToRun = 0

	for i in range(0,10):
		dt +=1
		program.dt = dt

		stepsToRun = 10000/dt


		program.loadSavedState(path="../data/states/energyErrorState")
		program.runThermalize(timesteps = stepsToRun)
		errors[i,0] = dt
		errors[i,1] = relEnergyError()

	# print errors 

	plotEnergyVsTimestep(errors)
	
	return

# runThermalizeStates()
# uniformVelocityDistrution()
# energyError()
testRun()




# Now that it is thermalized, we can trust the states enough to do statistical sampling.
#program.runThermalize(timesteps = 2000)

# plot_function()
# program.saveState(path="states/testState")

# In your main function in your c++ program you can read in these variables as arguments like this: 
# First, we assume some default values, then if arguments are given, use them instead

# int numTimeSteps = 1000;
# double dt = UnitConverter::timeFromSI(1e-14); 
# int numUnitCells = 10;
# float latticeConstant = 5.26;
# bool loadState = false;
# bool thermostatEnabled = false;
# float temperature = 300;
# if(args>1) {
#     dt = UnitConverter::timeFromSI(atof(argv[1])*1e-15);
#     numTimeSteps = atoi(argv[2]);
#     numberOfUnitCells = atoi(argv[3]);
#     latticeConstant = atof(argv[4]);
#     loadState = atoi(argv[5]);
#     thermostatEnabled = atoi(argv[6]);
#     temperature = atof(argv[7]);
# }
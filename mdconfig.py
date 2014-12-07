import subprocess
import os as os
from datetime import datetime

class MD:
	def __init__(self, dt=10, name="molecular-dynamics-fys3150"):
		self.name = name
		self.dt = dt
		self.numberOfTimesteps = 1000
		self.numberOfUnitCells = 10
		self.FCCLatticeConstant = 5.26
		self.loadState = True # If false, the program should create a new state with FCC
		self.thermostatEnabled = False
		self.temperature = 1000

	def run_command(self, cmd):
		subprocess.call(cmd, shell=True)

	def run(self):
		if not os.path.isfile(self.name):
			print "Executable "+self.name+" does not exist, aborting!"
			exit()
		
		command = "./%s %f %d %d %f %d %d %f" % (self.name, self.dt, self.numberOfTimesteps, self.numberOfUnitCells, 
													self.FCCLatticeConstant, self.loadState, self.thermostatEnabled, 
													self.temperature)
		print "Running command: ", command

		now = datetime.now()
		self.run_command(command)
		t1 = (datetime.now() - now).seconds
		timeStepsPerSecond = self.numberOfTimesteps / max(t1,1)

		print "Process used %d seconds (%d timesteps per second)" % ( t1, timeStepsPerSecond )

	def loadSavedState(self, path):
		print "Got into here?"
		print "Loading state from "+str(path)
		self.run_command("cp -r "+path+"/state ./")

	def saveState(self, path):
		print "Saving state to "+str(path)
		self.run_command("mkdir -p "+path)
		self.run_command("cp -r statisticalResults/statisticalValues.csv "+path)
		self.run_command("cp -r Plots/* "+path)
		self.run_command("cp -r state "+path)

	def runNewSystem(self, numberOfUnitCells, FCCLatticeConstant = 5.26):
		self.loadState = False
		self.numberOfTimesteps = 1
		self.numberOfUnitCells = numberOfUnitCells
		self.FCCLatticeConstant = FCCLatticeConstant
		self.run()
		self.loadState = True

	def runThermostat(self, temperature, timesteps):
		self.temperature = temperature
		self.numberOfTimesteps   = timesteps
		self.thermostatEnabled = True
		self.run()
		self.thermostatEnabled = False

	def runThermalize(self, timesteps):
		self.thermostatEnabled = False
		self.numberOfTimesteps   = timesteps
		self.run()
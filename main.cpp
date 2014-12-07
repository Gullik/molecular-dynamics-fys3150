#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/eulercromer.h>
#include <integrators/velocityverlet.h>
#include <lists/neighborlist.h>
#include <thermostats/berendsen.h>

#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>
#include <time.h>

#include <sstream>
#include <stdlib.h>

using namespace std;



int main(int argc, char* argv[])
{

    //Setting a lot of variables and initiliazing different parts of the system and objects in it

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;
    cout << "One unit of energy is " << UnitConverter::energyToSI(1.0) << " J" << endl;

    int numTimeSteps = 2000;
    double dt = UnitConverter::timeFromSI(1e-15);
    int numUnitCells = 10;
    float latticeConstant = 5.26;
    bool loadState = false;
    bool thermostatEnabled = false;
    float temperature = 300;
    double sigma = UnitConverter::lengthFromAngstroms(3.405);
    double epsilon = UnitConverter::energyFromSI(119.8*UnitConverter::kb);
    double cutOffLength = 2*sigma;

    if(argc>1) {
        dt = UnitConverter::timeFromSI(atof(argv[1])*1e-16);
        numTimeSteps = atoi(argv[2]);
        numUnitCells = atoi(argv[3]);
        latticeConstant = atof(argv[4]);
        loadState = atoi(argv[5]);
        thermostatEnabled = atoi(argv[6]);
        temperature = atof(argv[7]);
    }


    //Let python and the framework handle the statefile instead, just store it as the statefile and saves and python
    //copies the file if it is needed
    string statefile = "state";


    //Initializing
    System system;
    //Creating a new system of atoms, or loading up a ready made state
    if(!loadState)
    {
        system.createFCCLattice(numUnitCells, latticeConstant);
        system.removeMomentum();
        system.setPotential(new LennardJones(sigma, epsilon));
        system.setIntegrator(new VelocityVerlet());
        system.setNeighborList(new NeighborList(cutOffLength, system.systemSize()));
        system.setThermostat(new Berendsen(UnitConverter::temperatureFromSI(temperature), 0.1));
    }
    else
    {
        cout << "Loaded statefile " << statefile << endl;
        system.load(statefile, &system);
        system.setPotential(new LennardJones(sigma, epsilon));
        system.setIntegrator(new VelocityVerlet());
        system.setNeighborList(new NeighborList(cutOffLength, system.systemSize()));
        system.setThermostat(new Berendsen(UnitConverter::temperatureFromSI(temperature), 0.1));
    }

    if(thermostatEnabled)
    {
        cout << "Turned on thermostat" << endl;
        system.thermostat()->turn_on();
    }


    vec3 systemSize = system.systemSize();

    cout << "The systemSize is: " << systemSize << endl;
    cout << "The cutOffLength is: " << system.list()->cutOffLength() << endl;
    cout << "No of neighbors in each direction is: " << system.list()->noOfNeighbors() << endl;

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    clock_t start, end;     //To keep track of the time

    start = clock();

    IO *movie = new IO();   // To write the state to file
    movie->open("../molecular-dynamics-fys3150/movie.xyz");

    movie->saveState(&system); //The initial state should also be recorded
//    statisticsSampler->sample(&system);

    for(int timestep=0; timestep < numTimeSteps ; timestep++) {
        system.step(dt);

        if( timestep % 10 == 0)
            statisticsSampler->sample(&system);
        movie->saveState(&system);

        if(timestep % 100 == 0)
            cout << "At timestep: " << timestep << endl;
    }

    movie->close();

    end = clock();

    double time = 1.0*(end - start)/CLOCKS_PER_SEC;

   cout << "The program spent " << time << " s to calculate the steps, with " << system.atoms().size() << " atoms" <<  endl;

   //Saving the state
   system.save(statefile, &system);
   cout << "Saved state as " << statefile << endl;



    return 0;
}


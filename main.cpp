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

using namespace std;

int main()
{
    double dt = UnitConverter::timeFromSI(1e-15); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;
    cout << "One unit of energy is " << UnitConverter::energyToSI(1.0) << " J" << endl;

    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26);  //Not really necessary since the programs units already uses Angrstoms, but
    double sigma = UnitConverter::lengthFromAngstroms(3.405);           //safer to keep in case it changes base units
    double epsilon = UnitConverter::energyFromSI(119.8*UnitConverter::kb);
    int gridNodes = 6;      //After the neighborlists implementation it needs enough gridnodes to have 3+ neighors to work properly
    double cutOffLength = 3*sigma;


    System system;
    system.createFCCLattice(gridNodes, latticeConstant);
    system.setPotential(new LennardJones(sigma, epsilon));
    system.setIntegrator(new VelocityVerlet());
    system.setNeighborList(new NeighborList(cutOffLength, system.systemSize()));
    system.setThermostat(new Berendsen(UnitConverter::temperatureFromSI(200), 0.1));

    system.removeMomentum();

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



    for(int timestep=0; timestep < 1000 ; timestep++) {
        system.step(dt);

        if(timestep == 300)
            system.thermostat()->turn_off();


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

    return 0;
}


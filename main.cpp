#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/eulercromer.h>
#include <integrators/velocityverlet.h>

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
    int gridNodes = 1;


    System system;
    system.createFCCLattice(gridNodes, latticeConstant);
    system.setPotential(new LennardJones(sigma, epsilon));
    system.setIntegrator(new VelocityVerlet());

    system.removeMomentum();





    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    clock_t start, end;     //To keep track of the time

    start = clock();

    IO *movie = new IO(); // To write the state to file
    movie->open("../molecular-dynamics-fys3150/movie.xyz");

    movie->saveState(&system); //The initial state should also be recorded
//    statisticsSampler->sample(&system);



    for(int timestep=0; timestep<1; timestep++) {
        system.step(dt);
        if( timestep % 10 == 0)
            statisticsSampler->sample(&system); //Kinetic part still not working properly
        movie->saveState(&system);

    }

    movie->close();

    end = clock();

    double time = 1.0*(end - start)/CLOCKS_PER_SEC;

   cout << "The program spent " << time << " s to calculate the steps, with " << system.atoms().size() << " atoms" <<  endl;


    return 0;
}

//        test = 0;

//        for(int i = 0 ; i < system.atoms().size(); i++)

//            if( system.atoms()[i]->position.x < system.systemSize().x   &&
//                system.atoms()[i]->position.y < system.systemSize().y   &&
//                system.atoms()[i]->position.z < system.systemSize().z   &&
//                system.atoms()[i]->position.x >= 0                      &&
//                system.atoms()[i]->position.y >= 0                      &&
//                system.atoms()[i]->position.z >= 0
//                    )
//                test += 1;

//        cout << test << endl;


//    for(int n=0; n<100; n++) {
//        // Add one example atom. You'll have to create many such atoms in the createFCCLattice function above.
//        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // Argon mass, see http://en.wikipedia.org/wiki/Argon
////        atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
////        atom->position.randomUniform(0, system.systemSize().x);
//        system.atoms().push_back(atom); // Add it to the list of atoms
//    }



//    cout << "Calculating force " << endl;

//    system.calculateForces();

//    return 0;

//    double test = 0;

//    for(int i = 0 ; i < system.atoms().size(); i++)

//        if( system.atoms()[i]->position.x < system.systemSize().x   &&
//            system.atoms()[i]->position.y < system.systemSize().y   &&
//            system.atoms()[i]->position.z < system.systemSize().z   &&
//            system.atoms()[i]->position.x >= 0                      &&
//            system.atoms()[i]->position.y >= 0                      &&
//            system.atoms()[i]->position.z >= 0
//                )
//            test += 1;

//    cout << test << endl;

//    vec3 size = system.systemSize();

//    cout << size << endl;

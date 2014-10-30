#include <statisticssampler.h>
#include <potentials/potential.h>
#include <fstream>
#include <unitconverter.h>

ofstream myfile;

StatisticsSampler::StatisticsSampler()
{   
    //Clear the old statistical file and prepares a new one
        ofstream("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.tsv",ios::trunc);

        fstream myfile;
        myfile.open("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.tsv");
        myfile << "Time (?)"  << "\t" "Temperature" << "\t" << "Kinetic energy (?)" << "\t" << "Potential Energy (?)" << endl;
        myfile.close();
}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);                //Also calculates the instanteneous temperature
    samplePotentialEnergy(system);
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    // This takes all the atoms in the system and calculates the Kinetic energy according to E_k = 1/2mv², then they are written to a file

    double kineticEnergy = 0;
    double temperature = 0;

    for(int i = 0 ; i < system->atoms().size() ; i++)
    {
        Atom *atom = system->atoms()[i];

        kineticEnergy += atom->mass()*atom->velocity.lengthSquared() / 2.0 ;
    }

    temperature = 2.0/3.0*kineticEnergy / system->atoms().size();

    myfile.open ("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.tsv", ios::app);
             myfile << endl <<system->currentTime() << "\t" <<  UnitConverter::temperatureToSI( temperature )  <<"\t" << kineticEnergy << "\t";
    myfile.close();


}

void StatisticsSampler::samplePotentialEnergy(System *system)
{
    myfile.open ("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.tsv", ios::app);
             myfile << system->potential()->potentialEnergy();
    myfile.close();
}


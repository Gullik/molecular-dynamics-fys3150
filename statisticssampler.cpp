#include <statisticssampler.h>
#include <potentials/potential.h>
#include <fstream>
#include <unitconverter.h>

ofstream myfile;

StatisticsSampler::StatisticsSampler()
{   
    //Clear the old statistical file and prepares a new one
        ofstream("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.csv",ios::trunc);

        fstream myfile;
        myfile.open("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.csv");
        myfile << "Time (?)"  << "\t" << "Temperature" << "\t" << "Kinetic energy (?)" << "\t" << "Potential Energy (?)" <<
                  "\t" << "Number Density" << "\t" << "Pressure" << endl;
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
    sampleNumberDensity(system);
    samplePressure(system);

    //Write the measured values to file

    myfile.open ("../molecular-dynamics-fys3150/statisticalResults/statisticalValues.csv", ios::app);
    myfile << system->currentTime() << "\t" << m_temperature << "\t" << m_kineticEnergy << "\t" << m_potentialEnergy <<
              "\t" << m_numberDensity << "\t" << m_pressure <<endl;
    myfile.close();


}

void StatisticsSampler::sampleKineticEnergy(System *system)
{
    // This takes all the atoms in the system and calculates the Kinetic energy according to E_k = 1/2mv², then they are written to a file

    double kineticEnergy = 0;
    double temperature = 0;

    for(int i = 0 ; i < int(system->atoms().size()) ; i++)
    {
        Atom *atom = system->atoms()[i];

        kineticEnergy += atom->mass()*atom->velocity.lengthSquared() / 2.0 ;
    }

    temperature = 2.0/3.0*kineticEnergy / system->atoms().size();

    //Stores it in the statisticsSampler
    m_kineticEnergy = kineticEnergy;
    m_temperature = temperature;


}

void StatisticsSampler::samplePotentialEnergy(System *system)
{


    m_potentialEnergy = system->potential()->potentialEnergy();
}

void StatisticsSampler::sampleNumberDensity(System *system)
{
    double volume = system->systemSize().length();
    int noOfAtoms = system->atoms().size();

    double numberDensity = 1.0 * noOfAtoms / volume;


    m_numberDensity = numberDensity;

}

void StatisticsSampler::samplePressure(System *system)
{

    //This adds together the pressure resulting from the temperature and the pressure calculated from
    //the force felt between particles calculated in the force loop in the potential.

    double pressure = m_numberDensity * m_temperature + system->potential()->pressure() ;

    m_pressure = pressure;
}

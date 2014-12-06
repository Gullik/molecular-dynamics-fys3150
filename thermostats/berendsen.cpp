#include "berendsen.h"
#include <math.h>


Berendsen::Berendsen(double temperature, double tau):
m_tau(tau)

{
    m_wantedTemperature = temperature;
}

void Berendsen::applyThermostat(System *system, double dt)
{
    //First it calculates the temperature
    double kineticEnergy = 0;
    double temperature = 0;
    double gamma = 0;

    for(int i = 0 ; i < int(system->atoms().size()) ; i++)
    {
        Atom *atom = system->atoms()[i];

        kineticEnergy += atom->mass()*atom->velocity.lengthSquared() / 2.0 ;
    }

    temperature = 2.0/3.0*kineticEnergy / system->atoms().size();

    //Calculates the factor all velocities should be multiplied with, gamma = sqrt( 1 + dt/tau (T_bath/T_wanted - 1)  )
    gamma = sqrt( 1 + dt/m_tau*(m_wantedTemperature/temperature - 1));

    //Adjust all the velocities for all the atoms
    for(int i = 0 ; i < int(system->atoms().size()); i++)
    {
        system->atoms()[i]->velocity = system->atoms()[i]->velocity * gamma;
    }

//    cout << temperature << "vs" << m_wantedTemperature <<endl;
//    cout << gamma << endl;

}



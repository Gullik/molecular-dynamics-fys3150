#ifndef BERENDSEN_H
#define BERENDSEN_H

#include <thermostats/thermostat.h>

class Berendsen : public Thermostat
{
private:
    double m_tau;           //Relaxation time

public:
    Berendsen(double temperature, double tau);
    virtual void applyThermostat(System *system, double dt);

};

#endif // BERENDSEN_H

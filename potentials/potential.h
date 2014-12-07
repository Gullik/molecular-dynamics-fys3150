#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    double m_potentialEnergy;
    double m_pressure;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    double potentialEnergy();
    void setPotentialEnergy(double potentialEnergy);

    double pressure() {return m_pressure;}
    void setPressure(double pressure) {m_pressure = pressure;}
};

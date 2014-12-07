#pragma once
#include <system.h>

class StatisticsSampler
{
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System *system);
    void sampleKineticEnergy(System *system);
    void samplePotentialEnergy(System *system);
    void sampleNumberDensity(System *system);
    void samplePressure(System *system);

    //Sampled properties
    double m_temperature;
    double m_numberDensity;
    double m_potentialEnergy;
    double m_kineticEnergy;
    double m_pressure;
};

#pragma once
#include <potentials/potential.h>
#include <math/vec3.h>

using CompPhys::vec3;

class LennardJones : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
public:
    LennardJones(double sigma, double epsilon);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
    void calculateForcesBetweenCells(vec3 neighborPositionA, vec3 neighborPositionB, System *system);
    void calculateForcesInsideCell(vec3 neighborPosition, System *system);
};

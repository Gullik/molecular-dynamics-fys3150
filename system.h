#pragma once
#include <vector>
#include <atom.h>
#include <math/vec3.h>
#include <iostream>

using namespace std;
class Potential; class Integrator; class NeighborList;
using std::vector;
using CompPhys::vec3;

class System
{
private:
    vec3 m_systemSize;
    vector<Atom*> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_currentTime;
    int m_steps;
    NeighborList *m_list;

public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces();
    void step(double dt);

    void setNeighborList(NeighborList *neighborlist) { m_list = neighborlist; }
    NeighborList *list() {return m_list;}

//    NeighborList *neighborList() {return m_neighborList;}

    // Setters and getters
    vector<Atom *> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }



};

#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <iostream>

using namespace std;

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention

    double xLength =  m_systemSize.x;     // Vectors that should shift the atom one systemsize to the side when it
    double yLength =  m_systemSize.y;    // when it goes over it.
    double zLength =  m_systemSize.z;

    int k = 0;

    for(int i = 0; i < m_atoms.size(); i++)
    {
        //Takes all the for atoms and shifts them to the other side if
        //if they overstep the boundary

        if(m_atoms[i]->position.x > xLength)
            m_atoms[i]->position.x -= xLength;

        if(m_atoms[i]->position.x < 0)
            m_atoms[i]->position.x += xLength;

        if(m_atoms[i]->position.y > yLength)
            m_atoms[i]->position.y -= yLength;

        if(m_atoms[i]->position.y < 0)
            m_atoms[i]->position.y += yLength;

        if(m_atoms[i]->position.z > zLength)
            m_atoms[i]->position.z -= zLength;

        if(m_atoms[i]->position.z < 0)
            m_atoms[i]->position.z += zLength;
    };

    return;
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
}

void System::resetForcesOnAllAtoms() {

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {

}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}

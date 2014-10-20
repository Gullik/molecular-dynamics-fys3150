#include <integrators/velocityverlet.h>
#include <iostream>
#include <system.h>

using namespace std;

VelocityVerlet::VelocityVerlet() :
    m_firstStep(false) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::firstKick(System *system, double dt)
{
    cout << "Firstkick got called" << endl;

}

void VelocityVerlet::halfKick(System *system, double dt)
{
 //   cout << "Halfkick got called" << endl;
    //Halfkick shall calculate the velocity half a step forward
    for(int i = 0; i < system->atoms().size(); i++)
    {
        Atom *atom = system->atoms()[i];
        double halfTimeStepByMass = dt/(2*atom->mass());
        atom->velocity.addAndMultiply(atom->force, halfTimeStepByMass); // v(t+dt/2)=v(t)+ F(t)/m*dt/2
    }

}

void VelocityVerlet::move(System *system, double dt)
{
 //   cout << "Move got called" << endl;

//    cout << "Before move: " << system->atoms()[0]->position << endl;
    for(int i = 0; i < system->atoms().size(); i++)
    {
        Atom *atom = system->atoms()[i];
        atom->position.addAndMultiply( atom->velocity, dt) ; // r(t+ dt)= r(t) + v(t+dt/2)*dt
    }
//    cout << "After move: " << system->atoms()[0]->position << endl;

}

void VelocityVerlet::integrate(System *system, double dt)
{
    system->calculateForces();

    if(m_firstStep) {           // Do not know the purpose of first kick thought?
        firstKick(system, dt);
        m_firstStep = false;
    }
    halfKick(system, dt);   //This implements the first step, v(t + dt/2)= v(t) + F(t)/m*dt/2

    move(system, dt);       //System is moved, r(t+ dt)= r(t) + v(t+dt/2)*dt

    system->applyPeriodicBoundaryConditions();  //Put's them inside the boundary again

    system->calculateForces();      //New forces should be calculated, according to F(t+dt) = - d/dr U(r(r+dt))

    halfKick(system, dt);   //Kicks the second half of the step v(t+dt)=v(t+dt/2) + F(t+dt)/m*dt/2

//    cout << system->atoms()[0]->force  << endl;
}

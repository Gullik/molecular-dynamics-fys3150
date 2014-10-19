#include <potentials/lennardjones.h>
#include <iostream>

using namespace std;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    cout << "Got into potential" << endl;
    //This should calculate the force between to particles according to the gradient
    //of the L-J potential U = -4 E ((o/r)^12-(o/r)^6).
    //
    //It first calculates the length between a pair of molecules, then the force between
    //them is calculated as according to the formuala given in the report
    //Then the force is assigned to both atoms

    double lengthBetween = 0;
    vec3 vecBetween = vec3();

    for(int i = 0; i < system->atoms().size(); i++)
    {
        Atom *atom = system->atoms()[i];
        for(int j = i+1; j < system->atoms().size(); j++)
        {
            cout << i << j << endl;
            vecBetween = atom[j].position - atom[i].position;
            cout << vecBetween << endl;
        }
    }

    m_potentialEnergy = 0; // Remember to compute this in the loop
}

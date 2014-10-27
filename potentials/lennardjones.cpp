#include <potentials/lennardjones.h>
#include <iostream>
#include <math.h>

using namespace std;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
//    cout << "Got into potential" << endl;
    //This should calculate the force between to particles according to the gradient
    //of the L-J potential U = -4 E ((o/r)^12-(o/r)^6).
    //
    //It first calculates the length between a pair of molecules, then the force between
    //them is calculated as according to the formuala given in the report
    //Then the force is assigned to both atoms

    double Lx = system->systemSize().x;
    double Ly = system->systemSize().y;
    double Lz = system->systemSize().z;

    int noOfAtoms = system->atoms().size();

    vector<Atom*> atoms = system->atoms();

    double lengthBetween = 0;
    vec3 vecBetween = vec3();
    vec3 force = vec3();
    double forceMagnitude = 0;
    double forceConstant = 24*m_epsilon/m_sigma;
    double helpTerm = 0;
    

    for(int i = 0; i < noOfAtoms; i++)
    {
        for(int j = i+1; j < noOfAtoms; j++)
        {
            //Because of the periodic boundary conditions, in which an atom is shifted to the other side if it reaches the boundary
            //we also want the same to be true for forces. So the distance between two particles is the one that is shortest of the
            //straight line, or through the walls. The line between A and B is the shortest line of the starred line and the dotted line in
            //the drawing below.
            //  |                   |
            //  |                   |
            //  |***A - - - - B*****|
            //  |                   |
            //  |                   |
            //  |                   |

            vecBetween = atoms[i]->position - atoms[j]->position;   //The r_ij vector

            if(vecBetween.x > Lx/2)
                vecBetween.x -= Lx;
            if(vecBetween.x < -Lx/2)
                vecBetween.x += Lx;

            if(vecBetween.y > Ly/2)
                vecBetween.y -= Ly;
            if(vecBetween.y < -Ly/2)
                vecBetween.y += Ly;

            if(vecBetween.z > Lz/2)
                vecBetween.z -= Lz;
            if(vecBetween.z < -Lz/2)
                vecBetween.z += Lz;

            lengthBetween = vecBetween.length(); //The length of the distance between the particles, vecBetween is used to give direction to the force.



            //Calculating the force felt from each atom and adding them to the force each atom feels

            forceMagnitude = forceConstant * (2*pow(m_sigma/lengthBetween,13) - pow(m_sigma/lengthBetween,7)) ;

            force = vecBetween/lengthBetween * forceMagnitude; //Normalizing the vector to correspond to k_ij / |r_ij|


            atoms[i]->force.add(force);
            atoms[j]->force.addAndMultiply(force, -1);

            //Calculating the potential energy given by L-J 12-6 since we anyway are in the loop and has distances calculated
            // U = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]

            helpTerm = pow(m_sigma/lengthBetween, 6);

            m_potentialEnergy += 4*m_epsilon*(helpTerm*helpTerm - helpTerm) ; //The potential Energy should be counted only once for each molecule pair?

        }
    }


}

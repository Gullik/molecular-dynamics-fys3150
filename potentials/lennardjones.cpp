#include <potentials/lennardjones.h>
#include <iostream>
#include <math.h>
#include <lists/neighborlist.h>
#include <math/vec3.h>

using namespace std;
using CompPhys::vec3;

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{

    //This should calculate the force between to particles according to the gradient
    //of the L-J potential U = -4 E ((o/r)^12-(o/r)^6).
    //
    //It first calculates the length between a pair of molecules, then the force between
    //them is calculated as according to the formuala given in the report
    //Then the force is assigned to both atoms

    //After some new modifications it should not calculate the force between all the atoms, but only those in bordering neighborhoods
    //so if every atom in a neighborhood calculates the force between them and the next neighborhoods, all the bordering neighborhoods should
    //get calculated.
    //At the endpoints it needs to be modified


    int noOfNeighbors = system->list()->noOfNeighbors();
    int next_i,next_j,next_k;
    int prev_i, prev_j, prev_k;


    //Time to go through the list, and gather all the atoms the force should be calculated between
    for(int i=0; i < noOfNeighbors; i++)
    {
        for(int j=0; j < noOfNeighbors ; j++)
        {
            for(int k=0; k < noOfNeighbors ; k++)
            {

                //Now atoms should be appended with the increasing neighboring cells
                //Could not find any neat way to do it, this is confusing and slow as is. Should work to do it a better way.
                next_i = i + 1;
                next_j = j + 1;
                next_k = k + 1;

                prev_i = i -1;
                prev_j = j -1;
                prev_k = k -1;

                if(i == noOfNeighbors - 1)
                    next_i = 0;

                if(j == noOfNeighbors - 1)
                    next_j = 0;

                if(k == noOfNeighbors - 1)
                    next_k = 0;

                if(i == 0)
                    prev_i = noOfNeighbors - 1;

                if(j == 0)
                    prev_j = noOfNeighbors - 1;

                if(k == 0)
                    prev_k = noOfNeighbors - 1;

                //With itself
                calculateForcesInsideCell(vec3(i,j,k),system);
//                calculateForcesBetweenCells(vec3(i,j,k), vec3(i,   j, k), system);

                //The top layer, so it cycles through i and j while k is constant

                calculateForcesBetweenCells(vec3(i,j,k), vec3(prev_i,   prev_j, next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(i,        prev_j, next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   prev_j, next_k), system);

                calculateForcesBetweenCells(vec3(i,j,k), vec3(prev_i,   j,      next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(i,        j,      next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   j,      next_k), system);

                calculateForcesBetweenCells(vec3(i,j,k), vec3(prev_i,   next_j, next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(i,        next_j, next_k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   next_j, next_k), system);


                //The front layer that is left
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   prev_j,      k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   j,           k), system);
                calculateForcesBetweenCells(vec3(i,j,k), vec3(next_i,   next_j,      k), system);

                //The rest of the left side
                calculateForcesBetweenCells(vec3(i,j,k),  vec3(i,        next_j,      k), system);


            }
        }
    }

}

void LennardJones::calculateForcesInsideCell(vec3 neighborPosition, System *system)
{
    vector<Atom*> atoms = system->list()->neigborList()[neighborPosition.x][neighborPosition.y][neighborPosition.z];

    int noOfAtoms = atoms.size();

    double lengthBetween = 0;
    vec3 vecBetween = vec3();
    vec3 force = vec3();
    double forceMagnitude = 0;
    double forceConstant = 24*m_epsilon/m_sigma;
    double helpTerm = 0;
    double volume = pow(system->systemSize().x,3);

    double cutOffPotentialEnergy = 4*m_epsilon*(pow(m_sigma/system->list()->cutOffLength(),12)
                                                - pow(m_sigma/system->list()->cutOffLength(),6));


    for(int i = 0; i < noOfAtoms; i++)
    {
        for(int j = i + 1; j < noOfAtoms; j++)
        {

        vecBetween = atoms[j]->position - atoms[i]->position;   //The r_ij vector


        lengthBetween = vecBetween.length(); //The length of the distance between the particles, vecBetween is used to give direction to the force.

        if(fabs(lengthBetween) < system->list()->cutOffLength())
        {

            //Calculating the force felt from each atom and adding them to the force each atom feels

            forceMagnitude = forceConstant * (2*pow(m_sigma/lengthBetween,13) - pow(m_sigma/lengthBetween,7)) ;

            force = vecBetween/lengthBetween * forceMagnitude; //Normalizing the vector to correspond to k_ij / |r_ij|

            atoms[i]->force.addAndMultiply(force, -1);
            atoms[j]->force.add(force);

            atoms[i]->m_forceCounter += 1;
            atoms[j]->m_forceCounter += 1;

            //Calculating the potential energy given by L-J 12-6 since we anyway are in the loop and has distances calculated
            // U = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]

            helpTerm = pow(m_sigma/lengthBetween, 6);

            m_potentialEnergy += 4*m_epsilon*(helpTerm*helpTerm - helpTerm) - cutOffPotentialEnergy; //The potential Energy should be counted only once for each molecule pair?

            // This calculates the pressure due to forces between the particles, the constant part of the pressure due to rho_n*k_b*T
            // is added in the statisicssampler function.
            //So this calculates one part one of the summation in sum_(i>j) F_ij*r_ij
            m_pressure += 1.0/(3.0*volume) * forceMagnitude;

        }


        }
    }


}

void LennardJones::calculateForcesBetweenCells(vec3 neighborPositionA, vec3 neighborPositionB, System *system)
{
    //This should only calculate the forces between the different particles sent into it, while calculateForces
    //sort which atoms should be calculated together and calls this when it knows which ones are together

//    cout << "Calculates ForcesBetweenCells" << endl;

    vector<Atom*> a_atoms = system->list()->neigborList()[neighborPositionA.x][neighborPositionA.y][neighborPositionA.z];
    vector<Atom*> b_atoms = system->list()->neigborList()[neighborPositionB.x][neighborPositionB.y][neighborPositionB.z];

    int noOfAtomsInA = a_atoms.size();
    int noOfAtomsInB = b_atoms.size();

    double lengthBetween = 0;
    vec3 vecBetween = vec3();
    vec3 force = vec3();
    double forceMagnitude = 0;
    double forceConstant = 24*m_epsilon/m_sigma;
    double helpTerm = 0;

    double Lx = system->systemSize().x;
    double Ly = system->systemSize().y;
    double Lz = system->systemSize().z;

    double volume = pow(system->systemSize().x,3);

    double cutOffPotentialEnergy = 4*m_epsilon*(pow(m_sigma/system->list()->cutOffLength(),12)
                                                - pow(m_sigma/system->list()->cutOffLength(),6));

    for(int i = 0; i < noOfAtomsInA; i++)
    {
        for(int j = 0; j < noOfAtomsInB; j++)
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

        vecBetween = b_atoms[j]->position - a_atoms[i]->position;   //The r_ij vector

        //Implements the periodic boundary conditions, forces should act over the boundary
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

        if(lengthBetween == 0)
            continue;

        if(fabs(lengthBetween) < system->list()->cutOffLength())
        {

            //Calculating the force felt from each atom and adding them to the force each atom feels

            forceMagnitude = forceConstant * (2*pow(m_sigma/lengthBetween,13) - pow(m_sigma/lengthBetween,7)) ;

            force = vecBetween/lengthBetween * forceMagnitude; //Normalizing the vector to correspond to k_ij / |r_ij|


            a_atoms[i]->force.addAndMultiply(force, -1);
            b_atoms[j]->force.add(force);

            a_atoms[i]->m_forceCounter += 1;
            b_atoms[j]->m_forceCounter += 1;

            //Calculating the potential energy given by L-J 12-6 since we anyway are in the loop and has distances calculated
            // U = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]

            helpTerm = pow(m_sigma/lengthBetween, 6);

            m_potentialEnergy += 4*m_epsilon*(helpTerm*helpTerm - helpTerm) - cutOffPotentialEnergy; //The potential Energy should be counted only once for each molecule pair?

        }


        }
    }



}

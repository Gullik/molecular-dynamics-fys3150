#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <unitconverter.h>
#include <math.h>
#include <lists/neighborlist.h>

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

    double xLength =  m_systemSize.x;       // Lengths that should shift the atom one systemsize to the side when it
    double yLength =  m_systemSize.y;       // goes over it.
    double zLength =  m_systemSize.z;


    for(int i = 0; i < int(m_atoms.size()); i++)
    {
        //Takes all the for atoms and shifts them to the other side if
        //if they overstep the boundary

        if(m_atoms[i]->position.x >= xLength)
            m_atoms[i]->position.x -= xLength;

        if(m_atoms[i]->position.x < 0)
            m_atoms[i]->position.x += xLength;

        if(m_atoms[i]->position.y >= yLength)
            m_atoms[i]->position.y -= yLength;

        if(m_atoms[i]->position.y < 0)
            m_atoms[i]->position.y += yLength;

        if(m_atoms[i]->position.z >= zLength)
            m_atoms[i]->position.z -= zLength;

        if(m_atoms[i]->position.z < 0)
            m_atoms[i]->position.z += zLength;
    };

    return;
}

void System::removeMomentum() {
    // The function adds together the total momentum from all the atoms in the system, then it removes it to avoid drift of the system since the inital condition
    // gaussian velocity distribution has a non-zero momentum. If several types of atoms are madeit need to be multiplied with atom mass which is not done for now since it is
    // unnessary

    vec3 momentum = vec3();
    int noOfAtoms = m_atoms.size();

    for(int i = 0; i < noOfAtoms; i ++)
        momentum = momentum + m_atoms[i]->velocity;

    for(int i = 0; i < noOfAtoms; i ++)
    {
        m_atoms[i]->velocity.x = m_atoms[i]->velocity.x - (momentum.x / noOfAtoms);
        m_atoms[i]->velocity.y = m_atoms[i]->velocity.y - (momentum.y / noOfAtoms);
        m_atoms[i]->velocity.z = m_atoms[i]->velocity.z - (momentum.z / noOfAtoms);
    }

}

void System::resetForcesOnAllAtoms() {

    for(int i = 0; i < int(m_atoms.size()); i++)
        m_atoms[i]->force.setToZero();

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant) {

    // This function creates Nx*Ny*Nz cells and places four atoms, for each cell in a crystalline lattice
    // in those cells

    //Setting the systemsize
    this->setSystemSize(vec3(numberOfUnitCellsEachDimension *latticeConstant ,
                             numberOfUnitCellsEachDimension *latticeConstant ,
                             numberOfUnitCellsEachDimension *latticeConstant ));



    //Initialize the atoms and store pointers to them in the system->atoms list.
    //We need four atoms for each cell
    for(int i = 0; i < 4*pow(numberOfUnitCellsEachDimension, 3); i++)
    {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atoms().push_back(atom);
    }


    cout << "The density of the system is " << (this->atoms().size()* UnitConverter::massFromSI(6.63352088e-26) )
            / this->systemSize().lengthSquared() << " in atomic units "<< endl;

    int latticeNumber = 0;

    vec3 r_1 = vec3(0,0,0);
    vec3 r_2 = vec3(latticeConstant/2 , latticeConstant/2.0 , 0 );
    vec3 r_3 = vec3(0 , latticeConstant/2.0 , latticeConstant/2 );
    vec3 r_4 = vec3(latticeConstant/2 , 0 , latticeConstant/2.0 );


    for(int i = 0; i < numberOfUnitCellsEachDimension ; i++)
    {
        for(int j = 0; j < numberOfUnitCellsEachDimension ; j++)
        {
            for(int k = 0; k < numberOfUnitCellsEachDimension ; k++)
            {
                //Four atoms need to be put in each lattice cell, they have the position:R_ij = R_i + r_j;
                //where R_ij are their position, R_i is the lattice position and r_j is the position in the lattice

                m_atoms[latticeNumber    ]->position = vec3(i,j,k)*latticeConstant + r_1;
                m_atoms[latticeNumber + 1]->position = vec3(i,j,k)*latticeConstant + r_2;
                m_atoms[latticeNumber + 2]->position = vec3(i,j,k)*latticeConstant + r_3;
                m_atoms[latticeNumber + 3]->position = vec3(i,j,k)*latticeConstant + r_4;

                latticeNumber += 4;

            }
        }
    }

}

void System::calculateForces() {
    resetForcesOnAllAtoms();



    m_list->sortAtoms(this);



    m_potential->setPotentialEnergy(0);
    m_potential->calculateForces(this);

}

void System::step(double dt) {

    m_integrator->integrate(this, dt);

    m_steps++;
    m_currentTime += dt;
}

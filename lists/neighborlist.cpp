#include "neighborlist.h"
#include "atom.h"
#include <system.h>

NeighborList::NeighborList(double cutOffLength, vec3 systemSize)
{
    //This should only contain atoms? and it needs a function to add and delete atoms as they leave the box territory,
    //or enter it
    //It also need an identification vector, like neighbor (i,j,k)

    //Setting up lists

    int noOfNeighbors = int(systemSize.x/cutOffLength);

    vector<Atom*> atomList;

    vector< vector<Atom*> >  rowList;

    vector< vector< vector<Atom*> > > gridList;

    vector< vector< vector< vector<Atom*> > > > neighborList;

    for(int i = 0; i < noOfNeighbors ; i++)
    {
        rowList.push_back(atomList);
    }

    for(int j = 0; j < noOfNeighbors; j++ )
    {
        gridList.push_back(rowList);
    }

    for(int j = 0; j < noOfNeighbors; j++ )
    {
        neighborList.push_back(gridList);
    }

    m_neighborList = neighborList;

    m_noOfNeighbors = noOfNeighbors;          //Storing it in the class for convenience

    m_cutOffLength = cutOffLength;

}

void NeighborList::sortAtoms(System *system)
{
    //This function clears the neighbor lists of atoms, and the assigns the atoms again in the correct neighborhoods

    clear(system);      //Clears the nighborLists of atoms

    int i,j,k;

//    cout << "systemsize = " << system->systemSize().x << "; #Nabolag = " << m_noOfNeighbors << "; neighborLength = "
//         << system->systemSize().x/m_noOfNeighbors << endl;

    for(int n = 0 ; n < int(system->atoms().size()) ; n++ )
    {
        //Each atom needs to be placed in it's neighborhood

        i = int (system->atoms()[n]->position.x / system->systemSize().x * m_noOfNeighbors );
        j = int (system->atoms()[n]->position.y / system->systemSize().y * m_noOfNeighbors );
        k = int (system->atoms()[n]->position.z / system->systemSize().z * m_noOfNeighbors );

//        cout << "Atom " << n << endl;
//        cout << "Har position " <<  system->atoms()[n]->position << " og er i box " << i << j << k << endl;


        m_neighborList[i][j][k].push_back(system->atoms()[n]);


    }





}

void NeighborList::clear(System *system)
{
  //Clear all neighborhoods for atoms, to be used at the start of assigning them

    for(int i = 0 ; i < int( m_neighborList.size() ); i++)
    {
        for(int j = 0 ; j < int( m_neighborList[0].size() ); j++)
        {
            for(int k = 0 ; k < int( m_neighborList[0][0].size() ); k++)
                m_neighborList[i][j][k].clear();
        }
    }

}



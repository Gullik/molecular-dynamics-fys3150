#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <math/vec3.h>
#include <system.h>

using CompPhys::vec3;

class NeighborList
{
private:
    vector<vector<vector<vector<Atom*> > > >  m_neighborList;
    double m_noOfNeighbors;
public:
    NeighborList(double neighborSize, vec3 systemSize);
    vector<vector<vector<vector<Atom*> > > >  &neigborList() {return m_neighborList;}
    double &noOfNeighbors() {return m_noOfNeighbors;};
    void sortAtoms(System *system);
    void clear(System *system);
    void ghostCells(System *system); //Not implemented yet, not sure if I will.


};

#endif // NEIGHBORLIST_H

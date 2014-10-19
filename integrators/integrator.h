#pragma once

//enum Integrators {VelocityVerlet = 0};        //Not sure why it is here and what it does, it stops i
                                                // from compiling

class System;
class Integrator
{
public:
    Integrator();
    virtual ~Integrator() { }
    virtual void integrate(System* system, double timestep) = 0;
};

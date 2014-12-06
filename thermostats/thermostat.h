#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include <string>
#include <vector>
#include <system.h>

class System;

class Thermostat
{
private:
    bool m_switch;
protected:
    double m_wantedTemperature;
public:
    Thermostat();
    virtual ~Thermostat() {}
    double wantedTemperature() {return m_wantedTemperature; }
    virtual void applyThermostat(System *system, double dt) = 0;
    bool status() {return m_switch;}
    void turn_off() {m_switch = false;}
    void turn_on() {m_switch = true;}
};

#endif // THERMOSTAT_H

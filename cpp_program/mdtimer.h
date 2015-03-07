#ifndef MDTIMER_H
#define MDTIMER_H

#include <time.h>

class MDTimer
{
public:
    MDTimer();
    void PrintTimes();

    //Moving timer
    double move_start;
    double move_end;
    double current_move_time;
    void move_start_timer() {move_start = clock(); }
    void move_end_timer() ;

    //Force timer
    double force_start;
    double force_end;
    double current_force_time;
    void force_start_timer() {force_start = clock(); }
    void force_end_timer() ;

    //Thermostat timer
    double thermo_start;
    double thermo_end;
    double current_thermo_time;
    void thermo_start_timer() {thermo_start = clock(); }
    void thermo_end_timer() ;

    //Write timer
    double write_start;
    double write_end;
    double current_write_time;
    void write_start_timer() {write_start = clock(); }
    void write_end_timer() ;

    //Sampler timer
    double sampler_start;
    double sampler_end;
    double current_sampler_time;
    void sampler_start_timer() {sampler_start = clock(); }
    void sampler_end_timer() ;

    //Total timer
    double total_start;
    double total_end;
    void start_timer() {total_start = clock();}
    void end_timer() ;

};

#endif // MDTIMER_H

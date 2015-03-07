#include "mdtimer.h"
#include <time.h>
#include <iostream>


using namespace std;


MDTimer::MDTimer():
current_move_time(0),
current_force_time(0),
current_thermo_time(0),
current_write_time(0)
{

}

void MDTimer::PrintTimes()
{
    double total = 1.*(total_end - total_start)/CLOCKS_PER_SEC;
    double forceTime = 1.*( current_force_time )/CLOCKS_PER_SEC;
    double moveTime = 1.*( current_move_time )/CLOCKS_PER_SEC;
    double thermoTime = 1.*( current_thermo_time )/CLOCKS_PER_SEC;
    double writeTime = 1.*( current_write_time )/CLOCKS_PER_SEC;
    double sampleTime = 1.*( current_sampler_time )/CLOCKS_PER_SEC;




    cout << endl << fixed << "   Total time: "  << total << " s" << endl <<  endl;
    cout << "   Force time  : "   << forceTime << " s \t" << 100*forceTime/total << " %" << endl;
    cout << "   Move time   : "   << moveTime  << " s \t" << 100*moveTime/total << " %" << endl;
    cout << "   Thermo time : "   << thermoTime  << " s \t" << 100*thermoTime/total << " %" << endl;
    cout << "   Write time  : "   << writeTime << " s \t" << 100*writeTime/total << " %" << endl;
    cout << "   Sample time : "   << sampleTime << " s \t" << 100*sampleTime/total << " %" << endl;


    return;
}

void MDTimer::move_end_timer()
{
    move_end = clock();

    current_move_time += move_end - move_start;
}

void MDTimer::force_end_timer()
{
   force_end = clock();

    current_force_time += force_end - force_start;
}


void MDTimer::write_end_timer()
{
   write_end = clock();

    current_write_time += write_end - write_start;
}


void MDTimer::thermo_end_timer()
{
   thermo_end = clock();

    current_thermo_time += thermo_end - thermo_start;
}

void MDTimer::sampler_end_timer()
{
   sampler_end = clock();

    current_sampler_time += sampler_end - sampler_start;
}



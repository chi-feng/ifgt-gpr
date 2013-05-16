#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER

#include <time.h> /* clock_t, clock, CLOCKS_PER_SEC */

clock_t t;

void timer_start() {
    t = clock();
}

void timer_stop() {
    t = clock() - t;
}

float timer_value() { 
    return ((float)t) / CLOCKS_PER_SEC;
}

#endif
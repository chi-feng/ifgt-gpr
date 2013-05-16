#ifndef EXECUTION_TIMER
#define EXECUTION_TIMER

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

clock_t t;

double previous_eta = 0; 

void print_eta(int completed, int total, double seconds) {
    double eta = seconds * (total - completed);
    if (previous_eta == 0) previous_eta = eta;
    if (eta > previous_eta && eta < previous_eta * 1.2) eta = previous_eta;
    if (eta > 3600) {
        printf("ETA: %5.1f h",  eta / 3600);
    } else if (eta > 60) {
        printf("ETA: %5.1f m",  eta / 60);
    } else {
        printf("ETA: %5.1f s",  eta);
    }
    printf(" [%05d/%05d]\r", completed, total);
    fflush(stdout);
    previous_eta = eta;
}

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
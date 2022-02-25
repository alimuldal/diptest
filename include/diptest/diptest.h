#ifndef INCLUDE_DIPTEST_DIPTEST_H_
#define INCLUDE_DIPTEST_DIPTEST_H_
#define UNUSED(expr) do { (void)(expr); } while (0)

#include <stdio.h>

/* Subroutine */
void diptst(
    const double x[],
    const int *n_,
    double *dip,
    int *lo_hi,
    int *ifault,
    int *gcm,
    int *lcm,
    int *mn,
    int *mj,
    const int *min_is_0,
    const int *debug
);
#endif  // INCLUDE_DIPTEST_DIPTEST_H_

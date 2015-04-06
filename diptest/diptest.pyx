

import numpy as np
cimport numpy as np


cdef extern from "_dip.h" nogil:

    void diptst(const double* x,
                const int *n_,
                double *dip,
                int *lo_hi,
                int *ifault,
                int *gcm,
                int *lcm,
                int *mn,
                int *mj,
                const int *min_is_0,
                const int *debug)


def diptest(x, min_is_0=True, full_output=False, debug=0):
    """
    Hartigan & Hartigan's dip test statistic for unimodality

    Arguments:
    -----------
    x:              [n,] array  containing the input data

    min_is_0:       boolean, if True the minimum value of the test statistic is
                    allowed to be zero in cases where n < 2 or all values in x
                    are identical

    full_output:    boolean, see below

    debug:          int, 0 <= debug <= 2, print debugging messages

    Returns:
    -----------
    dip:    double, the dip statistic

    [res]:  dict, returned if full_output == True. contains the following
            fields:

            xs:     sorted input data as doubles
            n:      len(x)
            dip:    dip statistic
            lo:     indices of lower end of modal interval
            hi:     indices of upper end of modal interval
            xl:     lower end of modal interval
            xu:     upper end of modal interval
            gcm:    (last-used) indices of the greatest concave majorant
            lcm:    (last-used) indices of the least concave majorant

    """

    cdef:
        double[:] x_
        int n = x.shape[0]

        double dip = np.nan
        int[:] lo_hi = np.empty(4, dtype=np.int32)
        int ifault = 0
        int[:] gcm = np.empty(n, dtype=np.int32)
        int[:] lcm = np.empty(n, dtype=np.int32)
        int[:] mn = np.empty(n, dtype=np.int32)
        int[:] mj = np.empty(n, dtype=np.int32)

        int min_is_0_ = min_is_0
        int debug_ = debug

    # cast to double, force a copy to avoid inplace sort
    x = np.array(x, dtype=np.double, copy=True)

    # input needs to be sorted in ascending order
    x.sort()

    x_ = x

    diptst(&x_[0], &n, &dip, &lo_hi[0], &ifault, &gcm[0], &lcm[0], &mn[0],
           &mj[0], &min_is_0_, &debug_)

    if full_output:
        res_dict = {
            'xs':np.array(x_),
            'n':n,
            'dip':dip,
            'lo':lo_hi[0],
            'hi':lo_hi[1],
            'xl':x_[lo_hi[0]],
            'xu':x_[lo_hi[1]],
            'gcm':np.array(gcm[:lo_hi[2]]),
            'lcm':np.array(lcm[:lo_hi[3]]),
            'mn':np.array(mn),
            'mj':np.array(mj),
        }

        return dip, res_dict

    else:

        return dip
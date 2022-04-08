/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */

#include <diptest/wrapper.hpp>

namespace py = pybind11;

namespace diptest {

namespace details {

std::unique_ptr<double[]> std_uniform(const int64_t n_boot, const int64_t n, const int64_t seed) {
    std::mt19937_64 rng;
    rng.seed(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int64_t size = n_boot * n;
    std::unique_ptr<double[]> sample(new double[size]);
    for (int64_t i = 0; i < size; i++) {
        sample[i] = dist(rng);
    }
    return sample;
}

double diptest(const double* x_ptr, int N, int allow_zero, int debug) {
    int ifault = 0;
    int lo_hi[4] = {0, 0, 0, 0};
    std::unique_ptr<int[]> gcm(new int[N]);
    std::unique_ptr<int[]> lcm(new int[N]);
    std::unique_ptr<int[]> mn(new int[N]);
    std::unique_ptr<int[]> mj(new int[N]);

    double dip = diptst(x_ptr, N, &lo_hi[0], &ifault, gcm.get(), lcm.get(), mn.get(), mj.get(), allow_zero, debug);

    if (ifault == 1) {
        throw std::runtime_error("N must be >= 1.");
    } else if (ifault == 2) {
        throw std::runtime_error("x must be sorted in ascending error.");
    }
    return dip;
}  // diptest

}  // namespace details

double diptest(const py::array_t<double>& x, int allow_zero, int debug) {
    return details::diptest(x.data(), x.size(), allow_zero, debug);
}  // diptest

py::dict diptest_full(const py::array_t<double>& x, int allow_zero, int debug) {
    const double* x_ptr = x.data();
    int N = x.size();
    int ifault = 0;
    int lo_hi[4] = {0, 0, 0, 0};

    auto gcm = py::array_t<int>(N);
    auto lcm = py::array_t<int>(N);
    int* gcm_ptr = gcm.mutable_data();
    int* lcm_ptr = lcm.mutable_data();

    std::unique_ptr<int[]> mn(new int[N]);
    std::unique_ptr<int[]> mj(new int[N]);
    int* mn_ptr = mn.get();
    int* mj_ptr = mj.get();

    double dip = diptst(x_ptr, N, &lo_hi[0], &ifault, gcm_ptr, lcm_ptr, mn_ptr, mj_ptr, allow_zero, debug);

    if (ifault == 1) {
        throw std::runtime_error("N must be >= 1.");
    } else if (ifault == 2) {
        throw std::runtime_error("x must be sorted in ascending error.");
    }

    using namespace pybind11::literals;  // to bring in the `_a` literal NOLINT
    return py::dict(
        "dip"_a = dip,
        "lo"_a = lo_hi[0],
        "hi"_a = lo_hi[1],
        "xl"_a = x.at(lo_hi[0]),
        "xu"_a = x.at(lo_hi[1]),
        "_gcm"_a = gcm,
        "_lcm"_a = lcm,
        "_lh_2"_a = lo_hi[2],
        "_lh_3"_a = lo_hi[3]);
}  // diptest_full

double
diptest_pval(const double dipstat, const int64_t n, const int64_t n_boot, int allow_zero, int debug, int64_t seed) {
    std::random_device rd;
    std::mt19937_64 rng;
    if (seed == 0) {
        rng.seed(rd());
    } else {
        rng.seed(seed);
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double dip;
    int ifault = 0;
    int lo_hi[4] = {0, 0, 0, 0};
    std::unique_ptr<int[]> gcm(new int[n]);
    std::unique_ptr<int[]> lcm(new int[n]);
    std::unique_ptr<int[]> mn(new int[n]);
    std::unique_ptr<int[]> mj(new int[n]);
    std::unique_ptr<int[]> dips(new int[n_boot]);
    std::unique_ptr<double[]> sample(new double[n]);

    double* r_sample = sample.get();
    double* sample_end = r_sample + n;

    for (int64_t i = 0; i < n_boot; i++) {
        for (int64_t i = 0; i < n; i++) {
            r_sample[i] = dist(rng);
        }
        std::sort(r_sample, sample_end);
        dip = diptst(r_sample, n, &lo_hi[0], &ifault, gcm.get(), lcm.get(), mn.get(), mj.get(), allow_zero, debug);
        dips[i] = dipstat <= dip;
    }
    double p_val = std::accumulate(dips.get(), dips.get() + n_boot, 0.0) / n_boot;
    return p_val;
}  // diptest_pval

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
double diptest_pval_mt(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    int64_t seed,
    size_t n_threads) {
    std::random_device rd;
    if (seed == 0) {
        seed = rd();
    }

    // allocate whole sample block in one go.
    std::unique_ptr<double[]> sample = details::std_uniform(n_boot, n, seed);
    std::unique_ptr<bool[]> dips(new bool[n_boot]);

#pragma omp parallel num_threads(n_threads) shared(dips)
    {
        int ifault = 0;
        double* p_sample_end;
        double* p_sample;
        std::unique_ptr<int[]> lo_hi(new int[n]);
        std::memset(lo_hi.get(), 0, 4);
        std::unique_ptr<int[]> gcm(new int[n]);
        std::unique_ptr<int[]> lcm(new int[n]);
        std::unique_ptr<int[]> mn(new int[n]);
        std::unique_ptr<int[]> mj(new int[n]);

#pragma omp for
        for (int64_t i = 0; i < n_boot; i++) {
            // each thread get a memory block of size `n`
            // the below is bookkeeping to assign the correct block to each thread
            p_sample = sample.get() + (i * n);
            p_sample_end = p_sample + n;
            // sort the allocated block for this bootstrap sample
            std::sort(p_sample, p_sample_end);
            dips[i] = dipstat <= diptst(
                p_sample, n, lo_hi.get(), &ifault, gcm.get(), lcm.get(), mn.get(), mj.get(), allow_zero, debug
            );
        }
    }  // pragma parallel
    double p_val = std::accumulate(dips.get(), dips.get() + n_boot, 0.0) / n_boot;
    return p_val;
}  // diptest_pval_mt
#endif

namespace bindings {

void bind_diptest(py::module& m) {
    m.def("diptest", &diptest::diptest, py::arg("x"), py::arg("allow_zero") = 1, py::arg("debug") = 0);
}

void bind_diptest_full(py::module& m) {
    m.def("diptest_full", &diptest::diptest_full, py::arg("x"), py::arg("allow_zero") = 1, py::arg("debug") = 0);
}

void bind_diptest_pval(py::module& m) {
    m.def(
        "diptest_pval",
        &diptest::diptest_pval,
        py::arg("dipstat"),
        py::arg("n"),
        py::arg("n_boot") = 10000,
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0,
        py::arg("seed") = 0);
}

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
void bind_diptest_pval_mt(py::module& m) {
    m.def(
        "diptest_pval_mt",
        &diptest::diptest_pval_mt,
        py::arg("dipstat"),
        py::arg("n"),
        py::arg("n_boot") = 10000,
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0,
        py::arg("seed") = 0,
        py::arg("n_threads") = 4);
}
#endif

}  // namespace bindings
}  // namespace diptest

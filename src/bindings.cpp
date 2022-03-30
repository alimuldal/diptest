/* bindings.cpp -- Python bindings for diptest
 * Copyright 2022 Ralph Urlus
 */
#include <pybind11/pybind11.h>
#include <diptest/wrapper.hpp>

namespace diptest {
namespace bindings {

PYBIND11_MODULE(EXTENSION_MODULE_NAME, m) {
    bind_diptest(m);
    bind_diptest_full(m);
    bind_diptest_pval(m);
#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
    bind_diptest_pval_mt(m);
#endif

#ifndef OS_WIN
#ifdef DIPTEST_VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(DIPTEST_VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
#endif
}

}  // namespace bindings
}  // namespace diptest

#pragma once

#include <nanobind/nanobind.h>

#include "eckit/config/Configuration.h"
#include "atlas/util/Config.h"

namespace nb = ::nanobind;

namespace atlas4py {

    // Bindings for atlas::util::Config
    void bind_Config( nb::module_& m );

    // Convert eckit::Configuration to a corresponding Python object
    nb::object make_object( eckit::Configuration const& v );

} // namespace atlas4py

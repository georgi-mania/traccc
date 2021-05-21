#pragma once

#include "edm/cell.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"

#include <map>

namespace traccc {

    struct result {
        traccc::host_measurement_container measurements;
        traccc::host_spacepoint_container spacepoints;

        result() = delete;
        result(traccc::host_measurement_container mc, traccc::host_spacepoint_container sc) : measurements(mc), spacepoints(sc) {};
    };

    using geometry = std::map<traccc::geometry_id, traccc::transform3>;
    using demonstrator_input = std::map<size_t , traccc::host_cell_container>;
    using demonstrator_result = std::map<size_t, result>;
}

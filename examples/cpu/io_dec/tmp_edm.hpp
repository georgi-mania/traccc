#pragma once

#include "edm/cell.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"

#include <map>
#include <vecmem/memory/host_memory_resource.hpp>

namespace traccc {

    vecmem::host_memory_resource resource;

    struct result {
        traccc::host_measurement_container measurements;
        traccc::host_spacepoint_container spacepoints;

        void setMeasurements(const host_measurement_container &measurements) {
            result::measurements = measurements;
        }

        void setSpacepoints(const host_spacepoint_container &spacepoints) {
            result::spacepoints = spacepoints;
        }
    };

    using geometry = std::map<traccc::geometry_id, traccc::transform3>;
    using demonstrator_input = vecmem::vector<traccc::host_cell_container>;
    using demonstrator_result = vecmem::vector<traccc::result>;
}

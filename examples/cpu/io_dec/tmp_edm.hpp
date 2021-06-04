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

        void addMeasurement(const cell_module& module, const std::optional<measurement>& m) {
            if (m.has_value()) {
                traccc::host_measurement_collection holder;
                holder.push_back(std::move(m.value()));

                #pragma omp critical (results_measurement)
                {
                    measurements.headers.push_back(module);
                    measurements.items.push_back(std::move(holder));
                }
            }
        }

        void addSpacepoint(const cell_module& module, const std::optional<spacepoint>& s) {
            if (s.has_value()) {
                traccc::host_spacepoint_collection holder;
                holder.push_back(std::move(s.value()));

                #pragma omp critical (results_spacepoints)
                {
                    spacepoints.headers.push_back(module.module);
                    spacepoints.items.push_back(std::move(holder));
                };
            }
        }
    };

    using geometry = std::map<traccc::geometry_id, traccc::transform3>;
    using demonstrator_input = vecmem::vector<traccc::host_cell_container>;
    using demonstrator_result = vecmem::vector<traccc::result>;
}

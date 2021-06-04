#pragma once

#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"
#include <optional>

namespace traccc {

    struct spacepoint_formation {

        /// Callable operator for the space point formation, based on one single measurement
        ///
        /// @param measurements is the input measurement
        /// @return a measurement
        std::optional<spacepoint> operator()(const std::optional<measurement>& measurement, const cell_module& module) const {
            std::optional<spacepoint> maybe_spacepoint;

            if (measurement.has_value()) {
                this->operator()(measurement.value(), module, maybe_spacepoint);
            } else  {
                maybe_spacepoint = std::nullopt;
            }

            return maybe_spacepoint;
        }

        /// Callable operator for the space point formation, based on a measurement
        ///
        /// @param measurement is the input measurement
        /// @return a measurement
        void operator()(const measurement& m, const cell_module& module, std::optional<spacepoint>& maybe_spacepoint) const {
            // Run the algorithm
            spacepoint s;
            point3 local_3d = {m.local[0], m.local[1], 0.};
            s.global = module.placement.point_to_global(local_3d);
            maybe_spacepoint = s;
        }

    };

}

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"
#include "edm/cell.hpp"
#include "edm/cluster.hpp"
#include "edm/measurement.hpp"
#include "definitions/algebra.hpp"

#include <optional>
#include <iostream>

namespace traccc {

    /// Connected component labeling.
    struct measurement_creation {

        /// Callable operator for the connected component, based on a single cluster from
        /// one  module
        ///
        /// @param cluster are the input cells into the connected component, per module, unordered
        ///
        /// @return a measurement (if it satisfies the threshold condition)
        std::optional<measurement> operator()(const cluster &cluster, const cell_module &module) const {
            std::optional<measurement> maybe_measurement;
            if (!cluster.cells.empty()) {
                this->operator()(cluster, module, maybe_measurement);
            } else {
                maybe_measurement = std::nullopt;
            }

            return maybe_measurement;
        }

        void
        operator()(const cluster &cluster, const cell_module &module, std::optional<measurement>& maybe_measurement) const {
            // Run the algorithm
            point2 p = {0., 0.};
            variance2 v = {0., 0.};
            scalar totalWeight = 0.;

            for (const auto &cell : cluster.cells) {
                scalar weight = traccc::signal(cell.activation);
                if (weight > cluster.threshold) {
                    totalWeight += cell.activation;

                    auto cell_position = module.pixel(cell.channel0, cell.channel1);
                    p = p + weight * cell_position;
                    point2 square_pos = {cell_position[0] * cell_position[0],
                                         cell_position[1] * cell_position[1]};
                    v = v + weight * square_pos;
                }
            }

            if (totalWeight > 0.) {
                measurement m;
                // normalize the cell position
                m.local = 1. / totalWeight * p;
                // normalize the variance
                m.variance = 1. / totalWeight * v;
                auto pitch = module.pixel.get_pitch();
                // plus pitch^2 / 12
                m.variance = m.variance + point2{pitch[0] * pitch[0] / 12,
                                                           pitch[1] * pitch[1] / 12};
                // minus <x>^2
                m.variance = m.variance - point2{m.local[0] * m.local[0],
                                                   m.local[1] * m.local[1]};
                // @todo add variance estimation
                //measurements.push_back(std::move(m));
                maybe_measurement = m;
            } else {
                maybe_measurement = std::nullopt;
            }

        }
    };

} // namespace traccc

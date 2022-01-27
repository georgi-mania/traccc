/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"
#include <vector>
#include <limits>

namespace traccc {

    using channel_id = unsigned int;

    /// A cell definition: 
    ///
    /// maximum two channel identifiers
    /// and one activiation value, such as a time stamp
    struct cell {
        channel_id channel0 = 0;
        channel_id channel1 = 0;
        scalar activation = 0.;
        scalar time = 0.;
    };

    /// A cell collection: 
    ///
    /// it remembers the moduleentifier and also 
    /// keeps track of the cell ranges for chosing optimal
    /// algorithm.
    struct cell_collection { 

        event_id event = 0;
        geometry_id module = 0;
        transform3 placement = transform3{};

        std::vector<cell> items{};
        std::vector<channel_id> range0 = {std::numeric_limits<channel_id>::max(), 0};
        std::vector<channel_id> range1 = {std::numeric_limits<channel_id>::max(), 0};
    };

    using cell_container = std::vector<cell_collection>;

}


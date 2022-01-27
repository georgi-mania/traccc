/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <stdint.h>

namespace traccc {

    using scalar = float;
    using geometry_id = uint64_t;
    using event_id = uint64_t;

    using vector2 = std::vector<scalar>;
    using point2 = std::vector<scalar>;
    using variance2 = std::vector<scalar>;
    using point3 = std::vector<scalar>;
    using vector3 = std::vector<scalar>;
    using variance3 = std::vector<scalar>;

}

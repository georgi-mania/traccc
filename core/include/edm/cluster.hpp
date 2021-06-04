/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/cell.hpp"

#include <vector>
#include <functional>

namespace traccc {
    
    /// A cluster definition: 
    ///
    /// a list of cells that make up the cluster
    struct cluster {
        std::vector<cell> cells;
        float threshold = 0.;
    };

    using position_estimation = std::function<vector2(channel_id,channel_id)>;
    using signal_modeling = std::function<float(float)>;

    /// A cluster collection which carries the geometry_id, the clusters
    /// and the additional information to create the cluster position 
    /// from the channel id;
    struct cluster_collection {

        event_id event = 0;
        geometry_id module = 0;
        transform3 placement = transform3{};

        std::vector<cluster> items;
        
        position_estimation position_from_cell 
            = [](channel_id ch0,channel_id ch1) -> vector2 
                {return {static_cast<float>(ch0), static_cast<float>(ch1)}; };

        float threshold = 0.;
        signal_modeling signal 
            = [](float signal_in) -> float {return signal_in; };
        
    };

    using cluster_container = std::vector<cluster_collection>;

    vector2 position_from_cell (channel_id ch0,channel_id ch1){
        return {static_cast<float>(ch0), static_cast<float>(ch1)};
    }

    float signal (float signal_in) {
        return signal_in;
    }

}

//TODO Georgiana:
// 1. remove cluster_collection; functions exposed in namespace directly;
//  - in measurement_formation use info from the module and generic fns
//  - in spacepoint_formation use info from module
// 2. measurement- & spacepoint formation operator() -> works on 1 element instead of collection
//  -  aggregation in calling function
// 3. pass module info as fn closure to measurement- and spacepoint formation fns
// 4. validate results against sequential code

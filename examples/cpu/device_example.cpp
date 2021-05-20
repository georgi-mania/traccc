/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "edm/cell.hpp"
#include "edm/cluster.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"
#include "geometry/pixel_segmentation.hpp"
#include "algorithms/component_connection.hpp"
#include "algorithms/measurement_creation.hpp"
#include "algorithms/spacepoint_formation.hpp"
#include "csv/csv_io.hpp"

#include <vecmem/memory/host_memory_resource.hpp>

#include <iostream>
#include <chrono>
#include "omp.h"

int seq_run(const std::string& detector_file, const std::string& cells_dir, unsigned int events)
{
    auto env_d_d = std::getenv("TRACCC_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set TRACCC_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d) + std::string("/");

    // Read the surface transforms
    std::string io_detector_file = data_directory + detector_file;
    traccc::surface_reader sreader(io_detector_file, {"geometry_id", "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw", "rot_zu", "rot_zv", "rot_zw"});
    auto surface_transforms = traccc::read_surfaces(sreader);

    // Algorithms
    traccc::component_connection cc;
    traccc::measurement_creation mt;
    traccc::spacepoint_formation sp;

    // Output stats
    uint64_t n_cells = 0;
    uint64_t m_modules = 0;
    uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_space_points = 0;

    // Memory resource used by the EDM.
    vecmem::host_memory_resource resource;

    std::string event_string = "000000000";
    std::string event_number = std::to_string(1);
    event_string.replace(event_string.size()-event_number.size(), event_number.size(), event_number);

    std::string io_cells_file = data_directory+cells_dir+std::string("/event")+event_string+std::string("-cells.csv");
    traccc::cell_reader creader(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});
    traccc::host_cell_container cells_per_event = traccc::read_cells(creader, resource, surface_transforms);

    traccc::cell_module module = cells_per_event.headers[0];
    module.pixel.min_center_y=4.5;
    module.pixel.pitch_x=1.0;

    int size = 3;

    std::array<traccc::scalar , 3> array = {0.1, 0.2, 0.3};
    traccc::scalar* arr_pointer = array.data();

    vecmem::vector<traccc::scalar> v({1.0, 2.0, 3.0}, &resource);
    traccc::scalar* v_pointer = v.data();

    vecmem::jagged_vector<traccc::scalar> jv({vecmem::vector<traccc::scalar>({1,2,3}),
                                                vecmem::vector<traccc::scalar>({10, 20, 30})},
                                             &resource);
    traccc::scalar* jv_pointer = &jv.data()->at(0);
    for (int i=0; i<6; i++)
        printf("%f ", jv_pointer[i]);


#pragma omp target map(to:module, arr_pointer[0:size], v_pointer[0:size], jv_pointer[0:6])
    {
           int device = omp_is_initial_device();
           if (!device) {
               printf("On GPU: \n");
               printf("Module = %d\n", module.module);
               printf("Event  = %d\n", module.event);
               printf("Range0 = %lu\n", module.range0);
               printf("Range1 = %lu\n", module.range1);
               printf("Pixelx = %f\n", module.pixel.pitch_x);
               printf("Pixely = %f\n", module.pixel.pitch_y);
               printf("Minx = %f\n", module.pixel.min_center_x);
               printf("Miny = %f\n", module.pixel.min_center_y);

               printf("\nstd::array:\n");
               for (int i=0; i<size; i++)
                   printf("%f ", arr_pointer[i]);

               printf("\nvecmem::vector:\n");
               for (int i=0; i<size; i++)
                   printf("%f ", v_pointer[i]);

               printf("\nvecmem::jagged_vector:\n");
               for (int i=0; i<6; i++)
                   printf("%f ", jv_pointer[i]);

               printf("\nSuccess!\n");
           //    printf("Placement size= %ld\n ", module.placement._data.size());

           } else {
               printf("Still on CPU \n");
           }
    }
    /*
#pragma omp parallel for
    // Loop over events
    for (unsigned int event = 0; event < events; ++event){

        // Read the cells from the relevant event file
        std::string event_string = "000000000";
        std::string event_number = std::to_string(event);
        event_string.replace(event_string.size()-event_number.size(), event_number.size(), event_number);

        std::string io_cells_file = data_directory+cells_dir+std::string("/event")+event_string+std::string("-cells.csv");
        traccc::cell_reader creader(io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});
        traccc::host_cell_container cells_per_event = traccc::read_cells(creader, resource, surface_transforms);
        m_modules += cells_per_event.headers.size();

        // Output containers
        traccc::host_measurement_container measurements_per_event;
        traccc::host_spacepoint_container spacepoints_per_event;
        measurements_per_event.headers.reserve(cells_per_event.headers.size());
        measurements_per_event.items.reserve(cells_per_event.headers.size());
        spacepoints_per_event.headers.reserve(cells_per_event.headers.size());
        spacepoints_per_event.items.reserve(cells_per_event.headers.size());

#pragma omp parallel for
        for (std::size_t i = 0; i < cells_per_event.items.size(); ++i )
        {
            auto& module = cells_per_event.headers[i];
            module.pixel = traccc::pixel_segmentation{-8.425, -36.025, 0.05, 0.05};

            // The algorithmic code part: start
            traccc::cluster_collection clusters_per_module = cc(cells_per_event.items[i], cells_per_event.headers[i]);
            clusters_per_module.position_from_cell = module.pixel;

            traccc::host_measurement_collection measurements_per_module = mt(clusters_per_module, module);
            traccc::host_spacepoint_collection spacepoints_per_module = sp(module, measurements_per_module);
            // The algorithmnic code part: end

            n_cells += cells_per_event.items[i].size();
            n_clusters += clusters_per_module.items.size();
            n_measurements += measurements_per_module.size();
            n_space_points += spacepoints_per_module.size();

            measurements_per_event.items.push_back(std::move(measurements_per_module));
            measurements_per_event.headers.push_back(module);

            spacepoints_per_event.items.push_back(std::move(spacepoints_per_module));
            spacepoints_per_event.headers.push_back(module.module);
        }

        traccc::measurement_writer mwriter{std::string("event")+event_number+"-measurements.csv"};
        for (size_t i=0; i<measurements_per_event.items.size(); ++i){
            auto measurements_per_module = measurements_per_event.items[i];
            auto module = measurements_per_event.headers[i];
            for (const auto& measurement : measurements_per_module){
                const auto& local = measurement.local;
                mwriter.append({ module.module, local[0], local[1], 0., 0.});
            }
        }

        traccc::spacepoint_writer spwriter{std::string("event")+event_number+"-spacepoints.csv"};
        for (size_t i=0; i<spacepoints_per_event.items.size(); ++i){
            auto spacepoints_per_module = spacepoints_per_event.items[i];
            auto module = spacepoints_per_event.headers[i];

            for (const auto& spacepoint : spacepoints_per_module){
                const auto& pos = spacepoint.global;
                spwriter.append({ module, pos[0], pos[1], pos[2], 0., 0., 0.});
            }
        }

    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_cells << " cells from " << m_modules << " modules" << std::endl;
    std::cout << "- created " << n_clusters << " clusters. " << std::endl;
    std::cout << "- created " << n_measurements << " measurements. " << std::endl;
    std::cout << "- created " << n_space_points << " space points. " << std::endl;
*/
    return 0;
}

// The main routine
//
int main(int argc, char *argv[])
{
    if (argc < 4){
        std::cout << "Not enough arguments, minimum requirement: " << std::endl;
        std::cout << "./device_example <detector_file> <cell_directory> <events>" << std::endl;
        return -1;
    }

    auto detector_file = std::string(argv[1]);
    auto cell_directory = std::string(argv[2]);
    auto events = std::atoi(argv[3]);

    std::cout << "Running ./device_exammple " << detector_file << " " << cell_directory << " " << events << std::endl;
    auto start = std::chrono::system_clock::now();
    auto result = seq_run(detector_file, cell_directory, events);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> diff = end-start;
    std::cout << "Execution time: " << diff.count() << " sec." << std::endl;
    return 0;
}


/*
 * # OMP device offload example
set(GCC_GPU_OFFLOAD "-foffload=nvptx-none")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_GPU_OFFLOAD}")
add_executable (device_example device_example.cpp)
target_link_libraries (device_example LINK_PUBLIC traccc::core traccc::io vecmem::core OpenMP::OpenMP_CXX)
 */
/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <chrono>
#include <iomanip>
#include <iostream>

// vecmem
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// io
#include "io/csv.hpp"
#include "io/utils.hpp"
#include "io/reader.hpp"
#include "io/writer.hpp"

// algorithms
#include "clusterization/clusterization_algorithm.hpp"
#include "cuda/track_finding/seeding_algorithm.hpp"
#include "track_finding/seeding_algorithm.hpp"

int seq_run(const std::string& detector_file, const std::string& cells_dir,
            unsigned int events, bool skip_cpu) {

    // Read the surface transforms
    auto surface_transforms = traccc::read_geometry(detector_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_internal_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_cuda = 0;

    // Elapsed time
    float wall_time(0);
    float file_reading_cpu(0);
    float clusterization_cpu(0);
    float seeding_cpu(0);
    float clusterization_cuda(0);
    float seeding_cuda(0);

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Memory resource used by the EDM.
    vecmem::cuda::managed_memory_resource mng_mr;

    traccc::clusterization_algorithm ca;
    traccc::cuda::seeding_algorithm sa_cuda(&mng_mr);
    traccc::seeding_algorithm sa(&host_mr);

    /*time*/ auto start_wall_time = std::chrono::system_clock::now();

    // Loop over events
    for (unsigned int event = 0; event < events; ++event) {

        // Read the cells from the relevant event file

        /*time*/ auto start_file_reading_cpu = std::chrono::system_clock::now();

        // Read the cells from the relevant event file
        std::string io_cells_file =
            traccc::data_directory() + cells_dir + "/" +
            traccc::get_event_filename(event, "-cells.csv");
        traccc::cell_reader creader(
            io_cells_file, {"geometry_id", "hit_id", "cannel0", "channel1",
                            "activation", "time"});
        traccc::host_cell_container cells_per_event =
            traccc::read_cells(creader, host_mr, &surface_transforms);
	
        /*time*/ auto end_file_reading_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_file_reading_cpu =
            end_file_reading_cpu - start_file_reading_cpu;
        /*time*/ file_reading_cpu += time_file_reading_cpu.count();

        /*-----------------------------
              Clusterization (cpu)
          -----------------------------*/

        /*time*/ auto start_clusterization_cpu =
            std::chrono::system_clock::now();
        auto ca_result = ca(cells_per_event);
        auto& measurements_per_event = ca_result.first;
        auto& spacepoints_per_event = ca_result.second;

        /*time*/ auto end_clusterization_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_clusterization_cpu =
            end_clusterization_cpu - start_clusterization_cpu;
        /*time*/ clusterization_cpu += time_clusterization_cpu.count();

        n_modules += cells_per_event.headers.size();
        n_cells += cells_per_event.total_size();
        n_measurements += measurements_per_event.total_size();
        n_spacepoints += spacepoints_per_event.total_size();

        /*----------------------------
          Seeding algorithm
          ----------------------------*/

        // cuda

        /*time*/ auto start_seeding_cuda = std::chrono::system_clock::now();

        auto sa_cuda_result = sa_cuda(spacepoints_per_event);
        auto& seeds_cuda = sa_cuda_result.second;
        n_seeds_cuda += seeds_cuda.headers[0];

        /*time*/ auto end_seeding_cuda = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_seeding_cuda =
            end_seeding_cuda - start_seeding_cuda;
        /*time*/ seeding_cuda += time_seeding_cuda.count();

        // cpu

        /*time*/ auto start_seeding_cpu = std::chrono::system_clock::now();

        traccc::host_seed_container seeds;
        traccc::host_internal_spacepoint_container internal_sp_per_event;
        if (!skip_cpu) {
            auto sa_result = sa(spacepoints_per_event);
            internal_sp_per_event = sa_result.first;
            seeds = sa_result.second;
            n_internal_spacepoints += internal_sp_per_event.total_size();
        }
        n_seeds += seeds.total_size();

        /*time*/ auto end_seeding_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_seeding_cpu =
            end_seeding_cpu - start_seeding_cpu;
        /*time*/ seeding_cpu += time_seeding_cpu.count();

        /*----------------------------------
          compare seeds from cpu and cuda
          ----------------------------------*/

        if (!skip_cpu) {
            int n_match = 0;
            for (auto seed : seeds.items[0]) {
                if (std::find(
                        seeds_cuda.items[0].begin(),
                        seeds_cuda.items[0].begin() + seeds_cuda.headers[0],
                        seed) !=
                    seeds_cuda.items[0].begin() + seeds_cuda.headers[0]) {
                    n_match++;
                }
            }
            float matching_rate = float(n_match) / seeds.headers[0];
            std::cout << "event " << std::to_string(event)
                      << " seed matching rate: " << matching_rate << std::endl;
        }

        /*------------
             Writer
          ------------*/

        if (!skip_cpu) {
	    traccc::write_measurements(event, measurements_per_event);
	    traccc::write_spacepoints(event, spacepoints_per_event);
	    traccc::write_internal_spacepoints(event, internal_sp_per_event);
	    traccc::write_seeds(event, seeds);
        }
    }

    /*time*/ auto end_wall_time = std::chrono::system_clock::now();
    /*time*/ std::chrono::duration<double> time_wall_time =
        end_wall_time - start_wall_time;

    /*time*/ wall_time += time_wall_time.count();

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_spacepoints << " spacepoints from "
              << n_modules << " modules" << std::endl;
    std::cout << "- created        " << n_cells << " cells           "
              << std::endl;
    std::cout << "- created        " << n_measurements << " meaurements     "
              << std::endl;
    std::cout << "- created        " << n_spacepoints << " spacepoints     "
              << std::endl;
    std::cout << "- created        " << n_internal_spacepoints
              << " internal spacepoints" << std::endl;

    std::cout << "- created (cpu)  " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (cuda) " << n_seeds_cuda << " seeds" << std::endl;
    std::cout << "==> Elpased time ... " << std::endl;
    std::cout << "wall time           " << std::setw(10) << std::left
              << wall_time << std::endl;
    std::cout << "file reading (cpu)       " << std::setw(10) << std::left
              << file_reading_cpu << std::endl;
    std::cout << "clusterization_time (cpu)" << std::setw(10) << std::left
              << clusterization_cpu << std::endl;
    std::cout << "seeding_time (cpu)       " << std::setw(10) << std::left
              << seeding_cpu << std::endl;
    std::cout << "seeding_time (cuda)      " << std::setw(10) << std::left
              << seeding_cuda << std::endl;

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Not enough arguments, minimum requirement: " << std::endl;
        std::cout << "./seq_example <detector_file> <hit_directory> "
                     "<events> <skip_cpu>"
                  << std::endl;
        return -1;
    }

    auto detector_file = std::string(argv[1]);
    auto hit_directory = std::string(argv[2]);
    auto events = std::atoi(argv[3]);
    bool skip_cpu = std::atoi(argv[4]);

    std::cout << "Running ./seq_example " << detector_file << " "
              << hit_directory << " " << " " << events << " "
              << skip_cpu << " " << std::endl;
    return seq_run(detector_file, hit_directory, events, skip_cpu);
}

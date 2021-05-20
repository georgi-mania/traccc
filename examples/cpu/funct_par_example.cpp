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
#include <functional>
#include <chrono>

using namespace std::placeholders;

// Memory resource used by the EDM.
vecmem::host_memory_resource resource;

struct Results {
    vecmem::vector<traccc::host_measurement_container> measurements;
    vecmem::vector<traccc::host_spacepoint_container> spacepoints;

    Results(int size) {
        measurements.reserve(size);
        spacepoints.reserve(size);
    }
};

const std::string data_directory() {
    auto env_d_d = std::getenv("TRACCC_TEST_DATA_DIR");
    if (env_d_d == nullptr) {
        throw std::ios_base::failure("Test data directory not found. Please set TRACCC_TEST_DATA_DIR.");
    }
    std::string data_dir = std::string(env_d_d);
    return data_dir.append("/");
}

std::string get_event_filename(unsigned int event) {
    std::string event_string{"000000000"};
    std::string event_number = std::to_string(event);
    event_string.replace(event_string.size() - event_number.size(), event_number.size(), event_number);
    return event_string;
}

std::map<traccc::geometry_id, traccc::transform3> read_geometry(const std::string &detector_file) {
    // Read the surface transforms
    std::string io_detector_file = data_directory() + detector_file;
    traccc::surface_reader sreader(io_detector_file,
                                   {"geometry_id", "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw", "rot_zu", "rot_zv",
                                    "rot_zw"});
    return traccc::read_surfaces(sreader);
}

traccc::host_cell_container read_cells_from_event(unsigned int event, const std::string &cells_dir,
                                                  std::map<traccc::geometry_id, traccc::transform3> surface_transforms) {
    // Read the cells from the relevant event file
    std::string io_cells_file =
            data_directory() + cells_dir + std::string("/event") + get_event_filename(event) +
            std::string("-cells.csv");
    traccc::cell_reader creader(io_cells_file,
                                {"geometry_id", "hit_id", "cannel0", "channel1", "activation", "time"});
    return traccc::read_cells(creader, resource, surface_transforms);
}

vecmem::vector<traccc::host_cell_container>
read(int events, const std::string &detector_file, const std::string &cell_directory) {
    auto geometry = read_geometry(detector_file);
    auto readFn = std::bind(read_cells_from_event, _1, cell_directory, geometry);

    vecmem::vector<traccc::host_cell_container> input_data;
    input_data.reserve(events);

#pragma omp parallel for
    for (int event = 0; event < events; ++event) {
        traccc::host_cell_container cells_per_event = readFn(event);
        input_data.push_back(std::move(cells_per_event));
    }

    return input_data;
}

void write_measurements(unsigned int event, traccc::host_measurement_container &measurements_per_event) {
    traccc::measurement_writer mwriter{std::string("event") + std::to_string(event) + "-measurements.csv"};
    for (size_t i = 0; i < measurements_per_event.items.size(); ++i) {
        auto measurements_per_module = measurements_per_event.items[i];
        auto module = measurements_per_event.headers[i];
        for (const auto &measurement : measurements_per_module) {
            const auto &local = measurement.local;
            mwriter.append({module.module, local[0], local[1], 0., 0.});
        }
    }
}

void write_spacepoints(unsigned int event, traccc::host_spacepoint_container &spacepoints_per_event) {
    traccc::spacepoint_writer spwriter{std::string("event") + std::to_string(event) + "-spacepoints.csv"};
    for (size_t i = 0; i < spacepoints_per_event.items.size(); ++i) {
        auto spacepoints_per_module = spacepoints_per_event.items[i];
        auto module = spacepoints_per_event.headers[i];

        for (const auto &spacepoint : spacepoints_per_module) {
            const auto &pos = spacepoint.global;
            spwriter.append({module, pos[0], pos[1], pos[2], 0., 0., 0.});
        }
    }
}

void write(int events, Results *result) {

#pragma omp parallel for
    for (int event = 0; event < events; ++event) {
        write_measurements(event, result->measurements[event]);
        write_spacepoints(event, result->spacepoints[event]);
    }
}

auto *run(int events, vecmem::vector<traccc::host_cell_container> input_data) {

    // Algorithms
    traccc::component_connection cc;
    traccc::measurement_creation mt;
    traccc::spacepoint_formation sp;

    auto startAlgorithms = std::chrono::system_clock::now();

    // Output stats
    int64_t n_modules = 0;
    uint64_t n_cells = 0;
    uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_space_points = 0;

    Results *results = new Results(events);

#pragma omp parallel for reduction(+:n_modules, n_cells, n_clusters, n_measurements, n_space_points)
    for (int event = 0; event < events; ++event) {
        traccc::host_cell_container cells_per_event = input_data[event];
        // Output containers
        traccc::host_measurement_container measurements_per_event;
        traccc::host_spacepoint_container spacepoints_per_event;
        measurements_per_event.headers.reserve(cells_per_event.headers.size());
        measurements_per_event.items.reserve(cells_per_event.headers.size());
        spacepoints_per_event.headers.reserve(cells_per_event.headers.size());
        spacepoints_per_event.items.reserve(cells_per_event.headers.size());

#pragma omp parallel for
        for (std::size_t i = 0; i < cells_per_event.items.size(); ++i) {
            auto &module = cells_per_event.headers[i];
            module.pixel = traccc::pixel_segmentation{-8.425, -36.025, 0.05, 0.05};

            // The algorithmic code part: start
            traccc::cluster_collection clusters_per_module = cc(cells_per_event.items[i], cells_per_event.headers[i]);
            clusters_per_module.position_from_cell = module.pixel;

            traccc::host_measurement_collection measurements_per_module = mt(clusters_per_module, module);
            traccc::host_spacepoint_collection spacepoints_per_module = sp(module, measurements_per_module);
            // The algorithmnic code part: end

            n_modules += 1;
            n_cells += cells_per_event.items[i].size();
            n_clusters += clusters_per_module.items.size();
            n_measurements += measurements_per_module.size();
            n_space_points += spacepoints_per_module.size();

#pragma omp critical
            {
                measurements_per_event.items.push_back(std::move(measurements_per_module));
                measurements_per_event.headers.push_back(module);

                spacepoints_per_event.items.push_back(std::move(spacepoints_per_module));
                spacepoints_per_event.headers.push_back(module.module);
            }
        }

#pragma omp critical
        {
            results->measurements.push_back(std::move(measurements_per_event));
            results->spacepoints.push_back(std::move(spacepoints_per_event));
        };

    }

    auto endAlgorithms = std::chrono::system_clock::now();
    std::chrono::duration<double> diffAlgo = endAlgorithms - startAlgorithms;
    std::cout << "Algorithms time: " << diffAlgo.count() << " sec." << std::endl;

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_cells << " cells from " << n_modules << " modules" << std::endl;
    std::cout << "- created " << n_clusters << " clusters. " << std::endl;
    std::cout << "- created " << n_measurements << " measurements. " << std::endl;
    std::cout << "- created " << n_space_points << " space points. " << std::endl;

    return results;
}


// The main routine
//
int main(int argc, char *argv[]) {

    if (argc < 4) {
        std::cout << "Not enough arguments, minimum requirement: " << std::endl;
        std::cout << "./funct_par_example <detector_file> <cell_directory> <events>" << std::endl;
        return -1;
    }

    auto detector_file = std::string(argv[1]);
    auto cell_directory = std::string(argv[2]);
    auto events = std::atoi(argv[3]);

    auto start = std::chrono::system_clock::now();

    write(events,
          run(events,
              read(events, detector_file, cell_directory)));

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Total execution time: " << diff.count() << " sec." << std::endl;
}

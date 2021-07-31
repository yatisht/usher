#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "../usher_graph.hpp"

namespace po = boost::program_options;

extern Timer timer;

po::variables_map check_options(po::parsed_options parsed);

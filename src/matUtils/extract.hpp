#include "common.hpp"
#include "convert.hpp"
#include "filter.hpp"
#include "describe.hpp"
#include "uncertainty.hpp"
#include "select.hpp"
#include "translate.hpp"

po::variables_map parse_extract_command(po::parsed_options parsed);
void extract_main (po::parsed_options parsed);
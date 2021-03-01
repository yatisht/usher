#include "annotate.hpp"
#include "mask.hpp"
#include "convert.hpp"
#include "filter.hpp"
#include "describe.hpp"
#include "uncertainty.hpp"
#include "select.hpp"
#include "summary.hpp"

Timer timer; 

int main (int argc, char** argv) {
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, mask, convert, prune and describe.")
        ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    
    try {
        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();

        po::store(parsed, vm);
        std::string cmd = vm["command"].as<std::string>();
        if (cmd == "annotate"){
            annotate_main(parsed);
        } else if (cmd == "convert"){
            convert_main(parsed);
        } else if (cmd == "mask"){
            mask_main(parsed); 
        } else if (cmd == "filter"){
            filter_main(parsed); 
        } else if (cmd == "describe"){
            describe_main(parsed); 
        } else if (cmd == "uncertainty") {
            uncertainty_main(parsed);
        } else if (cmd == "summary") {
            summary_main(parsed);
        } else if (cmd == "help" || cmd == "--help" || cmd == "-h") { 
            // TODO: improve this message
            fprintf(stderr, "matUtils has several major subcommands: annotate, mask, convert, filter, uncertainty, summary, and describe.\nIndividual command options can be accessed with matUtils command --help, e.g. matUtils annotate --help will show annotation-specific help messages.");
            exit(0);
        } else {
            fprintf(stderr, "Invalid command. Please choose from annotate, mask, convert, filter, describe, uncertainty, summary, or help and try again.\n"););
            exit(1);
        }
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Please choose from annotate, mask, convert, filter, describe, uncertainty, summary, or help and try again.\n");
        exit(0);
    }

    return 0;
}


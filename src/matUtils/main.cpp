#include "annotate.hpp"
#include "filter.hpp"
#include "convert.hpp"
#include "prune.hpp"
#include "describe.hpp"

Timer timer; 

int main (int argc, char** argv) {
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, filter, convert, prune and describe.")
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
        } else if (cmd == "filter"){
            filter_main(parsed); 
        } else if (cmd == "prune"){
            prune_main(parsed); 
        } else if (cmd == "describe"){
            describe_main(parsed); 
        } else if (cmd == "help" || cmd == "--help" || cmd == "-h") { 
            // TODO: improve this message
            fprintf(stderr, "matUtils has three major subcommands: annotate, filter, convert, prune and describe. All three take a MAT .pb file as input.\nAnnotate adds information to the MAT. Use annotate when you have information you want to calculate or incorporate into the MAT .pb. The command 'matUtils annotate --help' will describe related options.\nFilter removes nodes or samples from the MAT. Use filter when you want to strip out low-quality samples or mask samples you want to avoid. The command 'matUtils filter --help' will describe related options.\nConvert produces files which are not MAT .pb format, such as .vcf or newick text files. Use convert when you want other file types. The command 'matUtils convert --help' will describe related options.\n");
            exit(0);
        } else {
            fprintf(stderr, "Invalid command. Please choose from annotate, filter, convert, prune, describe or help and try again.\n");
            exit(1);
        }
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Please choose from annotate, filter, convert, prune, describe or help and try again.\n");
        exit(0);
    }

    return 0;
}


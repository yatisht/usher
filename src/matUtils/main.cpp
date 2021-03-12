#include "annotate.hpp"
#include "mask.hpp"
#include "summary.hpp"
#include "extract.hpp"

Timer timer; 

int main (int argc, char** argv) {
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, mask, extract, uncertainty, and summary.")
        ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    std::string cmd;
    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();
    //this help string shows up over and over, lets just define it once
    std::string helpstr = "matUtils has several valid subcommands: \n\n"
        "extract: subsets the input MAT on various conditions and/or converts to other formats (newick, VCF, JSON etc)\n\n"
        "summary: calculates basic statistics and counts members in the input MAT\n\n"
        "annotate: assigns clade identities to nodes, directly or by inference\n\n"
        "uncertainty: calculates sample placement uncertainty metrics and writes the results to tsv\n\n"
        "mask: masks the mutations specific to the restricted set of input samples\n\n"
        "Individual command options can be accessed with matUtils command --help, e.g. matUtils annotate --help will show annotation-specific help messages.\n\n";
    
    try {
        po::store(parsed, vm);
        cmd = vm["command"].as<std::string>();
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Help follows:\n\n");
        fprintf(stderr, "%s", helpstr.c_str());
        //0 when no command is selected because that's what passes tests.
        exit(0);
    }
    if (cmd == "extract") {
        extract_main(parsed);
    } else if (cmd == "annotate") {
        annotate_main(parsed);
    } else if (cmd == "mask"){
        mask_main(parsed); 
    } else if (cmd == "uncertainty") {
        uncertainty_main(parsed);
    } else if (cmd == "summary") {
        summary_main(parsed);
    } else if (cmd == "help") { 
        fprintf(stderr, "%s", helpstr.c_str());
        exit(0);
    } else {
        fprintf(stderr, "Invalid command. Help follows:\n\n");
        fprintf(stderr, "%s", helpstr.c_str());
        exit(1);
    }

    return 0;
}


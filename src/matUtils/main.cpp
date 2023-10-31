#include "annotate.hpp"
#include "mask.hpp"
#include "summary.hpp"
#include "extract.hpp"
#include "merge.hpp"
#include "introduce.hpp"
#include "fix.hpp"
#include "version.hpp"
#include <cstddef>

Timer timer;

int main (int argc, char** argv) {
    po::options_description global("Command options");
    global.add_options()
    ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, mask, extract, uncertainty, introduce, fix, merge, version, and summary.")
    ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    std::string cmd;
    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();
    //this help string shows up over and over, lets just define it once
    std::string cnames[] = {"COMMAND","summary","extract","annotate","uncertainty","introduce", "merge", "mask", "fix", "version"};
    std::string chelp[] = {
        "DESCRIPTION\n\n",
        "calculates basic statistics and counts samples, mutations, and clades in the input MAT\n\n",
        "subsets the input MAT on various conditions and/or converts to other formats (newick, VCF, JSON)\n\n",
        "assigns clade identities to nodes, directly or by inference\n\n",
        "calculates sample placement uncertainty metrics and writes the results to tsv\n\n",
        "given sample region information, heuristically identifies points of geographic introduction along the phylogeny\n\n",
        "merge all samples of two input MAT files into a single output MAT \n\n",
        "masks the input samples\n\n",
        "fixes grandparent-reversion structures\n\n",
        "display version number\n\n"
    };
    try {
        po::store(parsed, vm);
        cmd = vm["command"].as<std::string>();
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "\nNo command selected. Help follows:\n\n");
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        //0 when no command is selected because that's what passes tests.
        exit(0);
    }
    if (cmd == "extract") {
        extract_main(parsed);
    } else if (cmd == "annotate") {
        annotate_main(parsed);
    } else if (cmd == "uncertainty") {
        uncertainty_main(parsed);
    } else if (cmd == "merge") {
        merge_main(parsed);
    } else if (cmd == "summary") {
        summary_main(parsed);
    } else if (cmd == "introduce") {
        introduce_main(parsed);
    } else if (cmd == "mask") {
        mask_main(parsed);
    } else if (cmd == "fix") {
        fix_main(parsed);
    } else if (cmd == "version") {
        std::cerr << "matUtils (v" << PROJECT_VERSION << ")" << std::endl;
    } else if (cmd == "help") {
        fprintf(stderr, "\n");
        std::cerr << "matUtils (v" << PROJECT_VERSION << ")" << std::endl;
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        exit(0);
    } else {
        fprintf(stderr, "\nInvalid command. Help follows:\n\n");
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        exit(1);
    }

    return 0;
}


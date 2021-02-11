#include "annotate.hpp"
#include "filter.hpp"
#include "convert.hpp"

Timer timer; 

int main (int argc, char** argv) {
    /*
    The new design principle for organizing matUtils is to divide the overall structure into three options. All three take protobuf files as input.
    First there is annotate, which calculates and stores uncertainty metrics and other tree-based metadata into the protobuf file. This returns a protobuf that is larger than the input (in bytes).
    Second there is filter, which prunes the tree based on threshold arguments and the uncertainty metrics from annotate. This returns a protobuf that is smaller than the input (in bytes).
    Third there is convert, which produces a different file format than protobuf. This returns non-protobuf files.

    Generally a workflow will call annotate first, then filter, then convert. Or, if metadata is precalculated and saved on our updated global tree, they should be able to call filter then convert without an extended annotate command.
    Ideally these would be chainable on the command line, though I'm not sure if we can write a protobuf file to stdout/read from stdin. This is something we may want to investigate as a future option.
    For now, they can use && and intermediate file names to write a simple matUtils pipeline.
    */
    po::options_description global("Command options");
    global.add_options()
        ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, filter, convert.")
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
        } else if (cmd == "help" || cmd == "--help" || cmd == "-h") { //trying to catch some of the intuitive things people will try.
            fprintf(stderr, "matUtils has three major subcommands: annotate, filter, and convert. All three take a MAT .pb file as input.\nAnnotate adds information to the MAT. Use annotate when you have information you want to calculate or incorporate into the MAT .pb. The command 'matUtils annotate --help' will describe related options.\nFilter removes nodes or samples from the MAT. Use filter when you want to strip out low-quality samples or mask samples you want to avoid. The command 'matUtils filter --help' will describe related options.\nConvert produces files which are not MAT .pb format, such as .vcf or newick text files. Use convert when you want other file types. The command 'matUtils convert --help' will describe related options.\n"); //very open to alternative wording/formatting to make this nicer/clearer.
            exit(0);
        } else {
            fprintf(stderr, "Invalid command. Please choose from annotate, filter, convert, or help and try again.\n");
            exit(1);
        }
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "No command selected. Please choose from annotate, filter, convert, or help and try again.\n");
        exit(0);
    }

    return 0;
}

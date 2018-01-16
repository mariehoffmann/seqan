


#include <iostream>
#include <map>

#include <boost/filesystem.hpp>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

//#include <seqan/parallel.h>
//#include <seqan/index.h>
//#include <seqan/store.h>
//#include <seqan/seq_io.h>
//#include <seqan/version.h>

#define MS_BINARY_NAME "matching_statistics"

namespace fs = boost::filesystem;

struct MSOptions
{
    int numThreads;
    int verbosity;
    CharString inputFilename;
    CharString outputFilename;
    std::string tmp_dir;
    fs::path tmp_dir_path;  // full canonical path
	MSOptions()
	{
        numThreads = 1;
        verbosity = 0;
        tmp_dir = "./tmp/";
	}
};

// Parse the command line and return the status of the parsing.
seqan::ArgumentParser::ParseResult
parseCommandLine(MSOptions & options, int argc, char const ** argv)
{
    // Setup command line parser.
    seqan::ArgumentParser parser(MS_BINARY_NAME);

    // Set short description, version, and date.
    setShortDescription(parser, "Matching statistics");
    //setCategory(parser, "Indices");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-g\\fP \\fIIN.{fq,fa}\\fP \\fIOUT.{fq,fa}\\fP \\fITMP_DIR\\fP");
    addDescription(parser,
                   "matching_statistics is a tool to compute longest or shortest "
                   "repeated substrings.");
    addDescription(parser,
                   "You have to specify the fasta or fastq input file \\fIIN.{fq,fa}\\fP. "
                   "Output will be written to \\fIOUT.{fq,fa}\\fP."
                   "Temporary files are written to \\fITMP_DIR\\fP if provided or to ./tmp/.");

    // MS gets two parameters:  The paths to the input and the output files.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, 0, "fa fasta fq fastq");
    setHelpText(parser, 0, "An input file with reads to be processed.");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, 1, "fa fasta fq fastq");
    setHelpText(parser, 1, "An output file to store the corrected reads.");

    // General Options.
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose output."));
    addOption(parser, seqan::ArgParseOption("t", "tmp_dir", "Folder for temporary output."));

    // Usage Examples.
    addTextSection(parser, "Examples");
    addText(parser,
            "Most users will only have to specify the genome length using the mandatory \\fB-g\\fP parameter and "
            "enable multi-threading using the \\fB-nt\\fP option.  For best performance, use as many threads "
            "as you have (virtual) cores in your machine.");

    std::string toolName = "\\fB" MS_BINARY_NAME "\\fP";
    addListItem(parser, toolName + " IN.fq OUT.fq",
                "Process reads in \\fIIN.fq\\fP with one thread and write the results to \\fIOUT.fq\\fP.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine().
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    if (isSet(parser, "verbose"))
        options.verbosity = 1;

    if (options.verbosity >= 1)
        std::cerr << "MS - Matching Statistics\n"
                  << "=======================\n\n";

    // Extract Argument and Option Value
    getArgumentValue(options.inputFilename, parser, 0);
    getArgumentValue(options.outputFilename, parser, 1);
    if (options.verbosity >= 1)
    {
        std::cerr << "Loading reads from " << options.inputFilename << std::endl;
        std::cerr << "Writing statistics to " << options.outputFilename << std::endl;
    }
    // TODO: when calling with -t my_tmp_dir there is error too many arguments
    getOptionValue(options.tmp_dir, parser, "tmp_dir");
    // Create tmp directory if not existing
    std::cout << "tmp_dir = " << options.tmp_dir << std::endl;
    fs::path full_path = fs::path(options.tmp_dir);
    std::cout << "full_path: " << full_path.string() << std::endl;
    std::cout << fs::exists(full_path) << std::endl;
    if (!fs::exists(full_path))
    {
        if (boost::filesystem::create_directory(dir))
            if (options.verbosity >= 1)
                std::cout << "Create tempory directory: " << full_path.string() << "\n";
        else
        {
            if (options.verbosity >= 1)
                std::cerr << "Failed to create temporary directory: " << full_path.string() << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }

    }
    //full_path = fs::absolute(full_path);
    //std::cout << "full_path: " << full_path.string() << std::endl;

    /*bool dir_exists = fs::exists(full_path);
    std::cout << "dir exists: " << dir_exists << std::endl;
    bool dir_is = fs::is_directory(full_path);
    std::cout << "dir is: " << dir_is << std::endl;
*/
//    fs::path canonicalPath = fs::canonical(full_path, relativeTo);
//    std::cout << "full_path: " << full_path.string() << std::endl;


    return seqan::ArgumentParser::PARSE_OK;
}

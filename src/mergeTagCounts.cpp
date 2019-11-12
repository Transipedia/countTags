#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "optionparser.h"
#include "dna.h"
#include "utils.h"
#include "CountsTable.h"

using namespace std;

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};
enum  optionIndex {UNKNOWN,HELP,PROBE_LENGTH,MIN_RECURRENCY};
const option::Descriptor usage[] =
{
  {UNKNOWN, 0,"" , ""    ,
    option::Arg::None, "USAGE: mergeTagCounts [options] tag-counts1.tsv tags-counts2.tsv\n\nOptions:" },
  {HELP,    0,"" , "help",
    option::Arg::None, "  --help  \tPrint usage and exit." },
  {PROBE_LENGTH, 0, "k","",
    Arg::Numeric,      "  -k=INT      \ttags length" },
  {MIN_RECURRENCY, 0, "m","",
    Arg::Numeric,      "  -m=INT      \tmin recurrency" },
  {0,0,0,0,0,0}
};

int main (int argc, char *argv[]) {
  uint tag_length = 30;
  uint min_recurrency = 1;
  uint nb_samples;
  std::string delimiter = "\t";
  CountsTable * tag_counts;

  /**********************************
   *
   *        Parsing options
   *
   *********************************/
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  option::Option options[stats.options_max], buffer[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error())
    return 1;

  if (options[HELP] || argc == 0) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  if (options[PROBE_LENGTH]) {
    tag_length = atoi(options[PROBE_LENGTH].arg);
  }

  if (options[MIN_RECURRENCY]) {
    min_recurrency = atoi(options[MIN_RECURRENCY].arg);
  }

  nb_samples = parse.nonOptionsCount();

  tag_counts = new CountsTable(nb_samples,tag_length);

  for (int s = 0; s < nb_samples; ++s) {
    cerr << "Reading file: " << parse.nonOption(s) << endl;
    string line;
    ifstream counts_file (parse.nonOption(s));
    if (counts_file.is_open())
    {
      // Reading the header line
      getline(counts_file,line);
      uint n = line.find(delimiter);
      string sample_name = line.substr(n+1);
      tag_counts->setSampleName(s,sample_name.c_str());

      while ( getline (counts_file,line) )
      {
        uint n = line.find(delimiter);
        string tag = line.substr(0,n);
        // BUG: segfault if name if present = text
        uint64_t count = stoi(line.substr(n+1));
        tag_counts->setCount(tag.c_str(),s,count);
      }
      counts_file.close();
    }
  }

  /* Filter-out counts */
  tag_counts->recurrencyFilter(min_recurrency);

  /* Print Counts */
  tag_counts->printCounts();

  delete tag_counts;
}

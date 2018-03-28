#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iterator>
#include <string.h>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

#include "optionparser.h"
#include "dna.h"
//#include <zlib.h>

#define MILLION 1000000
#define VERS "0.3"

//return the minumum value of the k-mer at pos p between strand rev and stran fwd
//TODO add a function that get a DNA string and a k value, and return a array of vector values
inline uint64_t valns(uint32_t p, char *dna,uint32_t k,int64_t *last,uint64_t *valfwd,uint64_t *valrev, bool stranded = false){
  int e=p-*last;
  if(e!=1){
    *last=p;
    *valfwd=DNAtoInt(&dna[p], k, true);
    *valrev=intRevcomp(*valfwd, k);
  }else{
    // Compute the new value from the previous one.
    uint64_t m=1;
    *valfwd%=m<<(2*k-2);
    *valfwd<<=2;
    int new_nuc = convNuc(dna[p+k-1]);
    *valfwd += new_nuc;
    *last=p;
    *valrev/=1<<(2);
    *valrev+=(uint64_t)compNuc(new_nuc)<<(2*k-2);
  }
  if(stranded || *valfwd < *valrev) {
    return *valfwd;
  } else {
    return *valrev;
  }
}

std::string join( const std::vector<std::string>& elements, const char* const separator)
{
  switch (elements.size())
  {
    case 0:
      return "";
    case 1:
      return elements[0];
    default:
      std::ostringstream os;
      std::copy(elements.begin(), elements.end()-1, std::ostream_iterator<std::string>(os, separator));
      os << *elements.rbegin();
      return os.str();
  }
}

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
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

enum  optionIndex {UNKNOWN,HELP,VERBOSE,VERSION,KMER_LENGTH,TAG_FILE,STRANDED,MAX_READS,NB_THREADS,NORMALIZE,TAG_NAMES,MERGE_COUNTS};
const option::Descriptor usage[] =
{
  {UNKNOWN,      0, "" , "",
    option::Arg::None, "The purpose of countTags is to count occurences of few tags in large set of fastq files.\n\n"
                       "USAGE: countTags [options] -i tags.file seq.fastq[.gz]\n"
                       "======\n"
                       " * Tags file format: fasta, tsv (tag[ \\t]name) or raw (tag)\n"
                       " * Use '-' for reading tags from STDIN\n"
                       "\nOptions:\n"
                       "========" },
  {TAG_FILE, 0, "i","",
    Arg::Required,     "  -i Tag_FileName      \ttag filename, or '-' for STDIN (MANDATORY)." },
  {HELP,         0, "h", "help",
    option::Arg::None, "  -h|--help  \tPrint usage and exit." },
  {VERBOSE,      0, "v", "verbose",
    option::Arg::None, "  -v|--verbose  \tPrint progress status on STDERR\n"
                       "  -vvv          \tPrint debug informations on STDERR." },
  {VERSION,      0, "V", "version",
    option::Arg::None, "  -V|--version  \tPrint version and exit." },
  {KMER_LENGTH, 0, "k","",
    Arg::Numeric,      "  -k INT      \ttag length [default: 22]." },
  {MAX_READS,    0, "m", "",
    Arg::Numeric,      "  -m INT      \tmax number of reads [default: UINT32_MAX]." },
  //{NB_THREADS, 0, "t","",
  //  Arg::Numeric,      "  -t=INT      \tnumber of threads" },
  {STRANDED,     0, "" , "stranded",
    option::Arg::None, "  --stranded  \tstrand-specific protocol." },
  {NORMALIZE,    0, "" , "normalize",
    option::Arg::None, "  --normalize  \tnormalize counts." },
  {TAG_NAMES,    0, "t", "tag-names",
    option::Arg::None, "  -t|--tag-names  \tprint tag names in the output." },
  {MERGE_COUNTS,    0,"" , "merge-counts",
    option::Arg::None, "  --merge-counts  \tmerge counts from all input FASTQs" },
  {UNKNOWN,      0, "" , "",
    option::Arg::None, "\nExamples:\n"
                       "=========\n"
                       " * countTags -k 30 --stranded -t MyBestTags.tsv MyAllFastq*.gz > MyCount.tsv\n"
                       " * countTags -k 30 -t - MyAllFastq*.gz < MyBestTags.raw\n" },
  {0,0,0,0,0,0}
};

class Tag {
public:
  uint64_t tag;
};

int main (int argc, char *argv[]) {

  // Config vars
  const char * tags_file;
  const char * seq_file;
  uint32_t tag_length = 22;
  bool stranded = false;
  bool normalize = false;
  bool print_tag_names = false;
  bool merge_counts = false;
  uint32_t nb_tags = 0;
  uint32_t max_reads = UINT32_MAX;
  int nb_threads = 1;
  uint32_t nb_samples;
  uint32_t i;
  uint32_t read_id;

  char * seq;
  std::string tag_name;
  uint32_t seq_length;
  uint64_t tag;
  uint64_t valrev,valfwd;
  uint64_t nb_factors;
  int64_t last;
  std::unordered_map<uint64_t,double*>  tags_counts;
  std::unordered_map<uint64_t,std::vector<std::string>> tags_names;
  std::unordered_map<uint64_t,double*>::iterator it_counts;
  std::vector<uint64_t> nb_factors_by_sample;

  // File vars
  std::string gzip_pipe = "gunzip -fc ";
  std::string tmp;
  FILE * file;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  uint32_t line_id = 0;

  /**********************************
   *
   *        Parsing options
   *
   *********************************/
  int verbose = 0; // verbose level

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

  if (options[VERSION]) {
    std::cout << VERS;
    return 0;
  }

  if (options[VERBOSE]) {
    verbose = options[VERBOSE].count();
  }

  // TODO we should check the length to choose the appropriate
  // integer size to use
  if (options[KMER_LENGTH]) {
    tag_length = atoi(options[KMER_LENGTH].arg);
    if (tag_length > 32) {
      std::cout << "ERROR: For now, K-mer length has to be < 32" << "\n\n";
      option::printUsage(std::cout, usage);
      return 1;
    }
  }

  if (options[MAX_READS]) {
    max_reads = atoi(options[MAX_READS].arg);
  }

  //if (options[NB_THREADS]) {
  //  nb_threads = atoi(options[NB_THREADS].arg);
  //}

  if (options[STRANDED]) {
    stranded = true;
  }

  if (options[NORMALIZE]) {
    normalize = true;
  }

  if (options[TAG_NAMES]) {
    print_tag_names = true;
  }

  if (options[MERGE_COUNTS]) {
    merge_counts = true;
  }

  if (options[UNKNOWN]) {
    for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
        std::cout << "Unknown option: " << opt->name << "\n";
    option::printUsage(std::cout, usage);
    return 1;
  }

  if(parse.nonOptionsCount() < 1) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  tags_file = options[TAG_FILE].arg;
  nb_samples = parse.nonOptionsCount();

  if (verbose>2) {
    std::cout << "File to analyse: " << std::to_string(parse.nonOptionsCount()) << std::endl;
    for (int i = 0; i < parse.nonOptionsCount(); ++i)
      fprintf(stdout, "Non-option argument #%d is %s\n", i, parse.nonOption(i));
  }

  /**********************************
   *
   *   Create tags counting table
   *
   *********************************/
  line_id = 0;

  if (verbose)
    std::cerr << "Counting k-mers" << std::endl;
  // Create hash table of k-mer counts
  line_id = 0;
  nb_tags = 0;

  std::istream *filein;
  std::ifstream fileif;
  if (strlen(tags_file) == 1) {
    // use stdin
    filein = &std::cin;
  } else {
    fileif.open(tags_file, std::ifstream::in);
    // check if not read error
    if (fileif.fail()) {
      std::cerr << "Error: Can't read "<< tags_file << std::endl;
      return 1;
    }
    filein = &fileif;
  }
  bool no_name = 1;
  // Parse file and detect file format (fas, raw or tsv)
  for (std::string lines; std::getline(*filein, lines); ) {
    if (lines.find(">") != std::string::npos) {
      // we got a fasta line
      no_name = 0;
      tag_name = lines;
      tag_name.erase(0,1); // remove the fasta ">" prefix
      if (verbose>2)
        std::cerr << "Find fasta line: " << tag_name << std::endl;
    } else {
      // Check if we have a raw or tsv line: tag\tname
      std::size_t found = lines.find_first_of(" \t");
      if (found != std::string::npos) {
        no_name = 0;
        tag_name = lines.substr(found+1);
        lines.erase(found, lines.length());
      }
      // insert nb_tags as tag_name if none read
      if (no_name)
        tag_name = std::to_string(nb_tags+1);
      // Take at least tag_length for each tags
      if (lines.length() < tag_length)
        continue;
      // convert tag to Int
      tag = DNAtoInt(lines.c_str(),tag_length,stranded);
      if (verbose>2)
        std::cerr << "tag: " << lines << ", name:" << tag_name << ", tagInt: " << tag;
      tags_counts[tag] = new double[nb_samples]();
      tags_names[tag].push_back(tag_name);
      nb_tags++;
      if (verbose>2)
       std::cerr << ", nb_tag: " << nb_tags << std::endl;
    }
    line_id++;
  }

   if (verbose)
     std::cerr << "Finished indexing tags" << std::endl;


  /**********************************
   *
   *            First pass
   *
   *********************************/
//#pragma omp parallel num_threads(nb_threads)
  for (int s = 0; s < nb_samples; ++s) {
    if (verbose)
       std::cerr << "Counting tags for file: " << parse.nonOption(s) << "\n";

    line = NULL;
    len = 0;
    line_id = 0;
    tmp = "";
//    file = popen(tmp.append(gzip_pipe).append(parse.nonOption(s)).c_str(),"r");
    nb_factors = 0;

    std::ifstream filein(parse.nonOption(s), std::ios_base::in | std::ios_base::binary);
    try {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(filein);
        for(std::string seq; std::getline(in, seq); )
        {
            // If this line is a sequence
            if(line_id % 4 == 1) {
              read_id = ((int)((double)line_id*0.25) + 1);
              if(read_id >= max_reads) {
                break;
              }
              // Print a user-friendly output on STDERR every each XXXX reads processed
              if (verbose && read_id % MILLION == 0) {
                std::cerr << (int)((double)line_id*0.25) + 1 << " reads parsed" << std::endl;
              }

              // Skip the sequence if the read length is < to the tag_length
              if(seq.length() < tag_length)
                continue;

              nb_tags = seq.length() - tag_length + 1;
              nb_factors += nb_tags;

              //uint64_t valrev,valfwd;
              last = -3;

              for(i = 0; i < nb_tags; i++) {
                it_counts = tags_counts.find(valns(i,(char *)seq.c_str(),tag_length,&last,&valfwd,&valrev,stranded));
                if(it_counts != tags_counts.end()) {
                  it_counts->second[s]++;
                }
              }
            }
            line_id++;
        }
    }
    catch(const boost::iostreams::gzip_error& e) {
         std::cout << e.what() << '\n';
    }

    nb_factors_by_sample.push_back(nb_factors);

    if(normalize && nb_factors > 0) {
      if (verbose)
        std::cerr << "Normalize counts" << std::endl;
      for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
        // TODO We should take into accout the error rate...
        if(it_counts->second[s] > 0)
          it_counts->second[s] = it_counts->second[s] * MILLION / nb_factors;
      }
    }

    // Close file and clear line buffer
    fclose(file);
    if (line)
      free(line);
  }

  /****
   * PRINT THE RESULTS
   */
  if (verbose)
    std::cerr << "tag_length: " << tag_length << "\n";
  // First print headers
  std::cout << "tag";
  if (print_tag_names)
    std::cout << "\ttag_names";
  if(!merge_counts) {
    for (int s = 0; s < nb_samples; ++s) {
      std::cout << "\t" << parse.nonOption(s);
    }
  } else {
    std::cout << "\tcounts";
  }
  std::cout << "\n";
  char *tag_seq = new char[tag_length+1];
  tag_seq[tag_length] = '\0';
  for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
    intToDNA(it_counts->first,tag_length,tag_seq);
    std::cout << tag_seq;
    if(print_tag_names) {
      std::cout << "\t" << join(tags_names[it_counts->first],",");
    }
    if(!merge_counts){
      for (int s = 0; s < nb_samples; ++s) {
        std::cout << "\t" << it_counts->second[s];
      }
    } else {
      double count_sum = 0;
      for (int s = 0; s < nb_samples; ++s) {
        count_sum += it_counts->second[s];
      }
      std::cout << "\t" << count_sum;
    }
    std::cout << std::endl;
  }

  std::cout << "total_factors";
  if (print_tag_names)
    std::cout << "\t*";
  if(!merge_counts) {
    for (int s = 0; s < nb_samples; ++s) {
      std::cout << "\t" << nb_factors_by_sample[s];
    }
  } else {
    uint64_t nb_factors_sum = 0;
    for (int s = 0; s < nb_samples; ++s) {
      nb_factors_sum += nb_factors_by_sample[s];
    }
    std::cout << "\t" << nb_factors_sum;
  }
  std::cout << std::endl;

}

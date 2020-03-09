/******************************************************************************\
*                                                                              *
*  Copyright © 2015-2018 -- IRMB/INSERM                                        *
*                           (Institute for Regenerative Medicine & Biotherapy  *
*                           Institut National de la Santé et de la Recherche   *
*                           Médicale)                                          *
*                                                                              *
*  Copyright (C) 2015-2018  Jérôme Audoux                                      *
*  Copyright (C) 2018-      Anthony Boureux                                    *
*                                                                              *
*                                                                              *
*  This file is part of countTags program.                                     *
*                                                                              *
*  The purpose of countTags is to count occurences of few tags in large set of *
*  fastq files.                                                                *
*                                                                              *
*   countTags is free software: you can redistribute it and/or modify          *
*   it under the terms of the GNU General Public License as published by       *
*   the Free Software Foundation, either version 3 of the License, or          *
*   (at your option) any later version.                                        *
*                                                                              *
*   countTags is distributed in the hope that it will be useful,               *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
*   GNU General Public License for more details.                               *
*                                                                              *
*   You should have received a copy of the GNU General Public License          *
*   along with countTags program. If not, see <http://www.gnu.org/licenses/>.  *
*                                                                              *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <string.h>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

// local .h
#include "optionparser.h"
#include "dna.h"
#include "version.h"
//#include <zlib.h>

#define MILLION 1000000
#define BILLION 1000000000

//return the minumum value of the k-mer at pos p between strand rev and stran fwd
//TODO add a function that get a DNA string and a k value, and return a array of vector values
inline uint64_t valns(uint32_t p, char *dna,uint32_t k,int64_t *last,uint64_t *valfwd,uint64_t *valrev, bool isstranded = false, bool getrev = false){
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
  // when paired and read are reverse
  if (getrev)
    return *valrev;
  // otherwise
  if(isstranded || *valfwd < *valrev) {
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
  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
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

enum  optionIndex {
  KMER_LENGTH,  // kmer length to use
  TAG_FILE,     // file with the tag
  TAG_NAMES,    // output the tag name
  MAX_READS,    // count on max reads
  STRANDED,     // count in stranded mode
  NOSTRANDED,   // count in no-stranded mode, usefull instead of having an empty option
  PAIRED,       // count in paired mode
  SUMMARY,      // output summary count in a file instead of console
  NORMALIZE,    // normalize count on kmer factor
  BILLIONOPT,      // normalize count on kmer factor by billion instead of million
  MERGE_COUNTS, // sum all column into one
  MERGE_COUNTS_COLNAME, // give a name for the merge column instead of 'count'
  READS_WRFILE, // output reads with tag in a file
  NB_THREADS,   // number of threads to use (TODO)
  UNKNOWN,HELP,VERBOSE,VERSIONOPT
  };
const option::Descriptor usage[] =
{
 /* const option::Descriptor usage[] = {
 *   { CREATE,                                            // index
 *     OTHER,                                             // type
 *     "c",                                               // shortopt
 *     "create",                                          // longopt
 *     Arg::None,                                         // check_arg
 *     "--create  Tells the program to create something." // help
 * }
 */
  {UNKNOWN,      0, "" , "",
    option::Arg::None, "The purpose of countTags is to count occurences of few tags in large set of fastq files.\n\n"
                       "USAGE: countTags [options] -i tags.file seq.fastq[.gz] ...\n"
                       "\nVersion: " VERSION "\n"
                       "======\n"
                       " * Tags file format: fasta, tsv (tag[ \\t]name) or raw (tag)\n"
                       " * Use '-' for reading tags from STDIN\n"
                       " * for now, countTags can't read tag file in gzip format, uncompress before and pass to countTags via a pipe (see example hereafter)\n"
                       "\nArguments:\n"
                       "========\n" },
  {UNKNOWN,      0, "", "",
    option::Arg::None, "Mandatory:" },
  {TAG_FILE, 0, "i","",
    Arg::NonEmpty,     "  -i Tag_FileName      \ttag filename, or '-' for STDIN." },
  {UNKNOWN,      0, "", "",
    option::Arg::None, "Options:" },
  {KMER_LENGTH, 0, "k","",
    Arg::Numeric,      "  -k INT      \ttag length [default: 22]." },
  {MAX_READS,    0, "", "maxreads",
    Arg::Numeric,      "  --maxreads INT      \tmax number of reads to analyze [default: INT32_MAX]." },
  //{NB_THREADS, 0, "t","",
  //  Arg::Numeric,      "  -t=INT      \tnumber of threads" },
  {STRANDED,     0, "" , "stranded",
    option::Arg::None, "  --stranded  \tanalyse only the strand of the read and tag (no reverse-complement)." },
  {NOSTRANDED,     0, "" , "nostranded",
    option::Arg::None, "  --nostranded  \tturn off stranded mode, do not care about strand." },
  {PAIRED,       0, "" , "paired",
    Arg::NonEmpty, "  --paired rf|fr|ff \tstrand-specific protocol (can use only 2 fastq with _1.fastq and _2.fastq in filename)." },
  {NORMALIZE,    0, "n" , "normalize",
    Arg::None, "  -n|--normalize  \tnormalize count on total of million of kmer present in each sample." },
  {BILLIONOPT,    0, "b" , "kbpnormalize",
    Arg::None, "  -b|--kbp \tnormalize count by billion of kmer instead of million (impliy -n|--normalize)." },
  {TAG_NAMES,    0, "t", "tag-names",
    Arg::None, "  -t|--tag-names  \tprint tag names in the output." },
  {MERGE_COUNTS, 0,"" , "merge-counts",
    Arg::None, "  --merge-counts  \tmerge counts from all input FASTQs" },
  {MERGE_COUNTS_COLNAME,    0,"" , "merge-counts-colname",
    Arg::NonEmpty, "  --merge-counts-colname  \tcolumn name when merge counts is used" },
  {READS_WRFILE, 0, "r","reads",
    Arg::NonEmpty,     "  -r|--reads fileName      \twrite reads matching kmer in a fileName (for now store tag and read only, it's not a fastq file" },
  {SUMMARY,    0,"" , "summary",
    Arg::NonEmpty, "  --summary file    \tprint statistic in a file" },
  {VERBOSE,      0, "v", "verbose",
    option::Arg::None, "  -v|--verbose  \tPrint statistic on STDERR\n"
                       "  -vv           \tPrint progress status on STDERR\n"
                       "  -vvv          \tPrint debug informations on STDERR." },
  {VERSIONOPT,      0, "V", "version",
    option::Arg::None, "  -V|--version  \tPrint version and exit." },
  {HELP,         0, "h", "help",
    option::Arg::None, "  -h|--help  \tPrint usage and exit." },
  {UNKNOWN,      0, "" , "",
    option::Arg::None, "\nExamples:\n"
                       "=========\n"
                       " * countTags -k 30 --stranded -t -i MyBestTags.tsv MyAllFastq*.gz > MyCount.tsv\n"
                       " * countTags -k 30 -i MyBestTags.tsv --paired rf  MyAllFastq_1.fastq.gz MyAllFastq_2.fastq.gz > MyCount.tsv\n"
                       " * countTags -k 30 -t -i - MyAllFastq*.gz < MyBestTags.raw\n"
                       " * zcat MyBestTags.raw.gz | countTags -k 30 -t -i - --summary mystats.summary MyAllFastq*.gz > MyCount_table.tsv\n"
                       },
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
  bool isstranded = false;
  bool ispaired = false;
  std::string paired;
  bool normalize = false;
  double normalize_factors; // normalize factor MILLION or BILLION
  bool print_tag_names = false;
  bool merge_counts = false;
  std::string merge_counts_colname = "counts";
  uint32_t nb_tags = 0;
  uint32_t max_reads = UINT32_MAX;
  int nb_threads = 1;
  uint32_t nb_samples;
  uint32_t i;
  uint32_t read_id;

  char * seq;
  std::string tag_name;
  std::string output_read; // store filename to output read matching kmer
  std::string summary_file; // store summary information in a filename
  uint32_t seq_length;
  uint64_t tag;
  uint64_t valrev,valfwd;
  uint64_t nb_factors;
  int64_t last;
  std::unordered_map<uint64_t,double*>  tags_counts;
  std::unordered_map<uint64_t,std::vector<std::string>> tags_names;
  std::unordered_map<uint64_t,double*>::iterator it_counts;
  std::vector<uint64_t> nb_factors_by_sample;
  std::vector<uint64_t> nb_reads_by_sample;

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

  if (options[VERSIONOPT]) {
    std::cout << VERSION;
    return 0;
  }

  if (options[VERBOSE]) {
    verbose = options[VERBOSE].count();
  }

  // Test if we have a file in input
  // Not using Arg::Required from option parser, because
  // is not working with option -V -h
  if (!options[TAG_FILE].count()) {
    std::cerr << "ERROR : -i tag_file required\n\n";
    option::printUsage(std::cerr, usage);
    return 1;
  }

  // TODO we should check the length to choose the appropriate
  // integer size to use
  if (options[KMER_LENGTH]) {
    tag_length = atoi(options[KMER_LENGTH].arg);
    if (tag_length > 32) {
      std::cerr << "ERROR: For now, K-mer length has to be < 32" << "\n\n";
      option::printUsage(std::cerr, usage);
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
    isstranded = true;
  }

  if (options[PAIRED]) {
    // turn ON paired option
    ispaired = true;
    isstranded = true;
    // turn ON strander mode, because it means nothing without it
    if (options[PAIRED].count()) {
      paired = options[PAIRED].arg;
    } else {
      paired = "rf";
    }
    if (verbose>2) {
      std::cerr << "\tPaired mode turn ON, with option " << paired << ".\n";
    }
  }

  if (options[NOSTRANDED]) {
    isstranded = false;
    ispaired = false;
  }

  if (options[NORMALIZE]) {
    normalize = true;
    normalize_factors = MILLION;
  }

  if (options[BILLIONOPT]) {
    normalize = true;
    normalize_factors = BILLION;
  }

  if (options[TAG_NAMES]) {
    print_tag_names = true;
  }

  if (options[MERGE_COUNTS]) {
    merge_counts = true;
    if (options[MERGE_COUNTS_COLNAME]) {
      if (verbose>2) {
        std::cerr << "\tColumn name when merging:" << options[MERGE_COUNTS_COLNAME].arg << "\n";
      }
      merge_counts_colname = options[MERGE_COUNTS_COLNAME].arg;
    }
  }

  if (options[READS_WRFILE].count()) {
    output_read = options[READS_WRFILE].arg;
  }

  if (options[SUMMARY].count()) {
    summary_file = options[SUMMARY].arg;
    // turn verbose to 1 at least
    verbose = verbose ? verbose++ : 1;
  }

  if (options[UNKNOWN]) {
    for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
        std::cerr << "Unknown option: " << opt->name << "\n";
    option::printUsage(std::cerr, usage);
    return 1;
  }

  if(parse.nonOptionsCount() < 1) {
    std::cerr << "No fastq file provided ?" << "\n";
    option::printUsage(std::cerr, usage);
    return 0;
  }

  tags_file = options[TAG_FILE].arg;
  nb_samples = parse.nonOptionsCount();

  if (verbose>2) {
    std::cerr <<  "Version: " << VERSION << std::endl;
    std::cerr << "File to analyse: " << std::to_string(parse.nonOptionsCount()) << std::endl;
    for (int i = 0; i < parse.nonOptionsCount(); ++i)
      fprintf(stderr, "Non-option argument #%d is %s\n", i, parse.nonOption(i));

  }

  /**********************************
   *
   *   Create tags counting table
   *
   *********************************/
  line_id = 0;

  if (verbose > 1)
    std::cerr << "Counting k-mers" << std::endl;
  // Create hash table of k-mer counts
  line_id = 0;
  nb_tags = 0;

  std::istream *filein;
  std::ifstream fileif;
  if (*tags_file == '-') {
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
      if (lines.length() < tag_length) {
        std::cerr << "Error: tag lower than kmer: " << lines << std::endl;
        continue;
      }
      // convert tag to Int
      tag = DNAtoInt(lines.c_str(),tag_length,isstranded);
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

  // Bug: didn't test if tag file is empty
  if (line_id == 0) {
    std::cerr << "I did not get or understand your tag sequences" << std::endl;
    exit(2);
  }

   if (verbose > 1)
     std::cerr << "Finished indexing tags" << std::endl;


  /**********************************
   *
   *            First pass
   *
   *********************************/

  // open file tp output reads matching kmer
  std::ofstream hfile_read;
  if (output_read.length()) {
    hfile_read.open(output_read, std::ifstream::out);
    // check if not write error
    if (hfile_read.fail()) {
      std::cerr << "Error: Can't write to read file "<< output_read << std::endl;
      return 1;
    }
  }

  // open file tp output summary information
  std::ofstream hfile_summary;
  if (summary_file.length()) {
    hfile_summary.open(summary_file, std::ifstream::out);
    // check if not write error
    if (hfile_summary.fail()) {
      std::cerr << "Error: Can't write to summary file "<< summary_file << std::endl;
      return 1;
    }
    // send cerr to hfile_summary;
    std::cerr.rdbuf(hfile_summary.rdbuf());
  }
  // print arguments use to summary file
  if (verbose) {
    std::cerr << "CountTags version\t" << VERSION << "\n";
    std::cerr << "Kmer_size\t" << tag_length << "\n";
    std::cerr << "Tag file in\t" << tags_file << "\n";
    if (max_reads < UINT32_MAX)
      std::cerr << "Maximun reads analyzed\t" << max_reads << "\n";
    std::cerr << "Normalize\t" << (normalize ? "Yes" : "No") << "\n";
    if (normalize) {
      std::cerr << "Normalize by " << normalize_factors << " of kmer." << "\n";
    }
    std::cerr << "Stranded\t" << (isstranded ? "Yes" : "No") << "\n";
    std::cerr << "Paired\t" << (ispaired ? paired : "No") << "\n";
    std::cerr << "Merge count\t" << (merge_counts ? "Yes" : "No") << "\n";
    std::cerr << "Write matched read in file\t" << (output_read.length() ? output_read : "None") << "\n";
  }
//#pragma omp parallel num_threads(nb_threads)
  for (int sample = 0; sample < nb_samples; ++sample) {
    if (verbose > 1)
       std::cerr << "Counting tags for file: " << "\t" << parse.nonOption(sample) << "\n";

    line = NULL;
    len = 0;
    line_id = 0;
    tmp = "";
    // Test if fastq file is present, otherwise exit with error 10
    std::ifstream testfile(parse.nonOption(sample));
    if (!testfile.good()) {
      std::cerr << "Error: Can't read fastq file " << parse.nonOption(sample) << std::endl;
      return 10;
    }

    // open file via pipe, so using another thread to gunzip the file
    file = popen(tmp.append(gzip_pipe).append(parse.nonOption(sample)).c_str(), "r");
    nb_factors = 0;

    // ispaired: have to get reverse complement for reverse pair
    // use getrev to get the reverse complement when rf/fr/ff
    bool getrev = false;

    if (ispaired) {
      if ( std::string(parse.nonOption(sample)).find("_1.fastq") != std::string::npos) {
        // we got the first pair
        if (paired.compare("rf") == 0) {
          getrev = true;
        }
      } else {
        // we got the second pair
        if (paired.compare("fr") == 0) {
          getrev = true;
        }
      }
    }
    if (verbose > 2)
      std::cerr << "Paired mode ON, getrev = " << std::to_string(getrev) << ", for file " << parse.nonOption(sample) << std::endl;

    while ((read = getline(&line, &len, file)) != -1) {
      // If this line is a sequence
      if(line_id % 4 == 1) {
        read_id = ((int)((double)line_id*0.25) + 1);
        if(read_id >= max_reads) {
          break;
        }
        // Print a user-friendly output on STDERR every each XXXX reads processed
        if (verbose > 1 && read_id % MILLION == 0) {
          std::cerr << (int)((double)line_id*0.25) + 1 << " reads parsed" << std::endl;
        }
        // set seq to line
        seq = line;
        seq_length = strlen(seq) - 1; // Minus 1 because we have a new line

        // Skip the sequence if the read length is < to the tag_length
        if(seq_length < tag_length)
          continue;

        nb_tags = seq_length - tag_length + 1;
        nb_factors += nb_tags;

        //uint64_t valrev,valfwd;
        last = -3;

        for(i = 0; i < nb_tags; i++) {
          it_counts = tags_counts.find(valns(i, seq, tag_length, &last, &valfwd, &valrev, isstranded, getrev));
          if(it_counts != tags_counts.end()) {
            it_counts->second[sample]++;
            // output read in a file if required
            // TODO: move new before loop, to avoid creation each time when option is on
            if (output_read.length()) {
              char *tag_seq = new char[tag_length+1];
              tag_seq[tag_length] = '\0';
                intToDNA(it_counts->first, tag_length,tag_seq);
                hfile_read << tag_seq;
                if(print_tag_names) {
                  hfile_read << "\t" << join(tags_names[it_counts->first], ",");
                }
              hfile_read << "\t" << parse.nonOption(sample) << "\t" << seq;
            }
          }
        }
      }
      line_id++;
    }

    // store statistic
    nb_factors_by_sample.push_back(nb_factors);
    nb_reads_by_sample.push_back(line_id);

//    if(normalize && nb_factors > 0) {
//      if (verbose > 1 )
//        std::cerr << "Normalize counts" << std::endl;
//      for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
//        // TODO We should take into accout the error rate...
//        if(it_counts->second[sample] > 0) {
//          it_counts->second[sample] = it_counts->second[sample] * normalize_factors / nb_factors;
//        }
//      }
//    }

    // Close file and clear line buffer
    pclose(file);
    if (line)
      free(line);
  }

  /****
   * PRINT THE RESULTS
   */
  // First print headers
  std::cout << "tag";
  if (print_tag_names)
    std::cout << "\ttag_names";
  if(!merge_counts) {
    for (int sample = 0; sample < nb_samples; ++sample) {
      std::cout << "\t" << parse.nonOption(sample);
    }
  } else {
    std::cout << "\t" << merge_counts_colname;
  }
  std::cout << "\n";
  // print tag + value
  char *tag_seq = new char[tag_length+1];
  tag_seq[tag_length] = '\0';
  for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
    intToDNA(it_counts->first,tag_length,tag_seq);
    std::cout << tag_seq;
    if(print_tag_names) {
      std::cout << "\t" << join(tags_names[it_counts->first],",");
    }
    if(!merge_counts){
      for (int sample = 0; sample < nb_samples; ++sample) {
        std::cout << "\t" << it_counts->second[sample];
      }
    } else {
      double count_sum = 0;
      for (int sample = 0; sample < nb_samples; ++sample) {
        count_sum += it_counts->second[sample];
      }
      std::cout << "\t" << count_sum;
    }
    std::cout << std::endl;
  }

  // print statistic
  std::cerr << "# Total statistic per file\n";
  std::cerr << "File\t";
  if(!merge_counts) {
    for (int sample = 0; sample < nb_samples; ++sample) {
      std::cerr << "\t" << parse.nonOption(sample);
    }
  } else {
    std::cerr << "\t" << merge_counts_colname;
  }
  std::cerr << "\n";
  std::cerr << "total_factors";
  if(!merge_counts) {
    for (int sample = 0; sample < nb_samples; ++sample) {
      std::cerr << "\t" << nb_factors_by_sample[sample];
    }
  } else {
    uint64_t nb_factors_sum = 0;
    for (int sample = 0; sample < nb_samples; ++sample) {
      nb_factors_sum += nb_factors_by_sample[sample];
    }
    std::cerr << "\t" << nb_factors_sum;
  }
  std::cerr << std::endl;
  std::cerr << "total_reads";
  if(!merge_counts) {
    for (int sample = 0; sample < nb_samples; ++sample) {
      std::cerr << "\t" << nb_reads_by_sample[sample];
    }
  } else {
    uint64_t sum = 0;
    for (int sample = 0; sample < nb_samples; ++sample) {
      sum += nb_reads_by_sample[sample];
    }
    std::cerr << "\t" << sum;
  }
  std::cerr << std::endl;

}

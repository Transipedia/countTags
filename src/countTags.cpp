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
#include <memory>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>

// use namespace std by default
using namespace std;

// local .h
#include "optionparser.h"
#include "dna.h"
#include "version.h"
#include "zstr.hpp"

#define MILLION 1000000
#define BILLION 1000000000

//return the minumum value of the k-mer at pos p between strand rev and stran fwd
//TODO add a function that get a DNA string and a k value, and return a array of vector values
inline uint64_t valns(uint32_t p, char *dna, uint32_t k, int64_t *last, uint64_t *valfwd, uint64_t *valrev, bool isstranded = false, bool getrev = false){
  int e=p-*last;
  if (e!=1){
    *last=p;
    *valfwd=DNAtoInt(&dna[p], k, true);
    *valrev=intRevcomp(*valfwd, k);
  } else{
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
  if (isstranded || *valfwd < *valrev) {
    return *valfwd;
  } else {
    return *valrev;
  }
}

string join( const vector<string>& elements, const char* const separator)
{
  switch (elements.size())
  {
    case 0:
      return "";
    case 1:
      return elements[0];
    default:
      ostringstream os;
      copy(elements.begin(), elements.end()-1, ostream_iterator<string>(os, separator));
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
  ALLTAGS,      // use all tags if sequence length > kmer length
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
                       "USAGE: countTags [options] -i tags.file[.gz] seq.fastq[.gz] ...\n"
                       "\nVersion: " VERSION "\n"
                       "======\n"
                       " * Tags file format: fasta, tsv (tag[ \\t,;]name) or raw (tag)\n"
                       " * Use '-' for reading tags from STDIN\n"
                       " * the file can be in gzip format, otherwise uncompress before and pass to countTags via a pipe (see example hereafter)\n"
                       "\nArguments:\n"
                       "========\n" },
  {UNKNOWN,      0, "", "",
    option::Arg::None, "Mandatory:" },
  {TAG_FILE, 0, "i","",
    Arg::NonEmpty,     "  -i Tag_FileName      \ttag filename, or '-' for STDIN ( gzip or not)." },
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
  {ALLTAGS,       0, "a" , "alltags",
    Arg::None, "  -a|--alltags \tGenerate all tags from the sequence if its length is greater than kmer length." },
  {NORMALIZE,    0, "n" , "normalize",
    Arg::None, "  -n|--normalize  \tnormalize count on total of million of kmer present in each sample." },
  {BILLIONOPT,    0, "b" , "kbpnormalize",
    Arg::None, "  -b|--kbp \tnormalize count by billion of kmer instead of million (imply -n|--normalize)." },
  {TAG_NAMES,    0, "t", "tag-names",
    Arg::None, "  -t|--tag-names  \tprint tag names in the output." },
  {MERGE_COUNTS, 0,"" , "merge-counts",
    Arg::None, "  --merge-counts  \tmerge counts from all input FASTQs" },
  {MERGE_COUNTS_COLNAME,    0,"" , "merge-counts-colname",
    Arg::NonEmpty, "  --merge-counts-colname  \tcolumn name when merge counts is used (imply --merge-count)" },
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
                       " * countTags -k 30 -t -i MyBestTags.raw.gz --summary mystats.summary MyAllFastq*.gz > MyCount_table.tsv\n"
                       " * bzcat MyBestTags.raw.bz | countTags -k 30 -t -i - --summary mystats.summary MyAllFastq*.gz > MyCount_table.tsv\n"
                       },
  {0,0,0,0,0,0}
};

class Tag {
public:
  uint64_t tag;
};

int main (int argc, char *argv[]) {

  // Options vars
  int verbose = 0;                               // --verbose: verbose level
  uint32_t tag_length = 22;                      // -k: default kmer length
  uint32_t max_reads = UINT32_MAX;               // --maxreads: max number of reads to analyze
//int nb_threads = 1;                            // --threads: unused
  bool isstranded = false;                       // --stranded: is data stranded
  bool ispaired = false;                         // --paired: is data paired
  string paired = "rf";                          // store paired format ([rf], fr, ff)
  bool doalltags = false;                        // --alltags: generated all tags from sequence
  bool normalize = false;                        // --normalize: switch to normalize count
  double normalize_factors;                      // normalize factor MILLION or BILLION
  bool print_tag_names = false;                  // --tag-names: print tag name in sdtout
  bool merge_counts = false;                     // --merge-counts: merge all count in one column
  string merge_counts_colname = "counts";        // --merge-counts-colname: colname for merged column
  string output_read;                            // --reads: filename to output read matching kmer
  string summary_file;                           // --summary: store summary information in a filename
  uint32_t nb_samples = 0;                       // number of fastq file on argument line
  const char * tags_file;                        // -i: mandatory, tag filename

  /**********************************
   *
   *     Parsing options
   *
   *********************************/

  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  option::Option options[stats.options_max], buffer[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error())
    return 1;

  if (options[HELP] || argc == 0) {
    option::printUsage(cout, usage);
    return 0;
  }

  if (options[VERSIONOPT]) {
    cout << VERSION;
    return 0;
  }

  if (options[VERBOSE]) {
    verbose = options[VERBOSE].count();
  }

  // Test if we have a file in input
  // Not using Arg::Required from option parser, because
  // is not working with option -V -h
  if (!options[TAG_FILE].count()) {
    cerr << "ERROR : -i tag_file required\n\n";
    option::printUsage(cerr, usage);
    return 1;
  }

  // Check we can use the kmer size set by the user
  if (options[KMER_LENGTH]) {
    tag_length = atoi(options[KMER_LENGTH].arg);
    if (tag_length > 32) {
      cerr << "ERROR: For now, K-mer length has to be < 32" << "\n\n";
      option::printUsage(cerr, usage);
      return 1;
    }
  }

  // Maximun reads to analyze
  if (options[MAX_READS]) {
    max_reads = atoi(options[MAX_READS].arg);
  }

  //if (options[NB_THREADS]) {
  //  nb_threads = atoi(options[NB_THREADS].arg);
  //}

  // Are data stranded ?
  if (options[STRANDED]) {
    isstranded = true;
  }

  // Are data paired ?
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
      cerr << "\tPaired mode turn ON, with option " << paired << ".\n";
    }
  }

  // Data are nostranded, turn off ispaired and isstranded
  if (options[NOSTRANDED]) {
    isstranded = false;
    ispaired = false;
  }

  // Analyze all tags=kmer from the sequence
  // Generate all tags from submited sequences
  if (options[ALLTAGS]) {
    doalltags = true;
  }

  // We normalize count to million
  if (options[NORMALIZE]) {
    normalize = true;
    normalize_factors = MILLION;
  }

  // We normalize count to billion
  if (options[BILLIONOPT]) {
    normalize = true;
    normalize_factors = BILLION;
  }

  // Tag names are print in stdout
  if (options[TAG_NAMES]) {
    print_tag_names = true;
  }

  // We merge all count in one column
  if (options[MERGE_COUNTS]) {
    merge_counts = true;
  }
  // Give a name to the merged column instead of 'counts'
  if (options[MERGE_COUNTS_COLNAME]) {
    // can use only --merge-counts-colname option
    // but set merge_counts to use only one colname
    merge_counts = true;
    if (verbose>2) {
      cerr << "\tColumn name when merging:" << options[MERGE_COUNTS_COLNAME].arg << "\n";
    }
    merge_counts_colname = options[MERGE_COUNTS_COLNAME].arg;
  }

  // Write reads matching a kmer in a file
  if (options[READS_WRFILE].count()) {
    output_read = options[READS_WRFILE].arg;
  }

  // output the summary in a file
  if (options[SUMMARY].count()) {
    summary_file = options[SUMMARY].arg;
    // turn verbose to 1 at least
    verbose = verbose ? verbose++ : 1;
  }

  // Print help for unknown option
  if (options[UNKNOWN]) {
    for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
        cerr << "Unknown option: " << opt->name << "\n";
    option::printUsage(cerr, usage);
    return 1;
  }

  // Need at least one fastq file on argument line
  if(parse.nonOptionsCount() < 1) {
    cerr << "No fastq file provided ?" << "\n";
    option::printUsage(cerr, usage);
    return 0;
  }
  // We have fastq, store number of element in argument line = nb samples
  nb_samples = parse.nonOptionsCount();

  // Tag filename is require, so is always defined here
  tags_file = options[TAG_FILE].arg;

  // print info for verbose
  if (verbose>2) {
    cerr <<  "Version: " << VERSION << endl;
    cerr << "File to analyse: " << to_string(parse.nonOptionsCount()) << endl;
    for (int i = 0; i < parse.nonOptionsCount(); ++i)
      fprintf(stderr, "Non-option argument #%d is %s\n", i, parse.nonOption(i));
  }

  /**********************************
   *
   *    Variables
   *
   *********************************/

  uint32_t nb_tags = 0;                          // store number of tags read
  uint64_t tag;                                  // tag int converted from sequence

  // hash table of vector to store tag+count
  unordered_map<uint64_t,double*> tags_counts;
  // iterators for tags_counts
  unordered_map<uint64_t,double*>::iterator it_counts;
  // hash table to store tag_name
  unordered_map<uint64_t,vector<string>> tags_names;
  // vector to store nb factors = kmer per sample
  vector<uint64_t> nb_factors_by_sample;
  // vector to store nv reads per sample
  vector<uint64_t> nb_reads_by_sample;

  ofstream hfile_summary;                        // file handle to write summary

  /**********************************
   *
   *     Create hash tags count table
   *
   *********************************/

  // local vars
  bool no_name = 1;                              // do we find a tag_name
  string tag_name;                               // store tag_name read or generated
  unique_ptr< istream > filein;                  // filehandle for tag file
  uint32_t nline_tag = 0;                        // store number of line read in tag_file

  if (verbose > 1)
    cerr << "Counting k-mers" << endl;

  // read tags from stdin or file with zstr stream,
  // which manage if gz or not
  try {
    if (*tags_file == '-') {
      // use stdin
      filein = unique_ptr< istream >(new zstr::istream(cin));
    } else {
      filein = unique_ptr< istream >(new zstr::ifstream(tags_file));
    }
  }
  // check if no read error, catch all exceptions
  catch (...){
    cerr << "Error: Can't read "<< tags_file << endl;
    exit(10);
  }

  // Parse file and detect file format (fas, raw or tsv)
  for (string lines; getline(*filein, lines); ) {
    if (lines.find(">") != string::npos) {
      // we got a fasta line
      no_name = 0;
      tag_name = lines;
      tag_name.erase(0,1); // remove the fasta ">" prefix
      if (verbose>2)
        cerr << "Find fasta line: " << tag_name << endl;
    } else {
      // Check if we have a raw or tsv line: tag\tname
      size_t found = lines.find_first_of("\t,; ");
      if (found != string::npos) {
        no_name = 0;
        tag_name = lines.substr(found+1);
        lines.erase(found, lines.length());
      }
      // insert nb_tags as tag_name if none read
      if (no_name)
        tag_name = to_string(nb_tags+1);
      // Take at least tag_length for each tags
      if (lines.length() < tag_length) {
        cerr << "Error: tag lower than kmer: " << lines << endl;
        continue;                                // go to next tag
      }
      // convert tag to Int
      uint32_t max_kmer_toread = 1;              // at least we get the first kmer

      // take all kmer from lines if ALLTAGS option is set
      if (doalltags) {
        max_kmer_toread = lines.length() - tag_length + 1;
      }

      for (uint32_t i = 0; i < max_kmer_toread ; i++){
        tag = DNAtoInt(lines.substr(i, tag_length).c_str(), tag_length, isstranded);
        if (verbose>2)
          cerr << "tag: " << lines.substr(i, tag_length) << ", name:" << tag_name << ", tag nb: " << i << ", tagInt: " << tag;
        // Create vector for each tag
        tags_counts[tag] = new double[nb_samples]();
        // Add tag_name or tag_name.kmer if take all tags from sequence
        if (! doalltags)
          tags_names[tag].push_back(tag_name);
        else {
          // use a temp string to avoid multiple add of .1.2.3.4 ...
          string tempstr = tag_name;
          tags_names[tag].push_back(tempstr.append(".").append(to_string(i)));
        }
        nb_tags++;
      }
      if (verbose>2)
       cerr << ", nb_tag: " << nb_tags << endl;
    }
    nline_tag++;
  }

  // We test that tag file is not empty
  if (nline_tag == 0) {
    cerr << "I did not get or understand your tag sequences" << endl;
    exit(2);
  }

   if (verbose > 1)
     cerr << "Finished indexing tags" << endl;

  /**********************************
   *
   *       Count kmer in fastq
   *
   *********************************/

  // local vars
  ofstream hfile_read;                           // file handle to write read of matching kmer
  string gzip_pipe = "gunzip -fc ";              // string to do pipe easily to decrompress fastq.gz file

  // open file to output reads matching kmer
  if (output_read.length()) {
    hfile_read.open(output_read, ifstream::out);
    // check if not write error
    if (hfile_read.fail()) {
      cerr << "Error: Can't write read matching kmer in file "<< output_read << endl;
      return 1;
    }
  }

  // open file tp output summary information
  if (summary_file.length()) {
    hfile_summary.open(summary_file, ifstream::out);
    // check if not write error
    if (hfile_summary.fail()) {
      cerr << "Error: Can't write to summary file "<< summary_file << endl;
      return 1;
    }
    // send cerr to hfile_summary;
    cerr.rdbuf(hfile_summary.rdbuf());
  }

  // print arguments use to summary file
  if (verbose) {
    cerr << "CountTags version\t" << VERSION << "\n";
    cerr << "Kmer_size\t" << tag_length << "\n";
    cerr << "Tag file in\t" << tags_file << "\n";
    if (max_reads < UINT32_MAX)
      cerr << "Maximun reads analyzed\t" << max_reads << "\n";
    cerr << "Normalize\t" << (normalize ? "Yes" : "No") << "\n";
    if (normalize) {
      cerr << "Normalize by " << normalize_factors << " of kmer." << "\n";
    }
    cerr << "Stranded\t" << (isstranded ? "Yes" : "No") << "\n";
    cerr << "Paired\t" << (ispaired ? paired : "No") << "\n";
    cerr << "Merge count\t" << (merge_counts ? "Yes" : "No") << "\n";
    cerr << "Write matched read in file\t" << (output_read.length() ? output_read : "None") << "\n";
  }

  // Read the fastq
//#pragma omp parallel num_threads(nb_threads)
  for (uint32_t sample = 0; sample < nb_samples; ++sample) {
    if (verbose > 1)
       cerr << "Counting tags for file: " << "\t" << parse.nonOption(sample) << "\n";

    // local vars specific to each sample
    FILE * hfastq;                               // handle to fastq file
    uint32_t seq_length = 0;                     // length of the read
    uint32_t nread = 0;                          // number of read analyzed

    // use c funtionc getline which allocate memory if line=NULL & len=0,
    // line has to be freed at the end
    char * line = NULL;                          // char* to each line read in hfastq
    size_t len = 0;                              // length of buffer line read by getline
    ssize_t linelen;                             // length of line read by getline, include \0
    uint32_t nline_read = 0;                     // store number of line read in hfastq

    // Create the *char to store the tag sequence if output_read is yes and there is match
    // c style string
    char *tag_seq = new char[tag_length+1];
    tag_seq[tag_length] = '\0';

    uint64_t nb_factors = 0;                     // number of factors = kmer in a sample
    // values to DNAtoInt
    uint64_t valrev = 0;
    uint64_t valfwd = 0;
    int64_t last = 0;
    string cmdline = "";                         // command line to pass to pipe via popen

    // Test if fastq file is present, otherwise exit with error 10
    ifstream testfile(parse.nonOption(sample));
    if (!testfile.good()) {
      cerr << "Error: Can't read fastq file " << parse.nonOption(sample) << endl;
      exit(10);
    }

    // open file via pipe, so using another thread to gunzip the file
    hfastq = popen(cmdline.append(gzip_pipe).append(parse.nonOption(sample)).c_str(), "r");

    // ispaired: have to get reverse complement for reverse pair
    // use bool getrev to get the reverse complement when rf/fr/ff
    bool getrev = false;

    if (ispaired) {
      if ( string(parse.nonOption(sample)).find("_1.fastq") != string::npos) {
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
      cerr << "Paired mode ON, getrev = " << to_string(getrev) << ", for file " << parse.nonOption(sample) << endl;

    while ((linelen = getline(&line, &len, hfastq)) != -1) {
      // If this line is a sequence
      if (nline_read % 4 == 1) {
        // how many read are analyzed
        nread = ((int)((double)nline_read*0.25) + 1);
        if (nread >= max_reads) {
          break;
        }
        // Print a user-friendly output on STDERR every each XXXX reads processed
        if (verbose > 1 && nread % MILLION == 0) {
          cerr << nread << " reads parsed" << endl;
        }
        // set seq to line
        seq_length = linelen - 1; // Minus 1 because we have a new line character
        // Skip the sequence if the read length is < to the tag_length
        if (seq_length < tag_length)
          continue;

        nb_tags = seq_length - tag_length + 1;
        nb_factors += nb_tags;
        last = -3;

        for (uint32_t i = 0; i < nb_tags; i++) {
          it_counts = tags_counts.find(valns(i, line, tag_length, &last, &valfwd, &valrev, isstranded, getrev));
          // did we find a tag in the seq
          if (it_counts != tags_counts.end()) {
            it_counts->second[sample]++;
            // output read in a file if required
            if (output_read.length()) {
                intToDNA(it_counts->first, tag_length, tag_seq);
                hfile_read << tag_seq;
                if(print_tag_names) {
                  hfile_read << "\t" << join(tags_names[it_counts->first], ",");
                }
              hfile_read << "\t" << parse.nonOption(sample) << "\t" << line;
            }
          }
        }
      }
      nline_read++;
    }                                            // end of while reading hfastq file

    // store statistic
    nb_factors_by_sample.push_back(nb_factors);
    nb_reads_by_sample.push_back(nread);

    if (normalize && nb_factors > 0) {
      if (verbose > 1 )
        cerr << "Normalize counts" << endl;
      for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
        if (it_counts->second[sample] > 0) {
          it_counts->second[sample] = it_counts->second[sample] * normalize_factors / nb_factors;
        }
      }
    }

    // Close file and clear line buffer
    pclose(hfastq);
    if (line)
      free(line);
  }                                              // end of for each sample

  /**********************************
   *
   *     PRINT THE RESULTS
   *
   *********************************/

  // First print headers
  cout << "tag";
  if (print_tag_names)
    cout << "\ttag_names";
  if(!merge_counts) {
    for (uint32_t sample = 0; sample < nb_samples; ++sample) {
      // keep basename for fastq filename
      size_t found = string(parse.nonOption(sample)).find_last_of("/");
      if (found == string::npos)
        found = -1;      // set to the beginning if no path present
      cout << "\t" << string(parse.nonOption(sample)).substr(found+1);
    }
  } else {
    cout << "\t" << merge_counts_colname;
  }
  cout << "\n";
  // print tag + value
  char *tag_seq = new char[tag_length+1];
  tag_seq[tag_length] = '\0';
  for (it_counts=tags_counts.begin(); it_counts!=tags_counts.end(); ++it_counts) {
    intToDNA(it_counts->first,tag_length,tag_seq);
    cout << tag_seq;
    if(print_tag_names) {
      cout << "\t" << join(tags_names[it_counts->first],",");
    }
    if(!merge_counts){
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        cout << "\t" << it_counts->second[sample];
      }
    } else {
      double count_sum = 0;
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        count_sum += it_counts->second[sample];
      }
      cout << "\t" << count_sum;
    }
    cout << endl;
  }


  // print statistic if verbose or summary
  if (verbose) {
    cerr << "# Total statistic per file\n";
    cerr << "File\t";
    if(!merge_counts) {
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        cerr << "\t" << parse.nonOption(sample);
      }
    } else {
      cerr << "\t" << merge_counts_colname;
    }
    cerr << "\n";
    cerr << "total_factors";
    if(!merge_counts) {
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        cerr << "\t" << nb_factors_by_sample[sample];
      }
    } else {
      uint64_t nb_factors_sum = 0;
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        nb_factors_sum += nb_factors_by_sample[sample];
      }
      cerr << "\t" << nb_factors_sum;
    }
    cerr << endl;
    cerr << "total_reads";
    if(!merge_counts) {
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        cerr << "\t" << nb_reads_by_sample[sample];
      }
    } else {
      uint64_t sum = 0;
      for (uint32_t sample = 0; sample < nb_samples; ++sample) {
        sum += nb_reads_by_sample[sample];
      }
      cerr << "\t" << sum;
    }
    cerr << endl;
  }
}

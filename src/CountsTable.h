/******************************************************************************\
*                                                                              *
*  Copyright © 2015-2018 -- IRMB/INSERM                                        *
*                           (Institute for Regenerative Medicine & Biotherapy  *
*                           Institut National de la Santé et de la Recherche   *
*                           Médicale)                                          *
*                                                                              *
*  Copyright (C) 2015-2018  Jérôme Audoux                                      *
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

#ifndef COUNTS_TABLE_H
#define COUNTS_TABLE_H

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "utils.h"
#include "dna.h"

class CountsTable {

private:
  std::unordered_map<uint64_t,uint32_t*> tags_counts;
  std::vector<std::string> sample_names;
  uint nb_samples;
  uint kmer_length;
  bool stranded;

public:
  CountsTable(uint nb_samples, uint kmer_length, bool stranded = false);

  ~CountsTable();

  void setSampleName(uint sample_id, const char *sample_name);

  uint32_t getCount(const char *kmer, uint sample_id);
  uint32_t getCount(uint64_t kmer, uint sample_id);

  /*
  * Set the counts value for one kmer of one sample
  */
  void setCount(const char *kmer, uint sample_id, uint value);
  void setCount(uint64_t kmer, uint sample_id, uint value);

  /*
  * Increment the count for one kmer of one sample
  */
  void incrementCount(const char *kmer, uint sample_id, uint value);
  void incrementCount(uint64_t kmer, uint sample_id, uint value);

  /*
  * Filters counts based on their recurrency accross samples
  */
  void recurrencyFilter(uint recurrency_threshold = 2);

  /*
  * Print counts on STDOUT
  */
  void printCounts(char sep = '\t');

};

#endif // COUNTS_TABLE_H

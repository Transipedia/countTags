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

#include "CountsTable.h"

CountsTable::CountsTable(uint nb_samples, uint kmer_length, bool stranded)
  : nb_samples(nb_samples), kmer_length(kmer_length), stranded(stranded)
{
  for (size_t i = 0; i < nb_samples; i++) {
    sample_names.push_back("sample" + std::to_string(i));
  }
}

CountsTable::~CountsTable() {
  std::unordered_map<uint64_t,uint32_t*>::iterator it_counts;
  for (it_counts = tags_counts.begin(); it_counts != tags_counts.end(); ++it_counts) {
    delete[] it_counts->second;
  }
}

void CountsTable::setSampleName(uint sample_id, const char *sample_name) {
  this->sample_names[sample_id] = sample_name;
}

uint32_t CountsTable::getCount(const char *kmer, uint sample_id) {
  uint64_t kmer_int = DNAtoInt(kmer, this->kmer_length, this->stranded);
  return this->getCount(kmer_int,sample_id);
}

uint32_t CountsTable::getCount(uint64_t kmer, uint sample_id) {
  std::unordered_map<uint64_t,uint32_t*>::iterator it_counts = tags_counts.find(kmer);
  // If the tag is already in the hash, we set the count for this->sample
  if (it_counts != tags_counts.end()) {
    return it_counts->second[sample_id];
  } else {
    return 0;
  }
}

void CountsTable::setCount(const char *kmer, uint sample_id, uint value) {
  uint64_t kmer_int = DNAtoInt(kmer, this->kmer_length, this->stranded);
  this->setCount(kmer_int,sample_id,value);
}

void CountsTable::setCount(uint64_t kmer, uint sample_id, uint value) {
  std::unordered_map<uint64_t,uint32_t*>::iterator it_counts = tags_counts.find(kmer);
  // If the tag is already in the hash, we set the count for this->sample
  if (it_counts != tags_counts.end()) {
    it_counts->second[sample_id] = value;
  // Otherwise we create a new entry
  } else {
    uint32_t * tag_counts = new uint32_t[this->nb_samples]();
    tag_counts[sample_id] = value;
    tags_counts[kmer] = tag_counts;
  }
}

void CountsTable::incrementCount(const char *kmer, uint sample_id, uint value) {
  uint64_t kmer_int = DNAtoInt(kmer, this->kmer_length, this->stranded);
  this->incrementCount(kmer_int,sample_id,value);
}

void CountsTable::incrementCount(uint64_t kmer, uint sample_id, uint value) {
  uint32_t count = this->getCount(kmer,sample_id);
  this->setCount(kmer, sample_id, value + count);
}

void CountsTable::recurrencyFilter(uint recurrency_threshold) {
  std::unordered_map<uint64_t,uint32_t*>::iterator it_counts = this->tags_counts.begin();
  while (it_counts != this->tags_counts.end()) {
    uint recurrency = 0;
    for (uint s = 0; s < this->nb_samples; ++s) {
      if (it_counts->second[s] > 0) {
        recurrency++;
      }
    }
    if (recurrency >= recurrency_threshold) {
      it_counts++;
    } else {
      it_counts = this->tags_counts.erase(it_counts);
    }
  }
}

void CountsTable::printCounts(char sep) {
  char *tag_seq = new char[this->kmer_length];
  std::unordered_map<uint64_t,uint32_t*>::iterator it_counts;

  /* Print Headers */
  std::cout << "tag";
  for (uint s = 0; s < this->nb_samples; ++s) {
    std::cout << sep << this->sample_names[s];
  }
  std::cout << std::endl;

  /* Print Counts */
  for (it_counts = tags_counts.begin(); it_counts != tags_counts.end(); ++it_counts) {
    intToDNA(it_counts->first,this->kmer_length,tag_seq);
    std::cout << tag_seq;
    for (uint s = 0; s < this->nb_samples; ++s) {
      std::cout << "\t" << it_counts->second[s];
    }
    std::cout << std::endl;
  }
}

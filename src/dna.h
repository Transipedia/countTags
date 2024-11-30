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

#ifndef DNA_H
#define DNA_H

#include <cstdint>

#include "utils.h"

inline uint convNuc(char nuc){
  switch (nuc){
    case 'a' : case 'A' : return 0 ;
    case 'c' : case 'C' : return 1 ;
    case 'g' : case 'G' : return 2 ;
    case 't' : case 'T' : return 3 ;
    default : return 0;
  }
  return 0;
}

inline uint compNuc(uint nuc){
  switch (nuc){
    case 0 : return 3 ;
    case 1 : return 2 ;
    case 2 : return 1 ;
    case 3 : return 0 ;
    default  : return 0;
  }
  return 0;
}

inline char intToNuc(uint c) {
  switch(c) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  default : return 'A';
  }
}

//int comparUint(const void *a1, const void* a2){
//  return (int) (* (uint64_t *)a1 - * (uint64_t *) a2);
//}

inline void intToDNA(uint64_t code, uint dna_length, char *dna) {
  uint64_t mask = 3;
  for (uint i=0; i < dna_length; i++) {
    dna[dna_length-i-1] = intToNuc(code & mask);
    code >>=2;
  }
}

inline uint64_t intRevcomp(uint64_t factor, uint32_t length) {
  uint64_t mask;
  if (length == 32)
    mask = ~0;
  else
    mask =  ((uint64_t) 1 << (2*length)) - 1;

  factor ^= mask;

  uint64_t mask_lsb;
  // Corresponds to the rightmost nucleotide
  mask_lsb = 3;
  uint64_t shift = 0;
  uint64_t result=0;
  for (uint32_t j(0);j<length;j++){
    result <<= 2;
    // get the leftmost nucleotide and put it at the end
    result |= (factor & mask_lsb) >> shift;
    mask_lsb <<= 2;
    shift += 2;
  }

  return result;
}

inline uint64_t DNAtoInt(const char *dna, uint32_t dna_length, bool stranded = false){
  uint64_t dna_int = 0;
  for (uint32_t i = 0; i< dna_length ; i++){
    dna_int <<= 2;
    dna_int |= convNuc(dna[i]);
  }
  // If the conversion is not "strand-specific" we calculate the reverse DNA it
  // and return the one that has the smallest value
  if (!stranded) {
    uint64_t rev_comp = intRevcomp(dna_int,dna_length);
    if (rev_comp < dna_int) {
      return rev_comp;
    }
  }
  return dna_int;
}

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

#endif // DNA_H_

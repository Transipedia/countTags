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
#include <sstream>
#include <string>
#include <vector>

#ifndef UTILS_H
#define UTILS_H

#define uint unsigned int

string join(const vector<string>& elements, const char* const separator)
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

/* a read line from fastq and check if not at EOF */
bool read_aline (string &astr, FILE * hfile)
{
  // use c funtion getline, for popen function, which allocate memory if line=NULL & len=0,
  // line has to be freed at the end
  char * line = NULL;                          // char* to each line read in hfastq
  size_t len = 0;                             // length of buffer line read by getline
  ssize_t linelen;                             // length of line read by getline, include \0
  linelen = getline(&line, &len, hfile);
  // remove newline char
  //cerr << "\tgetline : " << line << ", len : " << linelen << endl;
  if (linelen == -1) {
    astr = "";
    return 0;
  }
  astr = line;
  astr.erase(astr.size() - 1);
  // free memory to avoid memory leak
  if (line)
    free(line);
  return 1;
}

#endif

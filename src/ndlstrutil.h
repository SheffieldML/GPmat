#ifndef NDLSTRUTIL_H
#define NDLSTRUTIL_H
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include "ndlexceptions.h"
namespace ndlstrutil
{
  // allocate doubles in 100k chunks.
  const int ALLOCATECHUNK=(int)(100*(double)1024/(double)sizeof(double));
  // split a string into tokes given a delimiter.
  void tokenise(std::vector<std::string>& tokens,
		const std::string& str,
		const std::string& delimiters = " ");
  // a version of getline which ignores lines starting with #.
  bool getline(std::istream& in, std::string& line);
  // take a string and sent it to an output stream placing carriage returns correctly.
  void wrapOutputText(std::ostream& out, const std::string description, const int width, const int padding);
  // convert an integer to a string.
  std::string itoa(long value, int base=10);
  std::string dirSep();
}

#endif

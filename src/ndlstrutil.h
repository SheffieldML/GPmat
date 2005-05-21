#ifndef NDLSTRUTIL_H
#define NDLSTRUTIL_H
#include <vector>
#include <string>
#include <iostream>

namespace ndlstrutil
{
  const int ALLOCATECHUNK=(int)(100*(double)1024/(double)sizeof(double));
  void tokenise(std::vector<std::string>& tokens,
		const std::string& str,
		const std::string& delimiters = " ");
  // my version of getline which ignores lines starting with #.
  void getline(std::istream& in, std::string& line);
  void wrapOutputText(std::ostream& out, const std::string description, const int width, const int padding);

  // allocate doubles in 100k chunks.
}

#endif

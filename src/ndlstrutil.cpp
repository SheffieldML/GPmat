#include "ndlstrutil.h"

std::string ndlstrutil::dirSep()
{
  return "/";
}
void ndlstrutil::tokenise(std::vector<std::string>& tokens,
			  const std::string& str,
			  const std::string& delimiters)
{
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos!=pos 
	 || std::string::npos!=lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

bool ndlstrutil::getline(std::istream& in, std::string& line)
{
  bool val = false;
  do
  {
    if(std::getline(in, line))
      val = true;
  }
  while(line[0]=='#');
  return val;
}
void ndlstrutil::wrapOutputText(std::ostream& out, const std::string description, const int width, const int padding)
{
  std::string padStr="";
  for(int i=0; i<padding; i++)
    padStr += " ";
  std::vector<std::string> tokens;
  ndlstrutil::tokenise(tokens, description);
  size_t remain = width-padding;
  out << padStr;
  for(size_t i=0; i<tokens.size(); i++)
  {
    remain -= (tokens[i].size()+1);
    if (remain<1)
    {
      out << std::endl;
      out << padStr;
      remain=width-padding-tokens[i].size()-1;
    }
    out << tokens[i] << " ";
  }
}
std::string ndlstrutil::itoa(long value, int base)
{
  enum { kMaxDigits = std::numeric_limits<long>::digits };
  std::string buf;
  buf.reserve( kMaxDigits ); // Pre-allocate enough space.
  if (base < 2 || base > 16) throw ndlexceptions::Error("Unrecognised base in ndlstrutil::itoa()");
  long quotient = value;
  do 
  {
    buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
    quotient /= base;
  } 
  while ( quotient );
  // Append the negative sign for base 10
  if ( value < 0 && base == 10) buf += '-';
  std::reverse( buf.begin(), buf.end() );
  return buf;
}

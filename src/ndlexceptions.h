// This is very basic exception handling. It should probably be better integrated with the C++ exceptions.
#ifndef NDLEXCEPTIONS_H
#define NDLEXCEPTIONS_H
#include <vector>
#include <string>
#include <iostream>

namespace ndlexceptions 
{
class Error {
 public:
  Error(){}
  virtual ~Error(){}
  Error(std::string message)
    {
      msg = message;
    }
  virtual std::string getMessage() const
    {
      return msg;
    }
  void setMessage(const std::string message)
    {
      msg = message;
    }
 private:
  std::string msg;
};
class NotImplementedError : public ndlexceptions::Error {
 public:
  NotImplementedError(){}
  NotImplementedError(std::string message)
    {
      setMessage(message);
    }
};
class FileError : public ndlexceptions::Error {
 public:
  FileError(){}
  FileError(std::string file) : fileName(file) {}
  std::string getMessage() const
    {
      return "File error in: " + getFileName();
    }
  virtual void setFileName(const std::string file)
    {
      fileName = file;
    }
  std::string getFileName() const
    {
      return fileName;
    }
 private:
  std::string fileName;
};
class FileReadError : public ndlexceptions::FileError {
  
 public:
  FileReadError(){}
  FileReadError(std::string file)
    {
      setFileName(file);
    }
  std::string getMessage() const
    {
      return "Unable to read file: " + getFileName();
    }
};
class FileWriteError : public ndlexceptions::FileError {
  
 public:
  FileWriteError(){}
  FileWriteError(std::string file)
    {
      setFileName(file);
    }

  std::string getMessage() const
    {
      return "Unable to write to file: " + getFileName();
    }
};
class FileFormatError : public ndlexceptions::FileReadError {
 public:
  FileFormatError(){}
  FileFormatError(std::string file) 
    {
      setFileName(file);
    }
  std::string getMessage() const
    {
      return "Incorrect file format in: " + getFileName();
    }
};
class FileVersionError : public ndlexceptions::FileFormatError {
 public:
  FileVersionError(){}
  FileVersionError(std::string file) 
    {
      setFileName(file);
    }
  std::string getMessage() const
    {
      return "Incorrect file version in: " + getFileName();
    }
};

class MatrixError : public Error {

 public:
  MatrixError(){}
  MatrixError(std::string message) 
    {
      setMessage(message);
    }

};
class MatrixNonPosDef : public MatrixError {
 public:
  MatrixNonPosDef()
    {
      setMessage("Matrix non positive definite.");
    }
};
class MatrixConditionError : public MatrixError {
 public:
  MatrixConditionError()
    {
      setMessage("Matrix has low condition number.");
    }
};
class MatrixSingular : public MatrixError {
 public:
  MatrixSingular()
    {
      setMessage("Matrix is singular.");
    }
};
 
class CommandLineError : public Error {
 public:
  CommandLineError(){}
  CommandLineError(std::string message) 
    {
      setMessage(message);
    }
};

}
#endif

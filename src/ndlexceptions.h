// This is very basic exception handling. It should probably be better integrated with the C++ exceptions.
#ifndef NDLEXCEPTIONS_H
#define NDLEXCEPTIONS_H
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

namespace ndlexceptions 
{
  class Error : public std::exception
  {
  public:
    Error(){}
    virtual ~Error() throw (){}
    Error(std::string message) 
    {
      msg = message;
    }
    virtual std::string getMessage() const
    {
      return msg;
    }
    virtual void setMessage(const std::string message)
    {
      msg = message;
    }
    virtual const char* what() const throw()
    {
      return "Unknown ndlexception";
    }

  private:
    std::string msg;
  };
  class RuntimeError : public ndlexceptions::Error
  {
  public:
  RuntimeError() : Error() {}
  RuntimeError(std::string message) : Error(message) {}
    ~RuntimeError() throw() {}
    virtual const char* what() const throw()
    {
      return "Runtime error";
    }
  };
  class TypeError : public ndlexceptions::Error
  {
  public:
  TypeError() : Error() {}
  TypeError(std::string message) : Error(message) {}
    ~TypeError() throw() {}
    virtual const char* what() const throw()
    {
      return "Type error";
    }
  };
  class MatlabInterfaceError : public ndlexceptions::Error 
  {
  public:
  MatlabInterfaceError() : Error() {}
  MatlabInterfaceError(std::string message) : Error(message) {}
    ~MatlabInterfaceError() throw() {}
    virtual const char* what() const throw()
    {
      return "MATLAB interface error";
    }

  };
  class MatlabInterfaceReadError : public ndlexceptions::MatlabInterfaceError 
  {
  public:
  MatlabInterfaceReadError() : MatlabInterfaceError() {}
  MatlabInterfaceReadError(std::string fieldName) : MatlabInterfaceError(fieldName)     
    {
      setFieldName(fieldName);
    }
    ~MatlabInterfaceReadError() throw() {}
    void setFieldName(const std::string name)
    {
      fieldName = name;
    }
    std::string getFieldName() const
    {
      return fieldName;
    }
    std::string getMessage() const
    {
      return "Error in MATLAB interface reading " + getFieldName() + ".";
    }
    virtual const char* what() const throw()
    {
      return "MATLAB interface read error";
    }
  private:
    std::string fieldName;
  };
  class NotImplementedError : public ndlexceptions::Error 
  {
  public:
  NotImplementedError() : Error() {}
  NotImplementedError(std::string message) : Error(message) {}
    virtual const char* what() const throw()
    {
      return "Functionality not implemented.";
    }
  };

  class FileError : public ndlexceptions::Error 
  {
  public:
  FileError() : Error() {}
  FileError(std::string file) : fileName(file), Error("File error " + file) {}
    virtual ~FileError() throw() {}
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
    virtual const char* what() const throw()
    {
      return "File error";
    }
  private:
    std::string fileName;
  };

  class FileReadError : public ndlexceptions::FileError 
  {  
  public:
  FileReadError() : FileError() {}
  FileReadError(std::string file) : FileError("(reading) " + file)
      {
	setFileName(file);
      }
    virtual ~FileReadError() throw() {}
    virtual  std::string getMessage() const
    {
      return "Unable to read file: " + getFileName();
    }
    virtual const char* what() const throw()
    {
      return "File read error";
    }
   };


  class StreamFormatError : public ndlexceptions::Error
  {
  public:
  StreamFormatError() : Error() {}
  StreamFormatError(std::string field) : Error("Error when expecting field: " + field)
    {
	setFieldName(field);
    }
  StreamFormatError(std::string field, std::string note) : Error(field + ": " + note)
    {
      setFieldName(field);
    }
    virtual ~StreamFormatError() throw() {}
    std::string getFieldName() const
    {
      return fieldName;
    }
    void setFieldName(const std::string name) 
    {
      fieldName = name;
    }
    virtual const char* what() const throw()
    {
      return "Stream format error";
    }
 
  private:
    std::string fieldName;
  };

  class StreamVersionError : public ndlexceptions::StreamFormatError 
  {
  public:
  StreamVersionError() : StreamFormatError("version") {}
    virtual std::string getMessage() const
      {
	return "Incorrect version.";
      }
    virtual const char* what() const throw()
    {
      return "Stream version error";
    }
 
  };

  class FileFormatError : public ndlexceptions::FileError
  {
  public:
  FileFormatError() : FileError() {}
  FileFormatError(std::string fileName, ndlexceptions::StreamFormatError& err) : FileError(err.getMessage() + " in " + fileName) {}
  FileFormatError(std::string fileName, ndlexceptions::MatlabInterfaceReadError& err) : FileError(err.getMessage() + " in " + fileName) {}
  FileFormatError(std::string fileName) : FileError("File format error in  " + fileName) {}
    virtual ~FileFormatError() throw() {}
    virtual const char* what() const throw()
    {
      return "File format error";
    }
 
  };
  class FileWriteError : public ndlexceptions::FileError 
  {  
  public:
  FileWriteError() : FileError() {}
  FileWriteError(std::string file) : FileError(file) {}
    
    std::string getMessage() const
    {
	return "Unable to write to file: " + getFileName();
    }
    virtual const char* what() const throw()
    {
      return "File write error";
    }
 
  };

  class MatrixError : public Error 
  {

  public:
  MatrixError() : Error(){}
  MatrixError(std::string message)  : Error(message) {}
    virtual const char* what() const throw()
    {
      return "Matrix error";
    }
 

  };
  class MatrixNonPosDef : public MatrixError 
  {
  public:
  MatrixNonPosDef() : MatrixError("Matrix non positive definite error") {}
    virtual const char* what() const throw()
    {
      return "Non positive definite matrix.";
    }
 
  };
  class MatrixConditionError : public MatrixError 
  {
  public:
  MatrixConditionError() : MatrixError("Matrix has low condition number.") {}
    virtual const char* what() const throw()
    {
      return "Matrix condition error";
    }
 
  };
  class MatrixSingular : public MatrixError 
  {
  public:
  MatrixSingular() : MatrixError("Matrix is singular.") {}
    virtual const char* what() const throw()
    {
      return "Matrix is singular error";
    }
 
  };
 
  class CommandLineError : public Error 
  {
  public:
  CommandLineError() : Error() {} 
  CommandLineError(std::string message) : Error(message) {}
    virtual const char* what() const throw()
    {
      return "Command line format error";
    }
    
  };

}
#endif

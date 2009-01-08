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

  class MatlabInterfaceError : public ndlexceptions::Error 
  {
  public:
    MatlabInterfaceError(){}
    MatlabInterfaceError(std::string message)
    {
      setMessage(message);
    }
    ~MatlabInterfaceError() throw() {}
    virtual const char* what() const throw()
    {
      return "MATLAB interface error";
    }

  };
  class MatlabInterfaceReadError : public ndlexceptions::MatlabInterfaceError 
  {
  public:
    MatlabInterfaceReadError(){}
    MatlabInterfaceReadError(std::string fieldName)
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
    NotImplementedError()
    {
      setMessage("");
    }
    NotImplementedError(std::string message)
    {
      setMessage(message);
    }
    virtual const char* what() const throw()
    {
      return "Functionality not implemented.";
    }
  };

  class FileError : public ndlexceptions::Error 
  {
  public:
    FileError(){}
    FileError(std::string file) : fileName(file) {}
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
    FileReadError(){}
    FileReadError(std::string file)
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
    StreamFormatError(){}
    StreamFormatError(std::string field) 
    {
      setFieldName(field);
      setMessage("Error when expecting field: " + field);
    }
    StreamFormatError(std::string field, std::string note) 
    {
      setFieldName(field);
      setMessage(field + ": " + note);
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
    StreamVersionError() 
    {
      setFieldName("version");
    }
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
    FileFormatError(){}
    FileFormatError(std::string fileName, ndlexceptions::StreamFormatError& err)
    {
      setMessage(err.getMessage() + " in " + fileName);
    }
    FileFormatError(std::string fileName, ndlexceptions::MatlabInterfaceReadError& err)
    {
      setMessage(err.getMessage() + " in " + fileName);
    }
    virtual ~FileFormatError() throw() {}
    virtual const char* what() const throw()
    {
      return "File format error";
    }
 
  };
  class FileWriteError : public ndlexceptions::FileError 
  {  
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
    virtual const char* what() const throw()
    {
      return "File write error";
    }
 
  };

  class MatrixError : public Error 
  {

  public:
    MatrixError(){}
    MatrixError(std::string message) 
    {
      setMessage(message);
    }
    virtual const char* what() const throw()
    {
      return "Matrix error";
    }
 

  };
  class MatrixNonPosDef : public MatrixError 
  {
  public:
    MatrixNonPosDef()
    {
      setMessage("Matrix non positive definite error");
    }
    virtual const char* what() const throw()
    {
      return "Non positive definite matrix.";
    }
 
  };
  class MatrixConditionError : public MatrixError 
  {
  public:
    MatrixConditionError()
    {
      setMessage("Matrix has low condition number.");
    }
    virtual const char* what() const throw()
    {
      return "Matrix condition error";
    }
 
  };
  class MatrixSingular : public MatrixError 
  {
  public:
    MatrixSingular()
    {
      setMessage("Matrix is singular.");
    }
    virtual const char* what() const throw()
    {
      return "Matrix is singular error";
    }
 
  };
 
  class CommandLineError : public Error 
  {
  public:
    CommandLineError(){}
    CommandLineError(std::string message) 
    {
      setMessage(message);
    }
    virtual const char* what() const throw()
    {
      return "Command line format error";
    }
 
  };

}
#endif

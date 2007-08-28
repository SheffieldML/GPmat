// This is very basic exception handling. It should probably be better integrated with the C++ exceptions.
#ifndef NDLEXCEPTIONS_H
#define NDLEXCEPTIONS_H
#include <vector>
#include <string>
#include <iostream>

namespace ndlexceptions 
{
  class Error 
  {
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
    virtual void setMessage(const std::string message)
    {
      msg = message;
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
  };
  class MatlabInterfaceReadError : public ndlexceptions::MatlabInterfaceError 
  {
  public:
    MatlabInterfaceReadError(){}
    MatlabInterfaceReadError(std::string fieldName)
    {
      setFieldName(fieldName);
    }
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
      return "Error in MATLAB interface reading " + getFieldName();
    }
  private:
    std::string fieldName;
  };
  class NotImplementedError : public ndlexceptions::Error 
  {
  public:
    NotImplementedError(){}
    NotImplementedError(std::string message)
    {
      setMessage(message);
    }
  };

  class FileError : public ndlexceptions::Error 
  {
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

  class FileReadError : public ndlexceptions::FileError 
  {  
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


  class StreamFormatError : ndlexceptions::Error
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
    std::string getFieldName() const
    {
      return fieldName;
    }
    void setFieldName(const std::string name) 
    {
      fieldName = name;
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
    std::string getMessage() const
    {
      return "Incorrect version.";
    }
  };

  class FileFormatError : public ndlexceptions::FileReadError
  {
  public:
    FileFormatError(){}
    FileFormatError(std::string fileName, ndlexceptions::StreamFormatError err)
    {
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
  };

  class MatrixError : public Error 
  {

  public:
    MatrixError(){}
    MatrixError(std::string message) 
    {
      setMessage(message);
    }

  };
  class MatrixNonPosDef : public MatrixError 
  {
  public:
    MatrixNonPosDef()
    {
      setMessage("Matrix non positive definite.");
    }
  };
  class MatrixConditionError : public MatrixError 
  {
  public:
    MatrixConditionError()
    {
      setMessage("Matrix has low condition number.");
    }
  };
  class MatrixSingular : public MatrixError 
  {
  public:
    MatrixSingular()
    {
      setMessage("Matrix is singular.");
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
  };

}
#endif

/**
* EasyReader
* Provides a simple interface for text stream input 
* @author Guy J Brown
* Contains some code from Jim McGregor's SimpleReader class (used with permission)
* First version 16/2/2000
* Modification history:
* Modified 19/1/2001 to fix bug with readChar, which now ignores carriage return character
*/

import java.io.*;

public class EasyReader extends BufferedReader {

    // instance variables
    
    private boolean endOfFile = false;
    private boolean startOfFile = true;
    private char bufferChar = ' ';
    
    // Constructors
    
    /**
     * Create a new EasyReader that reads from standard input
     */
    public EasyReader() {
	super(new InputStreamReader(System.in));
    }
    
    /**
     * Create a new EasyReader that reads from a text file
     * @param s the name of the text file
     */
    public EasyReader(String s) {
	super(new InputStreamReader(getFileInputStream(s)));
    }
    
    // Private methods
    
    private static FileInputStream getFileInputStream(String s) {
	FileInputStream f = null;
	try {
	    f = new FileInputStream(s);
	}
	catch (IOException e) {
	    error("file not found");
	}
	return f;
    }
    
    private static void error(String s) {
	System.err.println("Error in EasyReader - "+s);
	System.err.flush();
	System.exit(0);
    }
    
    private void prompt(String p) {
	System.err.print(p);
	System.err.flush();
    }
    
    private void getNextChar() {
	if (!endOfFile) {
	    int c=0;
	    try {
		c=read();
	    }
	    catch (IOException e) {
		error("could not read another character");
	    }
	    if (c==-1) 
		endOfFile=true;
	    else
		bufferChar=(char)c;
	}
    }
    
    // methods
    
    /**
     * Check whether the end of the input stream has been reached
     * @return true if the end of the stream is reached
     */
    public boolean eof() {
	return endOfFile;
    }
    
    /**
     * Read a character from the input stream
     * @return a character value
     */
    public char readChar() {
	if (startOfFile) {
	    startOfFile=false;
	    getNextChar();
	}
	char ch;
	do {
	    ch = bufferChar;
	    getNextChar();
	} while (ch=='\n');
	return ch;
    }
    
    /**
     * Read a character from the input stream, with a prompt
     * @param s the prompt
     * @return a character value
     */
    public char readChar(String s) {
	prompt(s);
	return readChar();
    }
    
    /**
     * Read a string from the input stream
     * @return a string value
     */
    public String readString() {
	while (!endOfFile && (bufferChar==' ' || bufferChar=='\r' || bufferChar=='\n'))
	    getNextChar();
	String s = "";
	while (!endOfFile && (bufferChar!='\r' && bufferChar!='\n')) {
	    s+=bufferChar;
	    getNextChar();
	}
	return s;
    }
    
    /**
     * Read a string from the input stream, with a prompt
     * @return a string value
     */
    public String readString(String s) {
	prompt(s);
	return readString();
    }
    
    /**
     * Read a double from the input stream
     * @return a double value
     */
    public double readDouble() {
	double x = 0;
	try {
	    x=(new Double(readString())).doubleValue();
	}
	catch (Exception e) {
	    error("invalid floating point number");
	}
	return x;
    }
    
    /**
     * Read a double from the input stream, with a prompt
     * @return a double value
     */
    public double readDouble(String s) {
	prompt(s);
	return readDouble(); 
    }
    
    /**
     * Read a float from the input stream
     * @return a float value
     */
    public float readFloat() {
	float x = 0;
	try {
	    x=(new Float(readString())).floatValue();
	}
	catch (Exception e) {
	    error("invalid floating point number");
	}
	return x;
    }
    
    /**
     * Read a float from the input stream, with a prompt
     * @return a float value
     */
    public float readFloat(String s) {
	prompt(s);
	return readFloat();
    }
    
    /**
     * Read an integer from the input stream
     * @return an integer value
     */
    public int readInt() {
	int x = 0;
	try {
	    x=(new Integer(readString())).intValue();
	}
	catch (Exception e) {
	    error("invalid integer number");
	}
	return x;
    }
    
    /**
     * Read an integer from the input stream, with a prompt
     * @return an integer value
     */
    public int readInt(String s) {
	prompt(s);
	return readInt();
    }
    
    /**
     * Read a boolean from the input stream
     * @return a boolean value
     */
    public boolean readBoolean() {
	char c = readString().charAt(0);
	return (c=='t' || c=='T');
    }
    
    /**
     * Read a boolean from the input stream, with a prompt
     * @return a boolean value
     */
    public boolean readBoolean(String s) {
	prompt(s);
	return readBoolean();
    }		
}

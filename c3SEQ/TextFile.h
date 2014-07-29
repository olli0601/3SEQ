/* 
 * File:   TextFile.h
 * Author: Ha Minh Lam
 *
 * Created on 15 June 2011, 17:11
 */

#ifndef FILE_H
#define	FILE_H

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "String.h"

using namespace std;

class TextFile {
public:
    static const string DEFAULT_COLUMN_DELIM;
    
    explicit TextFile(string newFilePath);

    TextFile(const TextFile& orig);

    TextFile& operator=(const TextFile& rhs);

    virtual ~TextFile();

    /**
     * Get the path to the file (including the file name).
     * @return  The path to the file (including the file name).
     */
    virtual const string& getPath(void) const {
        return filePath;
    };

    /**
     * Get the in/out stream of the file.
     * @return  The in/out stream of the file.
     */
    virtual fstream* getStreamPtr(void) const {
        assert(fStreamPtr->is_open());
        return fStreamPtr;
    };

    /**
     * Delete the file if it exists.
     * @return  <b>true</b> if the file does not exist or it is successfully
     *          deleted.
     */
    virtual const bool removeFile(void) {
        if (exists()) {
            return (remove(filePath.c_str()) == 0);
        }
        return true;
    }
    
    /**
     * Write a line into file.
     * @param   line
     * @note    If the file has not opened, this method will automatically
     *          open it.<br>
     *          An "end-line" character will be automatically added at the
     *          end of the given line.
     */
    virtual void writeLine(string line) {
        if (!fStreamPtr->is_open()) {
            openToWrite();
        }
        assert(!isReadOnly);
        (*fStreamPtr) << line << endl;
    };

    /**
     * Write an empty line into file.
     */
    virtual void writeLine(void) {
        writeLine("");
    }

    /**
     * Write a string into file.
     * @param str
     * @note    If the file has not opened, this method will automatically
     *          open it.
     */
    virtual void write(string str) {
        if (!fStreamPtr->is_open()) {
            openToWrite();
        }
        assert(!isReadOnly);
        (*fStreamPtr) << str;
    };

    /**
     * Get a line from this file.
     * @return  A line of text.
     * @note    This method is safer than ifstream.getline() method since it can
     *          handle multiple end-line characters (LF, CR, CRLF).
     */
    virtual string getLine(void);
    
    /**
     * Indicate whether the file has been opened yet.
     * @return <b>true</b> if the file has been opened.
     */
    virtual const bool isOpen(void) const {
        return (fStreamPtr != NULL && fStreamPtr->is_open());
    }

    /**
     * Open the file stream to read from the beginning of the file
     */
    virtual void openToRead(void);

    /**
     * Open the file stream to write from the beginning of the file
     */
    virtual void openToWrite(void);

    /**
     * Close the file stream
     */
    virtual void close(void);

    /**
     * Check to see whether the given file exists.
     * @param fileName The file name (including path).
     * @return <b>true</b> if the file exists; <b>false</b> otherwise.
     */
    virtual const bool exists(void);

    /**
     * Indicate that the file stream is still readable (end of file and 
     * bad bits have not occurred).
     * @return <b>true</b> if the file stream is still readable.
     */
    virtual const bool isStreamReadable(void) const;

    /**
     * Retrieve a vector that contains all lines of the file.
     * @return  A vector that contains all lines of the file.
     * @note    All lines are trimmed before storing into the vector.
     */
    virtual vector<string> readAllLines(void);


private:
    static const char INVALID_PATH_CHARS[];
    static const unsigned long INVALID_CHARS_COUNT;

    static const char PATH_SEPARATORS[];
    static const unsigned long PATH_SEPARATORS_COUNT;

    /** The character that will be used to replace invalid path charaters */
    static const char INVALID_PATH_REPLACEMENT;
    
    static bool const isInvalidPathChar(char testedChar) {
        for (unsigned long i = 0; i < INVALID_CHARS_COUNT; i++) {
            if (testedChar == INVALID_PATH_CHARS[i]) {
                return true;
            }
        }
        return false;
    };
    
    static bool const isPathSeparator(char testedChar) {
        for (unsigned long i = 0; i < PATH_SEPARATORS_COUNT; i++) {
            if (testedChar == PATH_SEPARATORS[i]) {
                return true;
            }
        }
        return false;
    };

    /** Path to the file (including file name) */
    string filePath;

    /** The pointer to the file stream. */
    fstream* fStreamPtr;

    /** The flag indicates that the file is opening to read only. */
    bool isReadOnly;
};

#endif	/* FILE_H */


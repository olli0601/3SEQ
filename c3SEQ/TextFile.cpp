/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   TextFile.cpp
 * CREATED ON:  15 June 2011, 17:11
 * AUTHOR:      Ha Minh Lam
 * 
 * DESCRIPTION: Implement methods to handle text files easier.
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-06-15   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   15 June 2011
 */

#include "TextFile.h"
#include "Interface.h"

////////////////////////////////////////////////////////////////////////////////
// Static definition
////////////////////////////////////////////////////////////////////////////////

const char TextFile::INVALID_PATH_CHARS[] = {'\"', '*', '<', '>', '?', '|'};
const unsigned long TextFile::INVALID_CHARS_COUNT
        = sizeof (INVALID_PATH_CHARS) / sizeof (char);

const char TextFile::PATH_SEPARATORS[] = {'\\', '/'};
const unsigned long TextFile::PATH_SEPARATORS_COUNT
        = sizeof (PATH_SEPARATORS) / sizeof (char);

const char TextFile::INVALID_PATH_REPLACEMENT = '_';

const string TextFile::DEFAULT_COLUMN_DELIM = "\t";

////////////////////////////////////////////////////////////////////////////////

TextFile::TextFile(string newFilePath) : filePath(newFilePath) {
    fStreamPtr = new fstream;
    isReadOnly = false;

    /* Test file path */
    try {
        if (filePath.length() <= 0) throw (filePath);

        for (unsigned long i = 0; i < filePath.length(); i++) {
            if (isInvalidPathChar(filePath[i])) {
                throw (filePath);
            }
        }

    } catch (string invalidFilePath) {
        Interface::instance() << "Invalid file name: \"" << invalidFilePath << "\".\n";
        Interface::instance().showError(true, true); // show error & exit
    }
}

TextFile::TextFile(const TextFile& orig)
: filePath(orig.filePath), fStreamPtr(orig.fStreamPtr), isReadOnly(orig.isReadOnly) {
}

TextFile& TextFile::operator=(const TextFile& rhs) {
    if (this != &rhs) {
        this->filePath = rhs.filePath;
        this->fStreamPtr = rhs.fStreamPtr;
        this->isReadOnly = rhs.isReadOnly;
    }

    return *this;
}

TextFile::~TextFile() {
    if (fStreamPtr) {
        delete fStreamPtr;
        fStreamPtr = NULL;
    }
}

string TextFile::getLine(void) {
    assert(isReadOnly && fStreamPtr->is_open());

    string newLine("");

    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.
    std::istream::sentry se(*fStreamPtr);
    
    std::streambuf* buffer = fStreamPtr->rdbuf();

    while (true) {
        int newChar = buffer->sbumpc();
        
        switch (newChar) {
            case '\r':
                newChar = buffer->sgetc();
                if (newChar == '\n')
                    buffer->sbumpc();
                return newLine;
                
            case '\n':
            case EOF:
                return newLine;
                
            default:
                newLine += (char) newChar;
        }
    }
}

void TextFile::openToRead(void) {
    if (fStreamPtr->is_open()) {
        fStreamPtr->close();
    }

    fStreamPtr->open(filePath.c_str(), ios::in);
    isReadOnly = true;
}

void TextFile::openToWrite(void) {
    if (fStreamPtr->is_open()) {
        fStreamPtr->close();
    }

    fStreamPtr->open(filePath.c_str(), ios::out);
    isReadOnly = false;
}

void TextFile::close(void) {
    if (fStreamPtr->is_open()) {
        fStreamPtr->close();
    }
    isReadOnly = false;
}

const bool TextFile::exists(void) {
    if (filePath.length() <= 0) {
        return false;
    }

    openToRead();

    if (isStreamReadable()) {
        close();
        return true;
    } else {
        return false;
    }
}

const bool TextFile::isStreamReadable(void) const {
    if (!fStreamPtr->is_open() || fStreamPtr->fail() || fStreamPtr->eof()) {
        return false;
    }
    return true;
}

vector<string> TextFile::readAllLines(void) {
    vector<string> lines;
    string newLine("");

    assert(exists());
    openToRead();

    /* Read through each line of the file */
    while (isStreamReadable()) {
        newLine = getLine();
        newLine = String::trim(newLine);

        if (newLine.length() > 0) {
            lines.push_back(newLine);
        }
    }

    close();
    return lines;
}
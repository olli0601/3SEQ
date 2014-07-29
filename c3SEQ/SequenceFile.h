/* 
 * File:   SequenceFile.h
 * Author: Ha Minh Lam
 *
 * Created on 03 August 2011, 15:33
 */

#ifndef SEQUENCEFILE_H
#define	SEQUENCEFILE_H

#include <string>
#include <iostream>
#include "TextFile.h"
#include "GenomeSet.h"
#include "Segment.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  FASTA FILE
////////////////////////////////////////////////////////////////////////////////

class FastaFile : public TextFile {
public:
    explicit FastaFile(string newFilePath);

    virtual ~FastaFile();

    virtual void read(GenomeSet* storage, string newSegmentName);

    virtual void write(GenomeSet* genomeSet);
    
    static void write(GenomeSet* genomeSet, ostream* streamPtr);
};


////////////////////////////////////////////////////////////////////////////////
//  PHYLIP FILE
////////////////////////////////////////////////////////////////////////////////

class PhylipFile : public TextFile {
public:
    explicit PhylipFile(string newFilePath);

    virtual ~PhylipFile();

    virtual void read(GenomeSet* storage, string newSegmentName);

    virtual void write(GenomeSet* genomeSet);
    
    static void write(GenomeSet* genomeSet, ostream* streamPtr);
};


////////////////////////////////////////////////////////////////////////////////
//  NEXUS FILE
////////////////////////////////////////////////////////////////////////////////

class NexusFile : public TextFile {
public:
    explicit NexusFile(string newFilePath);

    virtual ~NexusFile();

    virtual void read(GenomeSet* storage, string newSegmentName);

    virtual void write(GenomeSet* genomeSet);
    
    static void write(GenomeSet* genomeSet, ostream* streamPtr);
};


////////////////////////////////////////////////////////////////////////////////
//  SEQUENCE FILE
////////////////////////////////////////////////////////////////////////////////

class SequenceFile : public FastaFile, public PhylipFile, public NexusFile {
public:

    enum Type {
        UNKNOWN, PHYLIP, NEXUS, FASTA
    };

    explicit SequenceFile(string newFilePath, Type newFileType);

    virtual ~SequenceFile();

    virtual const bool exists(void) {
        return FastaFile::exists();
    };

    virtual const string& getPath(void) {
        return FastaFile::getPath();
    };
    
    virtual void setFileType(Type newType) {
        assert(fileType == UNKNOWN);
        fileType = newType;
    };

    virtual void detectFileType(void);

    virtual void read(GenomeSet* storage, string newSegmentName);

    virtual void write(GenomeSet* genomeSet);


private:
    Type fileType;
};



#endif	/* SEQUENCEFILE_H */


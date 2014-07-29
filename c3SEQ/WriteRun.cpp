/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   WriteRun.cpp
 * CREATED ON:  19 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-19  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   19 August 2011
 */

#include "WriteRun.h"

WriteRun::WriteRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Write Mode");
    outputFile = NULL;
}

WriteRun::WriteRun(const WriteRun& orig) {
    assert(false); // should never reach here
}

WriteRun& WriteRun::operator=(const WriteRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

WriteRun::~WriteRun() {
    if (outputFile) delete outputFile;
}

void WriteRun::parseCmdLine() {
    int argNum = getRunArgsNum();

    if (argNum < 1) {
        Interface::instance() << "Not enough parameter for write-mode.\n";
        Interface::instance().showError(true, true);
    }

    string inFilePath = argVector[2];
    SequenceFile inputFile(inFilePath, SequenceFile::UNKNOWN);

    if (argNum >= 2) {
        string outFilePath = argVector[3];
        outputFile = new SequenceFile(outFilePath, SequenceFile::UNKNOWN);
    }

    if (!inputFile.exists()) {
        Interface::instance() << "File \"" << inFilePath << "\" not found.\n";
        Interface::instance().showError(true, true);
    }

    inputFile.detectFileType();
    parentDataset = new GenomeSet();
    inputFile.read(parentDataset, "");

    deleteCmdLineArgs(argNum + 1); // --write  inputFile  outputFile
    Run::parseCmdLine(); // Get common run options
}

void WriteRun::perform() {
    Run::perform(); // process common data

    if (cloOutputFileType == SequenceFile::UNKNOWN) {
        cloOutputFileType = SequenceFile::PHYLIP;
    }

    if (outputFile != NULL) {
        outputFile->setFileType(cloOutputFileType);
        outputFile->write(childDataset);

        Interface::instance()
                << "All extracted sequences have been stored into file\n"
                << Interface::DEFAULT_INDENT << "\"" << outputFile->getPath() << "\"\n";
        Interface::instance().showLog(true);
        
    } else {
        switch (cloOutputFileType) {
            case SequenceFile::FASTA:
                FastaFile::write(childDataset, &std::cout);
                break;
                
            case SequenceFile::NEXUS:
                NexusFile::write(childDataset, &std::cout);
                break;
                
            case SequenceFile::PHYLIP:  // fall through
            case SequenceFile::UNKNOWN: // fall through
            default:
                PhylipFile::write(childDataset, &std::cout);
        }
    }
}
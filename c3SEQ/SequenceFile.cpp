/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   SequenceFile.cpp
 * CREATED ON:  03 August 2011, 15:33
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: Implement methods to read/write data from/into sequence files.
 *              Currently support Fasta, Phylip, and Nexus files.
 *              Multiple segments in one file have not been supported yet.
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-08-03   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   03 August 2011
 */

#include "SequenceFile.h"

SequenceFile::SequenceFile(string newFilePath, Type newFileType)
: FastaFile(newFilePath), PhylipFile(newFilePath), NexusFile(newFilePath), fileType(newFileType) {
}

SequenceFile::~SequenceFile() {
}

void SequenceFile::detectFileType(void) {
    /* This method should not be called if the file type has been known */
    assert(fileType == UNKNOWN);

    string firstString("");

    /* Let assume this is a Fasta file to get rid 
     * of the complexity of the multiple inheritance. */
    assert(FastaFile::exists());
    FastaFile::openToRead();

    fstream* fileStream = FastaFile::getStreamPtr();
    (*fileStream) >> firstString;

    if (FastaFile::isStreamReadable()) {
        if (firstString == "#NEXUS") {
            fileType = NEXUS;
        } else if (firstString[0] == '>') {
            fileType = FASTA;
        } else {
            string secondString;
            (*fileStream) >> secondString;
            if (FastaFile::isStreamReadable()
                    && String::isInteger(firstString)
                    && String::isInteger(secondString)) {
                fileType = PHYLIP;
            }
        }
    }

    FastaFile::close();
}

void SequenceFile::read(GenomeSet* storage, string newSegmentName) {
    switch (fileType) {
        case FASTA:
            FastaFile::read(storage, newSegmentName);
            return;

        case PHYLIP:
            return PhylipFile::read(storage, newSegmentName);

        case NEXUS:
            return NexusFile::read(storage, newSegmentName);

        default:
            Interface::instance() << "Cannot detect sequence file format.\n";
            Interface::instance().showError(true, true); // show error and exit
    }

    assert(false);
}

void SequenceFile::write(GenomeSet* genomeSet) {
    switch (fileType) {
        case FASTA:
            FastaFile::write(genomeSet);
            return;

        case PHYLIP:
            PhylipFile::write(genomeSet);
            return;

        case NEXUS:
            NexusFile::write(genomeSet);
            return;

        default:
            Interface::instance() << "Cannot detect output file format.\n";
            Interface::instance().showError(true, true); // show error and exit
            return;
    }
}



////////////////////////////////////////////////////////////////////////////////
//  FASTA FILE
////////////////////////////////////////////////////////////////////////////////

FastaFile::FastaFile(string newFilePath) : TextFile(newFilePath) {

}

FastaFile::~FastaFile() {

}

void FastaFile::read(GenomeSet* storage, string newSegmentName) {
    unsigned long lastSegmentNum = 0;
    if (storage->getSize() > 0) {
        lastSegmentNum = storage->getGenome(0)->getSegmentNum();
    }

    if (newSegmentName.length() <= 0) {
        stringstream tmpStream;
        tmpStream << (lastSegmentNum + 1);
        newSegmentName = tmpStream.str();
    }

    string accession("");
    string sequenceStr("");
    string line("");

    unsigned long validSeqCount = 0;
    unsigned long longestSeq = 0;
    unsigned long shortestSeq = numeric_limits<unsigned long>::max();

    vector<string> notFoundList;

    assert(exists());
    openToRead();

    /* Read through each line of the file */
    while (isStreamReadable()) {
        line = getLine();
        line = String::trim(line);

        /* If reach end of file or a new accession number is found */
        if (!isStreamReadable() || (line.length() > 0 && line[0] == '>')) {
            if (line.length() > 0 && line[0] != '>') {
                sequenceStr += String::deleteAllSpace(line);
            }

            /* If the previous sequence has not been stored into the set */
            if (sequenceStr.length() > 0) {
                /* Check if its accession is valid */
                if (accession.length() <= 0) {
                    delete storage;
                    Interface::instance() << "Invalid accession at sequence "
                            << validSeqCount + 1 << " in file " << getPath() << ".\n";
                    Interface::instance().showError(true, true);
                }

                unsigned long seqLength = sequenceStr.length();
                if (seqLength > longestSeq) {
                    longestSeq = seqLength;
                }
                if (seqLength < shortestSeq) {
                    shortestSeq = seqLength;
                }

                Segment* segment = new Segment(sequenceStr, newSegmentName);

                /* Try to find the genome */
                Genome* genome = storage->getGenome(accession);

                if (genome && genome->getSegmentNum() > lastSegmentNum) {
                    /* This means duplicate accession found */
                    Interface::instance() << "Duplicate accession number '"
                            << genome->getAccession() << "' in file '"
                            << getPath() << "'.\n";
                    Interface::instance().showError(true, true); // show error & exit
                }

                if (lastSegmentNum == 0) {
                    genome = new Genome(accession);
                    genome->addSegment(segment);
                    storage->addGenome(genome);
                    validSeqCount++;

                } else {
                    if (genome) {
                        /* Found */
                        genome->addSegment(segment);
                        validSeqCount++;
                    } else {
                        /* Not found */
                        notFoundList.push_back(accession);
                    }
                }
            }

            if (line.length() > 0 && line[0] == '>') {
                accession = line.substr(1, line.length() - 1);
                sequenceStr = "";
            }

        } else {
            sequenceStr += String::deleteAllSpace(line);
        }
    };

    close();
    assert(validSeqCount <= storage->getSize());


    /* Verify that the sequence in this file is aligned */
    if (storage->getSize() > 0 && longestSeq != shortestSeq) {
        delete storage;
        Interface::instance() << "Fasta file is not aligned. Sequences range from "
                << shortestSeq << "nt to " << longestSeq << "nt.\n";
        Interface::instance().showError(true, true); // show error and exit
    }


    /* Show warnings if there is any unfound sequence */
    if (notFoundList.size() > 0) {
        Interface::instance() << "Unable to find the following " << notFoundList.size()
                << " sequences to add segment " << newSegmentName << " :" << endl;
        for (unsigned long i = 0; i < notFoundList.size(); i++) {
            Interface::instance() << "\t" << notFoundList[i] << endl;
        }
        Interface::instance() << endl
                << "Do you want to ignore these sequences and continue the analysis?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);
        }
    }


    /* Show warning if there is lack of segment data */
    if (validSeqCount < storage->getSize()) {
        Interface::instance() << "No data of segment " << newSegmentName
                << " for the following sequences:" << endl;
        for (unsigned long i = 0; i < storage->getSize(); i++) {
            if (storage->getGenome(i)->getSegmentNum() < lastSegmentNum + 1) {
                Genome* curGenome = storage->removeGenome(i);
                Interface::instance() << "\t" << curGenome->getAccession() << endl;
                delete curGenome;
            }
        }
        Interface::instance() << endl
                << "Do you want to ignore these sequences and continue the analysis?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);
        }
    }
}

void FastaFile::write(GenomeSet* genomeSet) {
    if (exists()) {
        Interface::instance()
                << "The file \"" << getPath() << "\" has already existed.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);

        } else {
            removeFile();
        }
    }

    openToWrite();
    fstream* fstreamPtr = getStreamPtr();
    write(genomeSet, fstreamPtr);
    close();
}

void FastaFile::write(GenomeSet* genomeSet, ostream* streamPtr) {
    for (unsigned long i = 0; i < genomeSet->getSize(); i++) {
        (*streamPtr)
                << ">" << genomeSet->getGenome(i)->getAccession() << endl
                << genomeSet->getGenome(i)->toString() << endl;
    }
    
    (*streamPtr) << endl;
}



////////////////////////////////////////////////////////////////////////////////
//  PHYLIP FILE
////////////////////////////////////////////////////////////////////////////////

PhylipFile::PhylipFile(string newFilePath) : TextFile(newFilePath) {

}

PhylipFile::~PhylipFile() {

}

void PhylipFile::read(GenomeSet* storage, string newSegmentName) {
    //FIXME: cannot read Phylip interleaved files.

    unsigned long lastSegmentNum = 0;
    if (storage->getSize() > 0) {
        lastSegmentNum = storage->getGenome(0)->getSegmentNum();
    }

    if (newSegmentName.length() <= 0) {
        stringstream tmpStream;
        tmpStream << (lastSegmentNum + 1);
        newSegmentName = tmpStream.str();
    }

    string accession("");
    string segmentStr("");

    unsigned long seqNum = 0;
    unsigned long validSeqCount = 0;
    unsigned long seqLength = 0;
    vector<string> notFoundList;

    try {
        assert(exists());
        openToRead();
        fstream* fStreamPtr = getStreamPtr();

        (*fStreamPtr) >> seqNum;
        (*fStreamPtr) >> seqLength;

        /* Now read in accession numbers and nt-sequences seqNum times */
        for (unsigned long i = 0L; i < seqNum; i++) {
            (*fStreamPtr) >> accession;
            (*fStreamPtr) >> segmentStr;

            if (segmentStr.length() != seqLength) {
                delete storage;
                Interface::instance() << "Sequence " << accession
                        << " has length " << segmentStr.length()
                        << ".  Should be " << seqLength << ".\n";
                Interface::instance().showError(true, true); // show error & exit
            }

            Segment* segment = new Segment(segmentStr, newSegmentName);

            /* Try to find the genome */
            Genome* genome = storage->getGenome(accession);

            if (genome && genome->getSegmentNum() > lastSegmentNum) {
                /* This means duplicate accession found */
                Interface::instance() << "Duplicate accession number \""
                        << genome->getAccession() << "\" in file \""
                        << getPath() << "\".\n";
                Interface::instance().showError(true, true); // show error & exit
            }

            if (lastSegmentNum == 0) {
                genome = new Genome(accession);
                genome->addSegment(segment);
                storage->addGenome(genome);
                validSeqCount++;

            } else {
                if (genome) {
                    /* Found */
                    genome->addSegment(segment);
                    validSeqCount++;
                } else {
                    /* Not found */
                    notFoundList.push_back(accession);
                }
            }
        }

        close();
        assert(validSeqCount <= storage->getSize());

    } catch (...) {
        Interface::instance() << "An error occurs during reading progress. "
                << "Please, check your file: \"" << getPath() << "\".\n";
        Interface::instance().showError(true, true);
    }


    /* Show warnings if there is any unfound sequence */
    if (notFoundList.size() > 0) {
        Interface::instance() << "Unable to find the following " << notFoundList.size()
                << " sequences to add segment " << newSegmentName << " :" << endl;
        for (unsigned long i = 0; i < notFoundList.size(); i++) {
            Interface::instance() << "\t" << notFoundList[i] << endl;
        }
        Interface::instance() << endl
                << "Do you want to ignore these sequences and continue the analysis?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);
        }
    }


    /* Show warning if there is lack of segment data */
    if (validSeqCount < storage->getSize()) {
        Interface::instance() << "No data of segment " << newSegmentName
                << " for the following sequences:" << endl;
        for (unsigned long i = 0; i < storage->getSize(); i++) {
            if (storage->getGenome(i)->getSegmentNum() < lastSegmentNum + 1) {
                Genome* curGenome = storage->removeGenome(i);
                Interface::instance() << "\t" << curGenome->getAccession() << endl;
                delete curGenome;
            }
        }
        Interface::instance() << endl
                << "Do you want to ignore these sequences and continue the analysis?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);
        }
    }
}

void PhylipFile::write(GenomeSet* genomeSet) {
    if (exists()) {
        Interface::instance()
                << "The file \"" << getPath() << "\" has already existed.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);

        } else {
            removeFile();
        }
    }

    openToWrite();
    fstream* fstreamPtr = getStreamPtr();
    write(genomeSet, fstreamPtr);
    close();
}

void PhylipFile::write(GenomeSet* genomeSet, ostream* streamPtr) {
    (*streamPtr)
            << "  " << genomeSet->getSize()
            << "  " << genomeSet->getExtractedLength() << endl;

    for (unsigned long i = 0; i < genomeSet->getSize(); i++) {
        (*streamPtr)
                << genomeSet->getGenome(i)->getAccession() << "\t"
                << genomeSet->getGenome(i)->toString() << endl;
    }
    
    (*streamPtr) << endl;
}



////////////////////////////////////////////////////////////////////////////////
//  NEXUS FILE
////////////////////////////////////////////////////////////////////////////////

NexusFile::NexusFile(string newFilePath) : TextFile(newFilePath) {

}

NexusFile::~NexusFile() {

}

void NexusFile::read(GenomeSet* storage, string newSegmentName) {
    //TODO: Read Nexus file
    Interface::instance() << "Reading NEXUS file has not been supported yet.\n";
    Interface::instance().showError(true, true); // show error and exit
}

void NexusFile::write(GenomeSet* genomeSet) {
    if (exists()) {
        Interface::instance()
                << "The file \"" << getPath() << "\" has already existed.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);

        } else {
            removeFile();
        }
    }

    openToWrite();
    fstream* fstreamPtr = getStreamPtr();
    write(genomeSet, fstreamPtr);
    close();
}

void NexusFile::write(GenomeSet* genomeSet, ostream* streamPtr) {
    (*streamPtr)
            << "#NEXUS\n\n[!Data from\n\n]\n\nbegin taxa;\n\tdimensions ntax="
            << genomeSet->getSize()
            << ";\n\ttaxlabels";

    for (unsigned long i = 0; i < genomeSet->getSize(); i++) {
        (*streamPtr) << "\n\t\t" << genomeSet->getGenome(i)->getAccession();
    }

    (*streamPtr)
            << ";\nend;\n\nbegin characters;\n\tdimensions nchar=" << genomeSet->getExtractedLength()
            << ";\n\tformat missing=? gap=- matchchar=. datatype=dna;\n\toptions gapmode=missing;\n\tmatrix\n";

    for (unsigned long i = 0; i < genomeSet->getSize(); i++) {
        (*streamPtr)
                << "\n" << genomeSet->getGenome(i)->getAccession()
                << "    \t" << genomeSet->getGenome(i)->toString();
    }

    (*streamPtr) << ";\nend;\n\n";
}
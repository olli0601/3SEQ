/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   TripletRun.cpp
 * CREATED ON:  25 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "TripletRun.h"
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-25  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   25 August 2011
 */

#include "TripletRun.h"
#include "Triplet.h"

TripletRun::TripletRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Triplet Run");
}

TripletRun::TripletRun(const TripletRun& orig) {
    assert(false); // should never reach here
}

TripletRun& TripletRun::operator=(const TripletRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

TripletRun::~TripletRun() {
}

void TripletRun::parseCmdLine() {
    int argNum = getRunArgsNum();

    if (argNum < 4) {
        Interface::instance() << "Not enough parameter for triplet run.\n";
        Interface::instance().showError(true, true);
    }
    argNum = 4;

    string filePath = argVector[2];
    SequenceFile sequenceFile(filePath, SequenceFile::UNKNOWN);
    if (!sequenceFile.exists()) {
        Interface::instance() << "File \"" << filePath << "\" not found.\n";
        Interface::instance().showError(true, true);
    }
    sequenceFile.detectFileType();
    parentDataset = new GenomeSet();
    sequenceFile.read(parentDataset, "");

    for (int i = 0; i < SEQ_NUM; i++) {
        accessionList[i] = argVector[3 + i];
    }

    deleteCmdLineArgs(argNum + 1);
    Run::parseCmdLine(); // Get common run options
}

void TripletRun::perform() {
    Run::perform(); // process common data

    /* The 3 sequences in the triplet run:
     *     [0] First parent
     *     [1] Second parent
     *     [2] Child
     */
    Genome * ptrGenomes[SEQ_NUM];

    for (int i = 0; i < SEQ_NUM; i++) {
        ptrGenomes[i] = parentDataset->getGenome(accessionList[i]);

        if (ptrGenomes[i] == NULL) {
            Interface::instance() << "No sequence with accession number \"" << accessionList[i] << "\".\n";
            Interface::instance().showError(true, true);
        }
    }

    Genome* dad = ptrGenomes[0];
    Genome* mum = ptrGenomes[1];
    Genome* child = ptrGenomes[2];

    Triplet triplet(dad, mum, child);

    Interface::instance()
            << "Child accession :       " << child->getAccession() << "\n"
            << endl
            << "Parent P\n"
            << "   Accession :          " << dad->getAccession() << "\n"
            << "   Distance to child :  " << dad->ntDistanceTo(child, false)
            << "nt  (" << dad->ntDistanceTo(child, true) << "nt ignoring gaps)\n"
            << endl
            << "Parent Q\n"
            << "   Accession :          " << mum->getAccession() << "\n"
            << "   Distance to child :  " << mum->ntDistanceTo(child, false)
            << "nt  (" << mum->ntDistanceTo(child, true) << "nt ignoring gaps)\n";
    Interface::instance().showLog(true);
    
    Interface::instance()
            << "Informative sites of type P :  " << triplet.getUpStep() << "\n"
            << "Informative sites of type Q :  " << triplet.getDownStep() << "\n"
            << "Maximum Descent :              " << triplet.getMaxDecent() << "\n";
    Interface::instance().showOutput(true);

    vector<Triplet::BreakPointPair> bpPairList = triplet.getBreakPointPairs();
    Interface::instance() << "Breakpoint Ranges:\n";

    double totalBpPair = 0.0;

    for (unsigned long i = 0; i < bpPairList.size(); i++) {
        Triplet::BreakPointPair bpPair = bpPairList[i];
        Interface::instance()
                << Interface::DEFAULT_INDENT << bpPair.toString() << "\n";

        totalBpPair +=
                static_cast<double> (bpPair.firstBp->getRightBoundPos()
                - bpPair.firstBp->getLeftBoundPos() + 1)
                *
                static_cast<double> (bpPair.secondBp->getRightBoundPos()
                - bpPair.secondBp->getLeftBoundPos() + 1);
    }

    Interface::instance()
            << endl
            << "# Breakpoint Pairs: " << totalBpPair << "\n"
            << "# Breakpoint Ranges: " << bpPairList.size() << "\n"
            << endl
            << "Minimum length of recombinant segment is " << triplet.getMinRecLength()
            << "nt \n(excluding gaps and ambiguous nucleotides)\n";
    Interface::instance().showOutput(true);

    stringstream nameStream;
    nameStream << parentDataset->getIndex(dad) + 1 << "_"
            << parentDataset->getIndex(mum) + 1 << "_"
            << parentDataset->getIndex(child) + 1;

    string fileName = nameStream.str();
    if (id.length() > 0) {
        fileName += "." + id;
    }
    fileName += "." + DEFAULT_TRIPLET_PS_FILE_EXT;

    triplet.generatePSFile(fileName);
}

/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   InfoRun.cpp
 * CREATED ON:  18 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION:
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-18  created
 *
 * VERSION:     1.0
 * LAST EDIT:   18 August 2011
 */

#include "InfoRun.h"
#include <cmath>

InfoRun::InfoRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Info Mode");
}

InfoRun::InfoRun(const InfoRun& orig) {
    assert(false); // should never reach here
}

InfoRun& InfoRun::operator=(const InfoRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }

    return *this;
}

InfoRun::~InfoRun() {
}

void InfoRun::parseCmdLine() {
    if (getRunArgsNum() < 1) {
        Interface::instance() << "Not enough parameter for info-run mode.\n";
        Interface::instance().showError(true, true);
    }

    string filePath = argVector[2];

    SequenceFile inputFile(filePath, SequenceFile::UNKNOWN);
    if (!inputFile.exists()) {
        Interface::instance() << "File \"" << filePath << "\" not found.\n";
        Interface::instance().showError(true, true);

    }

    inputFile.detectFileType();
    parentDataset = new GenomeSet();
    inputFile.read(parentDataset, "");

    deleteCmdLineArgs(2); // --info inputFile
    Run::parseCmdLine(); // Get common run options
}

void InfoRun::perform() {
    Run::perform(); // process common data


    /* Declare all variable as double-type to avoid overflowing and
     * to make further calculation more convenient */
    double parentNum = static_cast<double> (parentDataset->getSize());

    double extractedLength = static_cast<double> (parentDataset->getExtractedLength());

    double activeLength = static_cast<double> (cloLastNucleotidePos - cloFirstNucleotidePos);
    if (cloIsFirstToLastCutOut) {
        activeLength = static_cast<double> (parentDataset->getOriginalLength()) - activeLength;
    }

    const GenomeSet::DistanceStats distanceStats = parentDataset->getDistanceStats();


    Interface::instance()
            << "Number of sequences          :  " << parentNum << "\n"
            << "Number of distinct sequences :  " << parentDataset->getDistinctGenomes(false) << "\n"
            << "Number of comparisons        :  " << getNumTotalTriplets() << "\n";
    Interface::instance().showOutput(true);

    Interface::instance()
            << "Original sequence length is " << parentDataset->getOriginalLength() << " nt.\n";

    if (cloIsFirstToLastCutOut) {
        Interface::instance()
                << "Cutting out " << cloLastNucleotidePos - cloFirstNucleotidePos
                << " positions, from " << cloFirstNucleotidePos + 1
                << " to " << cloLastNucleotidePos << " inclusive.\n";

    } else {
        Interface::instance()
                << "Using " << cloLastNucleotidePos - cloFirstNucleotidePos
                << " positions, from " << cloFirstNucleotidePos + 1
                << " to " << cloLastNucleotidePos << " inclusive.\n";
    }

    if (cloUseAllSites) {
        Interface::instance()
                << endl
                << "Using all " << extractedLength << " sites.\n"
                << "Counting gap/nucleotide mismatches as differences.\n"
                << "    gapped: " << parentDataset->getGapNum()
                << " - mono: " << parentDataset->getMonoallelicNum()
                << " - bi: " << parentDataset->getBiallelicNum()
                << " - tri: " << parentDataset->getTriallelicNum()
                << " - tetra: " << parentDataset->getTetrallelicNum() << "\n";

    } else {
        Interface::instance()
                << endl
                << "Number of polymorphic sites is " << extractedLength
                << ". Using only polymorphic sites:\n"
                << "    bi: " << parentDataset->getBiallelicNum()
                << " - tri: " << parentDataset->getTriallelicNum()
                << " - tetra: " << parentDataset->getTetrallelicNum() << "\n";
    }

    string strSites = "sites";
    if (!cloUseAllSites) {
        strSites = "polymorphic sites";
    }
    Interface::instance()
            << endl
            << "Number of " << strSites << " where there are no gaps :  "
            << parentDataset->getNoGapSitesNum() << "/" << extractedLength << "\n";

    if (!cloUseAllSites) {
        Interface::instance()
                << endl
                << "Watterson's Theta            :  " << extractedLength / log(parentNum) << "\n"
                << "Watterson's Theta (per site) :  " << extractedLength / (log(parentNum) * activeLength) << "\n";
    }

    Interface::instance()
            << endl
            << "Min  pairwise-distance between strains :  " << distanceStats.min << "\n"
            << "Max  pairwise-distance between strains :  " << distanceStats.max << "\n"
            << "Mean pairwise-distance between strains :  " << distanceStats.mean << "\n"
            << "Mean pairwise-distance per site (pi)   :  " << distanceStats.mean / activeLength << "\n";

    Interface::instance().showOutput(true);
}
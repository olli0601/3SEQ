/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   MatchRun.cpp
 * CREATED ON:  25 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "MatchRun.h"
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-25  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   25 August 2011
 */

#include "MatchRun.h"
#include <cstdlib>

MatchRun::MatchRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Match Run");

    queryAccession = "";
    minMatchDistance = 0;
    maxMatchDistance = 0;
}

MatchRun::MatchRun(const MatchRun& orig) {
    assert(false); // should never reach here
}

MatchRun& MatchRun::operator=(const MatchRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

MatchRun::~MatchRun() {
}

void MatchRun::parseCmdLine() {
    int argNum = getRunArgsNum();

    if (argNum < 4) {
        Interface::instance() << "Not enough parameter for match run.\n";
        Interface::instance().showError(true, true);
    }
    argNum = 4;

    string filePath = argVector[2];
    queryAccession = argVector[3];
    minMatchDistance = strtoul(argVector[4], NULL, 0);
    maxMatchDistance = strtoul(argVector[5], NULL, 0);

    SequenceFile sequenceFile(filePath, SequenceFile::UNKNOWN);
    if (!sequenceFile.exists()) {
        Interface::instance() << "File \"" << filePath << "\" not found.\n";
        Interface::instance().showError(true, true);
    }
    sequenceFile.detectFileType();
    parentDataset = new GenomeSet();
    sequenceFile.read(parentDataset, "");

    if (minMatchDistance > maxMatchDistance) {
        Interface::instance() << "The minimum distance cannot be greater than the maximum distance.\n";
        Interface::instance().showError(true, true);
    }

    deleteCmdLineArgs(argNum + 1);
    Run::parseCmdLine(); // Get common run options
}

void MatchRun::perform() {
    Run::perform(); // process common data
    
    Genome* querySeq = parentDataset->getGenome(queryAccession);

    if (querySeq == NULL) {
        Interface::instance() << "The sequence \"" << queryAccession << "\" does not exist in the given dataset.\n";
        Interface::instance().showError(true, true);
    }

    Interface::instance().initCounter("Analysing dataset", 0, parentDataset->getSize(), true);

    for (unsigned long i = 0; i < parentDataset->getSize(); i++) {
        Interface::instance().count(i, false);

        MatchResult currentResult;
        currentResult.genome = parentDataset->getGenome(i);
        if (querySeq == currentResult.genome) continue;

        currentResult.distance = querySeq->ntDistanceTo(currentResult.genome, true);
        if (currentResult.distance >= minMatchDistance && currentResult.distance <= maxMatchDistance) {
            results.push_back(currentResult);
        }
    }

    Interface::instance().finishCounting();

    displayOutput();
    saveMatchIntoFile();
}

void MatchRun::displayOutput(void) const {
    if (results.size() <= 0) {
        Interface::instance() << "No matched sequence found.\n";
    } else {
        Interface::instance() << "The following " << results.size()
                << " sequences are between " << minMatchDistance
                << " and " << maxMatchDistance << " nt different from the query\n"
                << "sequence \"" << queryAccession << "\":\n";
    }

    Interface::instance().showLog(true);

    char formatedDist[50];

    for (unsigned long i = 0; i < results.size(); i++) {
        sprintf(formatedDist, "%4lu", results[i].distance);
        Interface::instance() << formatedDist << "\t" << results[i].genome->getAccession() << "\n";
    }
    
    Interface::instance().showOutput(true);
}

void MatchRun::saveMatchIntoFile(void) const {
    if (results.size() <= 0) return;

    string fileName = DEFAULT_MATCH_FILE_NAME;
    if (id.length() > 0) {
        fileName = id + fileName;
    }

    TextFile matchFile(fileName);
    if (matchFile.exists()) {
        Interface::instance()
                << "The file \"" << matchFile.getPath() << "\" already exists.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* If the user say "No" */
            return;
            
        } else {
            matchFile.removeFile();
        }
    }

    for (unsigned long i = 0; i < results.size(); i++) {
        matchFile.writeLine(results[i].genome->getAccession());
    }

    Interface::instance() << "All matched sequences are stored in the file \""
            << matchFile.getPath() << "\".\n"
            << "This file can be used as a subset file for further analysis.\n";
    Interface::instance().showLog(true);
}

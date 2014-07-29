/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   GenPRun.cpp
 * CREATED ON:  19 August 2011
 * AUTHOR:      Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-19  created
 *
 * VERSION:     1.0
 * LAST EDIT:   19 August 2011
 */

#include "GenPRun.h"
#include "PTable.h"
#include <cstdlib>

GenPRun::GenPRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Generate P-Value Table");
}

GenPRun::GenPRun(const GenPRun& orig) {
    assert(false); // should never reach here
}

GenPRun& GenPRun::operator=(const GenPRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

GenPRun::~GenPRun() {
}

void GenPRun::parseCmdLine() {
    if (getRunArgsNum() < 2) {
        Interface::instance() << "Not enough parameter to generate P-value table.\n";
        Interface::instance().showError(true, true);
    }

    pTableFilePath = argVector[2];
    pTableSize = 0;
    pTableSize = atoi(argVector[3]);

    if (pTableSize < 2) {
        Interface::instance() << "Invalid table size.\n";
        Interface::instance().showError(true, true);
    }

    TextFile testFile(pTableFilePath);
    if (testFile.exists()) {
        Interface::instance()
                << "The file \"" << pTableFilePath << "\" already exists.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* The program will be exit if the user say "No" */
            Interface::instance().throwExitSignal(0);
        } else {
            testFile.removeFile();
        }
    }
}

void GenPRun::perform() {
    PTable::instance().generateTable(pTableSize, pTableSize, pTableSize);

    if (PTable::instance().saveToFile(pTableFilePath)) {
        Interface::instance()
                << "The new P-value table has been stored into file: \"" << pTableFilePath << "\".\n";
        Interface::instance().showLog(true);
    } else {
        Interface::instance()
                << "The new P-value table cannot be stored into file: \"" << pTableFilePath << "\".\n"
                << "An error occurred during saving progress.\n";
        Interface::instance().showError(true, true);
    }
}
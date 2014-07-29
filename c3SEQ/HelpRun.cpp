/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   HelpRun.cpp
 * CREATED ON:  12 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-12  created
 *
 * VERSION:     1.0
 * LAST EDIT:   12 August 2011
 */

#include "HelpRun.h"

HelpRun::HelpRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram(DEFAULT_PROGRAM_DESC);
}

HelpRun::HelpRun(const HelpRun& orig) {
    assert(false); // should never copy
}

HelpRun& HelpRun::operator=(const HelpRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never assign
    }

    return *this;
}

HelpRun::~HelpRun() {
}

void HelpRun::parseCmdLine() {
    // Do nothing
}

void HelpRun::perform() {
    Interface::instance()
            << "For help, email author: " << AUTHOR_EMAIL << "\n"
            << endl
            << "Input file may be PHYLIP (sequential) or FASTA format.\n"
            << endl
            << Interface::SEPARATOR << endl
            << endl;

    Interface::instance().catchUsageModes();
    Interface::instance() << endl
            << Interface::SEPARATOR << endl
            << endl;
    Interface::instance().catchOptions();

    Interface::instance().showLog(true);
}
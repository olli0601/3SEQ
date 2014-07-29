/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   VersionRun.cpp
 * CREATED ON:  17 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-17  created
 *
 * VERSION:     1.0
 * LAST EDIT:   17 August 2011
 */

#include "VersionRun.h"
#include "config.h"

VersionRun::VersionRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram(DEFAULT_PROGRAM_DESC);
}

VersionRun::VersionRun(const VersionRun& orig) {
    assert(false); // never copy
}

VersionRun& VersionRun::operator=(const VersionRun& rhs) {
    if (this != &rhs) {
        assert(false); // never assign
    }

    return *this;
}

VersionRun::~VersionRun() {
}

void VersionRun::parseCmdLine() {
    // Do nothing
}

void VersionRun::perform() {
    Interface::instance().catchVersion();
    Interface::instance().showLog(true);
}
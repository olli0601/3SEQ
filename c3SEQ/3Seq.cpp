/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   3Seq.cpp
 * CREATED ON:  02 June 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: The main source file for executable 3SEQ.
 *              This file contains main() function.
 *
 * HISTORY:     Version     Date            Description
 *              1.1104      2011-04-07      Original source from 3SEQ version 1.110407
 *              1.1106      2011-06-02      Reorganise the whole program (follow OOP design).
 *
 * VERSION:     1.1106
 * LAST EDIT:   02 June 2011
 */



//#define TESTING /* Use this flag whenever you want to run unit tests */
//#ifdef TESTING
//
//#include "gtest/gtest.h"
//#include "tests/testlist.h"
//
//int main(int argc, char** argv) {
//    testing::InitGoogleTest(&argc, argv);
//    return RUN_ALL_TESTS();
//}
//
//#else


#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <new>

#include "config.h"
#include "Interface.h"
#include "Run.h"

int main(int argc, char** argv) {
    Run* run = NULL;

    try {
        /* Use the name of executed file for the program name */
        if (argc >= 1) {
            Interface::instance().setProgramName(argv[0]);
        }

        /* Detect run-mode and return a run-object */
        run = Run::getRun(argc, argv);

        if (run == NULL) {
            /* Cannot detect run-mode, just show the program usage */
            Interface::instance().startProgram(DEFAULT_PROGRAM_DESC);
            Interface::instance().showPTableLink();
            Interface::instance().catchUsageModes();
            Interface::instance().showLog(true);

        } else {
            /* Here is the main flow of every analysis (of all kinds) */
            run->parseCmdLine();
            run->initLogFile();
            run->perform();
        }

        /* Throw a SUCCESS signal to indicate that the analysis has been done successfully.
         * Besides, this method is also necessary to save the log file. */
        Interface::instance().throwExitSignal(0);


    } catch (Interface::ExitSignal) {
        /* Do nothing */


    } catch (std::bad_alloc) {
        /* Reaching here means that no ExitSignal has been thrown and
         * the program is terminated due to out-of-memory. */
        Interface::instance() << "Out of memory.\n";
        try {
            Interface::instance().showError(true, true);
        } catch (...) {
            // Do nothing
        }


    } catch (...) {
        /* Reaching here means that no ExitSignal has been thrown and
         * the program is terminated by an unknown error. */
        Interface::instance()
                << "An unexplained exception occurs during the analysis.\n"
                << "For more help, email author: " << AUTHOR_EMAIL << "\n";
        try {
            Interface::instance().showError(true, true);
        } catch (...) {
            // Do nothing
        }
    }

    if (run != NULL) {
        delete run;
        run = NULL;
    }

    return 0;
}

//#endif /* TESTING */
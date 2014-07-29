/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   TestPRun.cpp
 * CREATED ON:  21 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: This run-type is used to test P-values. After the user enters
 *              values for M, N, K, the program will return the corresponding
 *              P-value.
 *
 * NOTES:       This run-type should only work in CONSOLE mode (@see Interface.h).
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-21  Created
 *                          2012-04-18  Limit number of digits of M, N, K to
 *                                      avoid integer overflow.
 *
 * VERSION:     1.0
 * LAST EDIT:   18 April 2012
 */

#include "CheckPRun.h"
#include "PTable.h"
#include "Util.h"
#include <cstdlib>

CheckPRun::CheckPRun(int argc, char** argv) : Run(argc, argv) {
    pTableFile = NULL;

    Interface::instance().startProgram("Test P-Values");
}

CheckPRun::CheckPRun(const CheckPRun& orig) {
    assert(false); // should never reach here
}

CheckPRun& CheckPRun::operator=(const CheckPRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

CheckPRun::~CheckPRun() {
}

void CheckPRun::parseCmdLine() {

    if (getRunArgsNum() >= 1) {
        string pTableFilePath = argVector[2];
        pTableFile = new PTableFile(pTableFilePath);

        if (!pTableFile->exists()) {
            Interface::instance() << "The file \"" << pTableFilePath << "\" does not exist.\n";
            Interface::instance().showError(true, true);
        }
    }
}

void CheckPRun::perform() {
    static const string PROMPT = Interface::DEFAULT_INDENT + "Enter M N K to test :  ";
    static const string VAL_SEPARATOR = "  ";
    static const int PARAMS_NUM = 3;

    vector<MNK> history;

    loadPTable(pTableFile);

    cout << Interface::DEFAULT_INDENT + "Press ESC when you want to exit.\n" << endl;

    while (true) {
        string lastBuf[PARAMS_NUM] = {"", "", ""};
        string currentBuf[PARAMS_NUM] = {"", "", ""};
        int bufCursor = 0;
        unsigned long historyCursor = history.size();

        cout << PROMPT;

        while (true) {
            int keyCode = Util::getKey();


            if (keyCode == 27) {
                /*ESC*/
                cout << endl << endl;
                return; // end run


            } else if (keyCode == 13 || keyCode == 10
                    || keyCode == 9 || keyCode == 32) {

                /* CR || NL || TAB || SPACE */
                if (currentBuf[bufCursor].length() > 0) {
                    bufCursor++;

                    if (bufCursor >= PARAMS_NUM) {
                        cout << endl;

                        MNK newMnk;
                        newMnk.m = atoi(currentBuf[0].c_str());
                        newMnk.n = atoi(currentBuf[1].c_str());
                        newMnk.k = atoi(currentBuf[2].c_str());
                        history.push_back(newMnk);

                        break; // finish inner while loop

                    } else {
                        cout << VAL_SEPARATOR;
                    }
                }


            } else if (keyCode >= '0' && keyCode <= '9') {
                if (currentBuf[bufCursor].length() < MAX_NUM_LENGTH) {
                    char key = static_cast<char> (keyCode);
                    currentBuf[bufCursor] += key;
                    cout << key;
                }


            } else {
                string tmpStr;

                /* Clear the line */
                tmpStr = "\r" + PROMPT;
                for (int i = 0; i <= bufCursor; i++) {
                    for (unsigned long k = 0; k < currentBuf[i].length(); k++) {
                        tmpStr += " ";
                    }
                    tmpStr += VAL_SEPARATOR;
                }
                cout << tmpStr;


                if (history.size() > 0
                        && (keyCode == -65 || keyCode == -66)) {

                    if (keyCode == -65) {
                        /* UP */
                        if (historyCursor <= 0) {
                            historyCursor = history.size();
                        } else {
                            historyCursor--;
                        }

                    } else if (keyCode == -66) {
                        /* DOWN */
                        if (historyCursor == history.size()) {
                            historyCursor = 0;
                        } else {
                            historyCursor++;
                        }
                    }

                    if (historyCursor == history.size()) {
                        bufCursor = 0;
                        for (int i = 0; i < PARAMS_NUM; i++) {
                            currentBuf[i] = lastBuf[i];
                            if (lastBuf[i].length() > 0) {
                                bufCursor = i;
                            }
                        }

                    } else {
                        bufCursor = 2;
                        MNK currentMNK = history[historyCursor];

                        stringstream ssM;
                        ssM << currentMNK.m;
                        currentBuf[0] = ssM.str();

                        stringstream ssN;
                        ssN << currentMNK.n;
                        currentBuf[1] = ssN.str();

                        stringstream ssK;
                        ssK << currentMNK.k;
                        currentBuf[2] = ssK.str();
                    }


                } else if (keyCode == 8 || keyCode == 127) {
                    /* BACKSPACE */
                    if (currentBuf[bufCursor].length() > 0) {
                        currentBuf[bufCursor] = currentBuf[bufCursor].substr(
                                0, currentBuf[bufCursor].length() - 1);

                    } else if (bufCursor > 0) {
                        bufCursor--;
                    }
                }


                /* Show the new line */
                tmpStr = "\r" + PROMPT;
                for (int i = 0; i < bufCursor; i++) {
                    tmpStr += currentBuf[i] + VAL_SEPARATOR;
                }
                tmpStr += currentBuf[bufCursor];
                cout << tmpStr;
            }


            /* If the last buffer and the current buffer are the same */
            if (historyCursor == history.size()) {
                for (int i = 0; i < PARAMS_NUM; i++) {
                    lastBuf[i] = currentBuf[i];
                }
            }
        }


        unsigned long historySize = history.size();
        if (historySize > 0) {
            MNK lastMNK = history[historySize - 1];
            int m = lastMNK.m;
            int n = lastMNK.n;
            int k = lastMNK.k;

            if (PTable::instance().canCalculateExact(m, n, k)) {
                Interface::instance()
                        << "P-value [" << m << ", " << n << ", " << k << "]  =  "
                        << PTable::instance().getExactPValue(m, n, k) << endl;
            } else {
                Interface::instance()
                        << "P-value [" << m << ", " << n << ", " << k << "]  ~  "
                        << PTable::instance().approxPValue(m, n, k)
                        << "  (approximated, H-S discrete)" << endl;
            }

            Interface::instance().showOutput(true);
        }
    }
}

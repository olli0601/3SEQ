/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   PTable.cpp
 * CREATED ON:  11 July 2011, 14:30
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: A singleton class that handles all the operations related to the
 *              P-value table (generate/save/load/search).
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-07-11   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   11 July 2011
 */

#include "PTable.h"
#include "PTableFile.h"


////////////////////////////////////////////////////////////////////////////////
//  STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

const double PTable::NaPVAL = -1.0;

/* 1 MB = 1048576 Byte - use double to be easy in calculating. */
const double PTable::BYTE_IN_MB = 1048576.0;
const double PTable::FLOAT_SIZE = static_cast<double> (sizeof (float));

////////////////////////////////////////////////////////////////////////////////

void PTable::clearTable(void) {
    if (indexArray) {
        delete[] indexArray;
        indexArray = NULL;
    }

    if (table) {
        delete[] table;
        table = NULL;
    }

    if (ykTableForLastK) {
        delete ykTableForLastK;
        ykTableForLastK = NULL;
    }
}

const int PTable::estimateMemNeededInMB(int newMSize, int newNSize, int newKSize) {
    unsigned long numOfStoredPVals = 0;
    /* Special cases will not be stored:
     *      m < 1
     *      n < 1      n < k
     *      k < 2      k < n - m + 1
     */
    for (int m = 1; m <= newMSize; m++) {
        for (int n = 1; n <= newNSize; n++) {
            int minK = (n - m + 1 > 2) ? n - m + 1 : 2;
            int maxK = (n < newKSize) ? n : newKSize;
            int numOfAcceptedK = (maxK >= minK) ? maxK - minK + 1 : 0;
            numOfStoredPVals += static_cast<unsigned long> (numOfAcceptedK);
        }
    }

    double memUsed = static_cast<double> (numOfStoredPVals) * FLOAT_SIZE / BYTE_IN_MB;

    int result = static_cast<int> (round(memUsed));
    if (result <= 0) result = 1;
    return result;
}

const unsigned long PTable::initIndexArray() {
    assert(indexArray == NULL);

    indexArray = new unsigned long[mSize * nSize];

    unsigned long totalNumOfStoredPVals = 0;

    /* Special cases will not be stored:
     *      m < 1
     *      n < 1      n < k
     *      k < 2      k < n - m + 1
     */
    for (int m = 1; m <= mSize; m++) {
        for (int n = 1; n <= nSize; n++) {
            indexArray[(m - 1) * nSize + n - 1] = totalNumOfStoredPVals;

            int minK = (n - m + 1 > 2) ? n - m + 1 : 2;
            int maxK = (n < kSize) ? n : kSize;
            int numOfAcceptedK = (maxK >= minK) ? maxK - minK + 1 : 0;
            totalNumOfStoredPVals += static_cast<unsigned long> (numOfAcceptedK);
        }
    }

    return totalNumOfStoredPVals;
}

void PTable::initialize(int newMSize, int newNSize, int newKSize) {
    clearTable();

    mSize = newMSize;
    nSize = newNSize;
    kSize = newKSize;

    assert(mSize > 0 && nSize > 0 && kSize > 0);
    numStoredVals = initIndexArray();

    table = new (nothrow) float[numStoredVals];

    if (!table) {
        Interface::instance() << "Unable to allocate memory for the P-value table.\n";
        Interface::instance().showError(true, true);
    }
}

void PTable::generateTable(int newMSize, int newNSize, int newKSize) {
    int memNeeded = estimateMemNeededInMB(newMSize, newNSize, newKSize);

    if (maxMemInMB < 0) {
        Interface::instance()
                << "The program needs ~" << memNeeded * 2 << "MB of RAM to generate the P-value table.\n"
                << "The real size of the table is ~"
                << memNeeded << "MB.\n"
                << "Do you want to continue?";

        bool yes = Interface::instance().yesNoQuestion(Interface::YES);
        if (!yes) Interface::instance().throwExitSignal(0);

    } else if (memNeeded > maxMemInMB) {
        Interface::instance() << "Not enough memory to generate the P-value table.\n"
                << "At least " << memNeeded << "MB of RAM needed.\n";
        Interface::instance().showError(true, true);
    }


    initialize(newMSize, newNSize, newKSize);

    ykTableForLastK = new YkTable(mSize, nSize, kSize);

    Interface::instance().initCounter("Generating P-value table", 2, kSize);

    /* k and n will be started from 2 since 0 and 1 are special cases */
    for (int k = 2; k <= kSize; k++) {
        Interface::instance() << "  -   Elapsed time :  " << Interface::instance().getElapsedTime();
        Interface::instance().count(k);

        ykTableForLastK->generateTable(k - 1);

        for (int m = 1; m <= mSize; m++) {

            /* Since n < k and n >= m + k are special cases */
            int maxN = (m + k < nSize) ? m + k : nSize;
            for (int n = k; n <= maxN; n++) {
                calculatePVal(m, n, k);
            }
        }

    }
    Interface::instance().finishCounting();

    delete ykTableForLastK;
    ykTableForLastK = NULL;
}

void PTable::calculatePVal(const int m, const int n, const int k) {
    assert(m >= 0 && n >= 0 && k >= 0);
    if (k > n
            || m == 0 /* k <= n */
            || n == 0 /* k == 0 */
            || k == 0 || k == 1 /* n >= 1 */
            || n - m >= k) {
        /* These cases are not necessary to be calculated */
        return;
    }

    assert(m <= mSize && n <= nSize && k <= kSize && m >= 1 && n >= 1 && k >= 2);

    unsigned long index = get1DIndex(m, n, k);
    assert(index >= 0 && index < numStoredVals);

    float fM = static_cast<float> (m);
    float fN = static_cast<float> (n);

    float pVal =
            (fM * getExactPValue(m - 1, n, k)
            + fN * getExactPValue(m, n - 1, k)
            + fN * ykTableForLastK->getYValue(m, n - 1, k - 1, k - 1)
            ) / (fM + fN);
    table[index] = pVal;
}

const bool PTable::saveToFile(const string& fileName) const {
    PTableFile pTableFile(fileName);
    return pTableFile.save(*this);
}


////////////////////////////////////////////////////////////////////////////////
//
//       Yk-Table class (private class)
//
////////////////////////////////////////////////////////////////////////////////

PTable::YkTable::YkTable(int newMSize, int newNSize, int newJSize) {
    mSize = newMSize;
    nSize = newNSize;
    jSize = newJSize;

    currentK = -1;

    lastYkk = new float[mSize * nSize];
    numStoredVals = initIndexArray();
    table = new (nothrow) float[numStoredVals];

    if (!table) {
        Interface::instance()
                << "Unable to allocate memory for Y-value table.\n"
                << "Try creating a smaller table.\n";
        Interface::instance().showError(true, true);
    }

    /* Initialise to reduce later calculation */
    for (unsigned long i = 0; i < numStoredVals; i++) {
        table[i] = 0.0f;
    }
}

PTable::YkTable::~YkTable(void) {
    if (indexArray) {
        delete[] indexArray;
        indexArray = NULL;
    }

    if (lastYkk) {
        delete[] lastYkk;
        lastYkk = NULL;
    }

    if (table) {
        delete[] table;
        table = NULL;
    }
}

const unsigned long PTable::YkTable::initIndexArray() {
    indexArray = new unsigned long[mSize * nSize];

    unsigned long totalNumOfValues = 0;

    /* m, n will begin from 1 since the cases where m, n = 0 are not 
     * needed to be stored */
    for (int m = 1; m <= mSize; m++) {
        for (int n = 1; n <= nSize; n++) {
            indexArray[(m - 1) * nSize + n - 1] = totalNumOfValues;

            int minJ = (n - m > 0) ? n - m : 0;
            int maxJ = (n < jSize) ? n : jSize;
            int numOfAcceptedJValues = (maxJ >= minJ) ? maxJ - minJ + 1 : 0;
            totalNumOfValues += static_cast<unsigned long> (numOfAcceptedJValues);
        }
    }

    return totalNumOfValues;
}

void PTable::YkTable::generateTable(const int k) {
    assert(k >= 1 && k <= jSize);

    if (currentK < 0 || k == 1) {
        /* Make sure this table is generated from k = 1 */
        assert(k == 1);

        /* Initialise lastYkk */
        for (int m = 1; m <= mSize; m++) {
            for (int n = 1; n <= nSize; n++) {
                /* Since k-1 = 0 and m, n > 0, all Y[m, n, k-1, k-1] = 0 */
                lastYkk[(m - 1) * nSize + n - 1] = 0.0f;
            }
        }

    } else {
        /* Make sure that the new k == currentK + 1 because the layer of k 
         * cannot be generated without the layer of k-1 */
        assert(k == currentK + 1);

        /* Update lastYkk */
        for (int m = 1; m <= mSize; m++) {
            int minN = k - 1; // k >= 2
            int maxN = k - 1 + m;
            if (maxN > nSize) maxN = nSize;

            /* Since k > n and k+m < n are special cases */
            for (int n = minN; n <= maxN; n++) {
                lastYkk[(m - 1) * nSize + n - 1] = getYValue(m, n, k - 1, k - 1);
            }
        }
    }

    currentK = k;

    /* Generate new layer */
    for (int m = 1; m <= mSize; m++) {
        int maxN = currentK + m;
        if (maxN > nSize) maxN = nSize;

        /* Since k > n and k+m < n are special cases */
        for (int n = currentK; n <= maxN; n++) {
            int minJ = (n - m > 0) ? n - m : 0;
            int maxJ = (n < jSize) ? n : jSize;

            for (int j = minJ; j <= maxJ; j++) {
                unsigned long ykIndex = get1DIndex(m, n, j);
                assert(ykIndex >= 0 && ykIndex < numStoredVals);

                float fM = static_cast<float> (m);
                float fN = static_cast<float> (n);

                if (j == 0) {
                    table[ykIndex] = (fM / (fM + fN))
                            * (getYValue(m - 1, n, k, 1) + getYValue(m - 1, n, k, 0));

                } else if (j == currentK) {
                    float Y_m_nPre_kPre_kPre = 0.0f; // Y[m, n-1, k-1, k-1]
                    if (n == 1) {
                        if (currentK == 1) {
                            /* n-1 = k-1 = j-1 = 0 */
                            Y_m_nPre_kPre_kPre = 1.0f;
                        }
                    } else {
                        Y_m_nPre_kPre_kPre = lastYkk[(m - 1) * nSize + (n - 1) - 1];
                    }

                    table[ykIndex] = (fN / (fM + fN))
                            * (Y_m_nPre_kPre_kPre + getYValue(m, n - 1, k, j - 1));

                } else {
                    table[ykIndex] = (fM * getYValue(m - 1, n, k, j + 1)
                            + fN * getYValue(m, n - 1, k, j - 1)
                            ) / (fM + fN);
                }
                
                
                ////////////////////////////////////////////////////////////////
                // @DEBUG: Test if Y[m, n, k, k-1] == Y[m, n, k-1, k-1]
                ////////////////////////////////////////////////////////////////
//                if (j == currentK - 1) {
//                    Interface::instance()
//                            << "Y[" << m << ", " << n << ", " << currentK << ", " << j << "] = "
//                            << table[ykIndex]
//                            << "   :   "
//                            << "Y[" << m << ", " << n << ", " << (currentK - 1) << ", " << (currentK - 1) << "] = "
//                            << lastYkk[(m - 1) * nSize + n - 1];
//                    Interface::instance().showLog(true);
//                }
            }
        }
    }
}


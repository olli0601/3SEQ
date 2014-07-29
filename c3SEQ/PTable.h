/* 
 * File:   PTable.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 11 July 2011, 14:30
 */

#ifndef PTABLE_H
#define	PTABLE_H

#include <fstream>
#include <cassert>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>

#include "Interface.h"
#include "stats.h"


using namespace std;

/**
 * Implement a singleton table that hold all necessary P-values.
 * @note For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class PTable {
public:
    /** Not a P-value. */
    static const double NaPVAL;

    ~PTable(void) {
        clearTable();
    }

    static PTable& instance(void) {
        static PTable instance;
        return instance;
    }

    /**
     * Pre-calculate the amount of RAM (in MB) that is needed to load a table of
     * the given size.
     * @param newMSize
     * @param newNSize
     * @param newKSize
     * @return RAM allotment (in MB).
     */
    static const int estimateMemNeededInMB(int newMSize, int newNSize, int newKSize);

    /**
     * Get the approximated RAM allotment (in MB) that the current table occupies.
     * @return The approximated RAM allotment (in MB) that the current table occupies.
     */
    const int getMemUsedInMB(void) const {
        double memUsed = static_cast<double> (numStoredVals) * FLOAT_SIZE / BYTE_IN_MB;

        int result = static_cast<int> (round(memUsed));
        if (result <= 0) result = 1;
        return result;
    }

    int getMSize(void) const {
        return mSize;
    }

    const int getNSize(void) const {
        return nSize;
    }

    const int getKSize(void) const {
        return kSize;
    }

    const unsigned long getNumStoredVals(void) const {
        return numStoredVals;
    }

    float* getDataPtr(void) const {
        return table;
    }

    /**
     * Test to see if the exact P-value can be calculated.
     * @param m
     * @param n
     * @param k
     * @return  <b>true</b> if the exact P-value can be calculated;
     *          <b>false</b>, otherwise.
     */
    const bool canCalculateExact(const int m, const int n, const int k) const {
        if (m < 0 || n < 0 || k < 0) return false;

        /* If k > n, the P-value is always zero */
        if (k > n) {
            return true;
        }

        /* In these case the p-value will be 1 */
        if (m == 0 /* k <= n */
                || n == 0 /* k == 0 */
                || k == 0 || k == 1 /* n >= 1 */
                || n - m >= k) {
            return true;
        }

        return (m <= mSize && n <= nSize && k <= kSize);
    }

    /**
     * Get the exact P-value.
     * @param m
     * @param n
     * @param k
     * @return  The exact P-value.
     * @note    To be safe, only execute this method when <b>canCalculateExact()</b>
     *          returns <b>true</b>; otherwise, false assertion may occur.
     */
    const float getExactPValue(const int m, const int n, const int k) const {
        assert(m >= 0 && n >= 0 && k >= 0);

        /* If k > n, the P-value is always zero */
        if (k > n) {
            return 0.0f;
        }

        /* In these case the p-value will be 1 */
        if (m == 0 /* k <= n */
                || n == 0 /* k == 0 */
                || k == 0 || k == 1 /* n >= 1 */
                || n - m >= k) {
            return 1.0f;
        }

        assert(m <= mSize && n <= nSize && k <= kSize);

        /* Reaching here means that m, n >= 1 and k >= 2 */
        unsigned long index = get1DIndex(m, n, k);
        assert(index >= 0 && index < numStoredVals);

        return table[index];
    }

    /**
     * Approximate P-value using Hogan Siegmund's methods.
     * @param m
     * @param n
     * @param k
     * @return The approximated P-value.
     *         The approximation methods are implemented with "long double" type
     *         for more precision.
     */
    const long double approxPValue(const int m, const int n, const int k) {
        /* Approximate using both continuous methods and discrete methods,
         * then take the greater value.
         * All the approximation methods are implemented with "long double" type
         * for more precision. */
        //long double contApprox = stats::siegmund::continuousApprox(m, n, k);
        //long double discreteApprox = stats::siegmund::discreteApprox(m, n, k);
        //return (contApprox > discreteApprox) ? contApprox : discreteApprox;

        /* Using only discrete approximation. */
        return stats::siegmund::discreteApprox(m, n, k);
    }

    /**
     * Set new size and allocate memory for the table.
     * @param newMSize
     * @param newNSize
     * @param newKSize
     */
    void initialize(int newMSize, int newNSize, int newKSize);

    /**
     * Allocate memory for the table and generate all P-values.
     * @param newMSize
     * @param newNSize
     * @param newKSize
     * @note    This P-value table is just temporarily stored on RAM.
     *          To save it into a file, call <b>saveToFile()</b> method.
     */
    void generateTable(int newMSize, int newNSize, int newKSize);

    /**
     * Save the current P-value table into file.
     * @param fileName
     * @return <b>true</b> if saving succeed, <b>false</b> otherwise.
     */
    const bool saveToFile(const string& fileName) const;



private:

    static const double BYTE_IN_MB;
    static const double FLOAT_SIZE;

    /* Disable constructor and assignment for singleton */
    PTable(void) {
        maxMemInMB = -1; // -1 means unlimited

        mSize = 0;
        nSize = 0;
        kSize = 0;
        numStoredVals = 0;

        indexArray = NULL;
        table = NULL;
        ykTableForLastK = NULL;
    }

    PTable(const PTable& orig);
    PTable operator=(const PTable& orig);

    /**
     * Delete all P-values in the table and release memory.
     */
    void clearTable(void);

    /**
     * Initialise the indexArray and return the space (maximum # of P-values) 
     * needed to be stored in this table.
     * @return The maximum # of P-values that are needed to be stored in this table.
     */
    const unsigned long initIndexArray();

    /**
     * Translate 3D-index (m, n, k) into 1D-index.<br>
     * This index is used to access the P-value table.
     * @param m
     * @param n
     * @param k
     * @return 1D-index
     */
    const unsigned long get1DIndex(int m, int n, int k) const {
        assert(m <= mSize && n <= nSize && k <= kSize);

        /* This function will only be called when m, n >= 1 and k >= 2 */
        assert(m >= 1 && n >= 1 && k >= 2);

        unsigned long mnIndex = indexArray[(m - 1) * nSize + n - 1];

        int minK = (n - m + 1 > 2) ? n - m + 1 : 2;
        int maxK = (n < kSize) ? n : kSize;
        assert(k >= minK && k <= maxK);
        unsigned long kIndex = static_cast<unsigned long> (k - minK);

        return mnIndex + kIndex;
    }

    /**
     * Calculate a P-value.
     * @param m
     * @param n
     * @param k
     */
    void calculatePVal(const int m, const int n, const int k);

    /**
     * Maximum memory allotment for storing this table
     */
    int maxMemInMB;

    /**
     * Total number of P-Value stored in this table.
     */
    unsigned long numStoredVals;

    int mSize;
    int nSize;
    int kSize;

    /**
     * The 2D array which will be used to find the 1D-index of P[m, n, k] 
     * in the P-value table.
     */
    unsigned long* indexArray;

    // &lt;     ~   <
    // &nbsp;   ~   SPACE
    /**
     * The 3D table that contains P-values P[m, n, k].                      <br>
     * Special cases will not be stored:                                    <br>
     * &nbsp;&nbsp;&nbsp;&nbsp;     m &lt; 1                                <br>
     * &nbsp;&nbsp;&nbsp;&nbsp;     n &lt; 1  &nbsp;&nbsp;  n &lt; k        <br>
     * &nbsp;&nbsp;&nbsp;&nbsp;     k &lt; 2  &nbsp;&nbsp;  k &lt;= n - m   <br>
     *                                                                      <br>
     * The definition of P-value table:                                     <br>
     * P[m, n, k] = sum[i:k->n] ( X[m, n, i] )                              <br>
     *                                                                      <br>
     * But P-value table will be calculated as the following formula:       <br>
     *                       P[m, n, k] = m * P[m-1, n, k] / (m + n) +      <br>
     * &nbsp;&nbsp;&nbsp;&nbsp;         + n * P[m, n-1, k] / (m + n) +      <br>
     * &nbsp;&nbsp;&nbsp;&nbsp;         + n * Y[m, n-1, k-1, k-1] / (m + n)
     */
    float* table;


    class YkTable;


    /**
     * The 3D table contains Y[m, n, k-1, j] where k is a fixed value. 
     * This table is just temporarily used to create the p-value table.
     */
    YkTable* ykTableForLastK;
};




////////////////////////////////////////////////////////////////////////////
//
//       Yk-Table class (inner class)
//
////////////////////////////////////////////////////////////////////////////

/**
 * In this Yk-table, Yk[m, n, j] = Y[m, n, k, j]  (k is a fixed number).
 * This means Yk-table is a layer (defined by k) of Y-table. <br>
 * This table is only used temporarily to calculate the P-value table.
 * 
 * @note For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class PTable::YkTable {
public:
    explicit YkTable(int newMSize, int newNSize, int newJSize);
    ~YkTable(void);

    void generateTable(const int k);

    /**
     * Get Y[m, n, k, j].
     * @param m
     * @param n
     * @param k This is given just to make sure that the currentK of this 
     *          table is the same as the given k.
     * @param j
     * @return  Y[m, n, k, j]
     */
    const float getYValue(int m, int n, int k, int j) const {
        assert(k == currentK);
        assert(m >= 0 && n >= 0 && j >= 0);

        if (n == 0) {
            return (k == 0 && j == 0) ? 1.0f : 0.0f;

        } else if (m == 0) {
            return (n == k && n == j) ? 1.0f : 0.0f;

        } else if (k == 0 && j == 0 /*n > 0*/) {
            return 0.0f;

        } else if (j > k
                || k > n || k < n - m
                || j > n || j < n - m) {
            return 0.0f;
        }

        assert(m <= mSize && n <= nSize && j <= jSize);

        unsigned long index = get1DIndex(m, n, j);
        assert(index >= 0 && index < numStoredVals);

        return table[index];
    };


private:

    YkTable(void) {
    };

    /* Disable assignment and copy constructor */
    YkTable(const YkTable& orig);
    YkTable operator=(const YkTable& orig);

    /**
     * Initialise the index array and return the space (maximum # of entities) 
     * needed for the Yk-table. <br>
     * This is used only in the constructor.
     * @return The space (maximum # of entities) needed for the Yk-table.
     */
    const unsigned long initIndexArray();

    /**
     * Translate 3D-index into 1D-index.<br>
     * This index is used to access the Yk-table.
     * @param m
     * @param n
     * @param j
     * @return 1D-index
     */
    const unsigned long get1DIndex(int m, int n, int j) const {
        assert(m >= 1 && m <= mSize && n >= 1 && n <= nSize && j >= 0 && j <= jSize);
        unsigned long mnIndex = indexArray[(m - 1) * nSize + n - 1];

        int minJ = (n - m > 0) ? n - m : 0;
        int maxJ = (n < jSize) ? n : jSize;
        assert(j >= minJ && j <= maxJ);
        unsigned long jIndex = static_cast<unsigned long> (j - minJ);

        return (mnIndex + jIndex);
    };


    int mSize;
    int nSize;
    int jSize;

    int currentK;

    /**
     * The maximum # of entities that will be stored in this Yk-table.
     */
    unsigned long numStoredVals;

    /**
     * The 3D table contains Y[m, n, k, j] where k is a fixed value. 
     */
    float* table;

    /**
     * The 2D table that contains the value of Yk'[m, n, k'] where k' = k-1 <br>
     * This table is needed to calculate Yk[m, n, k] (j == k) since:
     * Yk[m, n, k] = (n/(m+n)) * (Yk'[m, n-1, k'] + Yk[m, n-1, k-1])
     */
    float* lastYkk;

    /**
     * The 2D array which will be used to find the real index of Y[m, n, k, j]
     * in the yTableForLastK and P[m, n, k] in the p-value table.
     */
    unsigned long* indexArray;
};
////////////////////////////////////////////////////////////////////////////
//       Finish Yk-Table class
////////////////////////////////////////////////////////////////////////////




#endif	/* PTABLE_H */


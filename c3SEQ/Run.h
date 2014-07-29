/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   Run.h
 * CREATED ON:  12 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-12  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   12 August 2011
 */

#ifndef RUN_H
#define	RUN_H

#include <cassert>
#include <string>
#include <limits>
#include <sys/time.h>

//#include "tclap/CmdLine.h"
//#include "tclap/ValueArg.h"
//#include "tclap/SwitchArg.h"
//#include "tclap/ArgException.h"

#include "Interface.h"
#include "GenomeSet.h"
#include "String.h"
#include "SequenceFile.h"
#include "TextFile.h"
#include "PTableFile.h"


using namespace std;
//using namespace TCLAP;

class Run {
public:

    enum Mode {
        VERSION, HELP, GEN_P_TABLE, TEST_P_VALUES, INFO, SINGLE_BREAK_POINT,
        TRIPLET, FULL, WRITE, MATCH
    };

    virtual ~Run();

    static Run* getRun(int argc, char** argv);

    virtual void parseCmdLine(void);

    virtual void perform(void);
    
    virtual void initLogFile(void) const;
    
    virtual const bool isLogFileSupported(void) const = 0;

    virtual const Mode getMode(void) const = 0;


protected:
    /** Not an unsigned long number */
    static const unsigned long NaUL;

    /** Max value of unsigned long numbers */
    static const unsigned long MAX_UL;


    /**
     * This type of constructor should never be used. The reason for this 
     * to exist is just to make the inheritance less complex.
     */
    explicit Run();

    /**
     * This is the correct type of constructor.
     * @param argc
     * @param argv
     */
    explicit Run(int argc, char** argv);

    virtual const double getNumTotalTriplets(void) {
        if (numTotalTriplets < 0.0) {
            calculateTripletNum();
        }
        return numTotalTriplets;
    };

    virtual const double getNumTripletsForCorrection(void) {
        if (numTripletsForCorrection < 0.0) {
            calculateTripletNum();
        }
        return numTripletsForCorrection;
    };

    /**
     * Calculate the total number of triplets in this run and the number of
     * triplets that will be used for the statistical correction.
     */
    virtual void calculateTripletNum(void);

    /**
     * Delete some argument from <b>argVector</b>. The argument <b>1</b> to 
     * argument <b>numDelete</b> will be deleted (argument <b>0</b> will be kept).
     * @param numDelete
     */
    virtual void deleteCmdLineArgs(int numDelete);

    virtual const int getRunArgsNum(void) const;

    /**
     * Remove neighbour and identical genomes if necessary.
     */
    virtual void removeNeighborAndIdenticalGenomes(void);

    /**
     * Add a P-value into the histogram. This method should be use instead of
     * adding directly to avoid array-out-of-bound.
     * @param pValue
     */
    virtual void addPValIntoHistogram(double pValue);

    /**
     * Save the histogram of P-values into file.
     */
    virtual void savePValHistogram(void);
    
    /**
     * If the <b>pTableFile</b> is NOT <b>NULL</b> try to load that file;
     * otherwise, search for a P-value table file in current directory
     * (and load it into RAM if found).
     * @param pTableFile    A nullable pointer to a P-table file.
     * @note    If no P-value table is successfully loaded, the program will be
     *          terminated by this method.
     */
    virtual void loadPTable(PTableFile* pTableFile);

    /** Indicates if children & parents sequences are in different files */
    bool isChildParentInDiffFiles;

    GenomeSet* parentDataset;
    GenomeSet* childDataset;

    TextFile* fileSubset;

    ////////////////////////////////////////////////////////////////////////////
    // Command Line Options
    ////////////////////////////////////////////////////////////////////////////

    int argCount;

    char** argVector;

    /** 
     * The ID for this run. This ID will be added after the name of the output file.<b>
     * This is useful when running multiple analysis concurrently, so, the names of 
     * output files will not be duplicate.
     */
    string id;

    /**
     * Indicates if skipped triples will be recorded.
     * @default false
     */
    bool cloIsSkippedTriplesRecorded;

    /** 
     * Uses all sites rather than just polymorphic sites.<br>
     * Useful when using the -out option to create an output file of 
     * the some subsequences; and the -write mode as well. 
     * @default false
     */
    bool cloUseAllSites;

    /**
     * Indicates if algorithm should stop once a significant triple has been found.
     * @default false
     */
    bool cloStopAtFirstRecombinant;

    /**
     * Indicates if the program is called by a Python script.
     * @default false
     */
    bool cloIsStartedByPythonScript;

    /**
     * Rejection threshold
     * @default 0.05
     */
    double cloRejectThreshold;

    /** 
     * First nucleotide to be analyzed, e.g. <b>-f100 -l200</b> considers 
     * nucleotide positions 100-200 inclusive (first position is 1, not 0); 
     * default is all positions. 
     */
    unsigned long cloFirstNucleotidePos;

    /** 
     * Last nucleotide to be analyzed, e.g. <b>-f100 -l200</b> considers 
     * nucleotide positions 100-200 inclusive (first position is 1, not 0); 
     * default is all positions. 
     */
    unsigned long cloLastNucleotidePos;

    /**
     * This can be used in combination with the <b>-f</b> and <b>-l</b> switches 
     * which set the variables <b>firstNucleotidePos</b> and <b>lastNucleotidePos</b>.  
     * Normally, if you enter <b>-f20</b> and <b>-l40</b> then only the 
     * subsequence from nucleotides 20 to 40 is analysed. However, if you set 
     * <b>isFirstToLastCutOut = true</b>, then the sequence from nt 20-40
     * would be cut out and the remainder of the sequence would be analyzed.
     * @default false
     */
    bool cloIsFirstToLastCutOut;

    /**
     * This is the controlled by the command line argument <b>-#</b>.
     * When user selects <b>-#</b>, this is set to <b>false</b>, and multiple
     * comparisons corrections will use <b>uTestedTriplesNumInRun</b> instead 
     * of <b>uTestedTriplesNum</b> in the Dunn-Sidak correction. <br>
     * The purpose of this switch is to be able to query a sequence
     * e.g. <b>-b35 -e35</b> and then use the <b>-#</b> switch to correct for 
     * the comparisons done in that run only. <br>
     * When running in parallel, make sure <b>-#</b> switch is NOT used.
     * @default true
     */
    bool cloIsAllTriplesUsed;

    /** 
     * @default false
     */
    bool cloIsSubsetRemoved;

    /** 
     * @default false
     */
    bool cloIsIdenticalSeqRemoved;

    /** 
     * If the distance between 2 sequences is equal or smaller than this distance, 
     * they will be considered as neighbour sequences. Then, 1 of them will be 
     * remove from the dataset if the flag <b>isIdenticalSeqRemoved</b> is turned on.
     * @default 0
     */
    unsigned long cloRemoveDistance;

    /**
     * Indicate that all break point positions for all candidate recombinants
     * will be calculated; default is to calculate the best break point positions 
     * per candidate recombinant.
     * @default false
     * @note    command line switch <b>-bp-all</b>
     */
    bool cloAllBpCalculated;

    /**
     * Indicate that no break point positions will be calculated;
     * default is to calculate the best break point positions per 
     * candidate recombinant.
     * @default false
     * @note    command line switch <b>-bp-none</b>
     */
    bool cloNoBpCalculated;

    /**
     * You should be able to turn this on, i.e. suppress the
     * creation of the 3s.rec file, with the "-nr" option
     */
    bool cloNo3sRecFile;

    /** 
     * @default false
     */
    bool cloUseSiegmundApprox;

    /** 
     * The minimum length for a recombinant segment to be considered as a 
     * long recombinant.
     * @default 100
     */
    unsigned long cloMinLongRecombinantLength;

    /** 
     * @default Unknown
     */
    SequenceFile::Type cloOutputFileType;

    ////////////////////////////////////////////////////////////////////////////


private:
    Run(const Run& orig);

    Run& operator=(const Run& rhs);

    void generateAutoID(void);

    void verifyFirstAndLastNuPos(void);

    void verifyBeginAndEndSeqIndex(void);

    void applySubset(void);

    void extractGenomes(void);

    void removeNonPolymorphic(void);

    /* Histogram of p-values */
    unsigned long pValsHistogram[PVAL_HISTOGRAM_SIZE];

    /** Total number of triplets that will be tested by this run */
    double numTotalTriplets;

    /** Number of comparisons used in Dunn-Sidak corrections */
    double numTripletsForCorrection;

    string pValHistogramFileName;


    ////////////////////////////////////////////////////////////////////////////
    // Command Line Options
    ////////////////////////////////////////////////////////////////////////////

    /**
     * Indicates if the analysis just run on the 3rd nucleotide (in codon).
     * @default false
     */
    bool cloThirdPositionsOnly;

    /**
     * Indicates if the analysis just run on the 1st and 2nd nucleotide (in codon).
     * @default false
     */
    bool cloFirstAndSecondPositionsOnly;

    /**
     * This index describes what child sequence to start at 
     * when looping through child sequences
     * @default 0
     */
    unsigned long cloBeginSequence;

    /**
     * This index describes what child sequence to end at when looping through 
     * child sequences; if all child sequences are to be tested, 
     * you must set <b>endSequence = N</b>.
     * @default 0
     */
    unsigned long cloEndSequence;

    ////////////////////////////////////////////////////////////////////////////
};

#endif	/* RUN_H */


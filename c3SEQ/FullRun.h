/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   FullRun.h
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

#ifndef FULLRUN_H
#define	FULLRUN_H

#include <cassert>
#include <string>
#include <set>
#include "Run.h"
#include "SequenceFile.h"
#include "PTableFile.h"


class FullRun : public Run {
public:
    explicit FullRun(int argc, char** argv);

    virtual ~FullRun();

    virtual const bool isLogFileSupported(void) const {
        return true;
    };
    
    virtual const Mode getMode() const {
        return FULL;
    };
    
    virtual void parseCmdLine();
    
    virtual void perform();


private:
    FullRun(const FullRun& orig);

    FullRun& operator=(const FullRun& rhs);
    
    /**
     * Test the dataset, and show some info about the dataset.
     */
    void verifyData(void);
    
    /**
     * Prepare for the run.
     */
    void setup(void);
    
    /**
     * Run the analysis. This is the core of "full-run".
     */
    void progress(void);
    
    /**
     * Show the progress counter.
     * @param currentLoop
     * @param isFinish
     */
    void showProgress(double currentLoop, bool isFinish) const;
    
    /**
     * Show the result of the analysis.
     */
    void displayResult(void);
    
    void recordRecombinantTriplet(Triplet* triplet);    
    
    /** P-value table file */
    PTableFile* pTableFile;
    
    /** Minimum P-value that has been reached during this run */
    double minPVal;
    
    /**
     * Number of triplets where recombinations have been detected
     * (significant P-values were returned).
     */
    double numRecombinantTriplets;
    
    /** Number of triplets where P-values are computed exactly */
    double numComputedExactly;
       
    /** Number of triplets that where P-values were approximated */
    double numApproximated;
    
    /**
     * Number of triplets that have been skipped during the analysis
     * (due to P-value table out-of-bound & no use of approximation).
     */
    double numSkipped;

    /** The length of the longest recombinant segment */
    unsigned long longestRecombinantSegment;
    
    TextFile* fileSkippedTriplets;
    TextFile* fileRecombinants;
    TextFile* fileLongRecombinants;
    
    /** Number of random sequences in each random loop */
    unsigned long randomSeqNum;
    
    /** Number of random iteration */
    unsigned long randomLoopNum;
    
    /** Indicates the current run is in random mode. */
    bool isRandomMode;   
    
    /**
     * This variable is only used in random mode to indicate the number of
     * random iterations where there is at least 1 recombination event detected.
     */
    unsigned long recombinationSensitivity;
};

#endif	/* FULLRUN_H */


/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   MatchRun.h
 * CREATED ON:  25 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-25  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   25 August 2011
 */

#ifndef MATCHRUN_H
#define	MATCHRUN_H

#include <cassert>
#include <vector>
#include "Run.h"

class MatchRun : public Run {
public:
    explicit MatchRun(int argc, char** argv);

    virtual ~MatchRun();

    virtual const bool isLogFileSupported(void) const {
        return true;
    };

    virtual const Mode getMode() const {
        return MATCH;
    };

    virtual void parseCmdLine();

    virtual void perform();


private:
    MatchRun(const MatchRun& orig);

    MatchRun& operator=(const MatchRun& rhs);

    void displayOutput(void) const;
    
    void saveMatchIntoFile(void) const;
    
    string queryAccession;
    
    /** 
     * The minimum distance of sequence shown that match the query sequence. 
     * @default 0
     */
    unsigned long minMatchDistance;

    /**
     * The maximum distance of sequence shown that match the query sequence. 
     * @default 0
     */
    unsigned long maxMatchDistance;

    struct MatchResult {
        Genome* genome;
        unsigned long distance;
    };
    
    vector<MatchResult> results;
};

#endif	/* MATCHRUN_H */


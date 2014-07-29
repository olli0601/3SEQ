/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   TripletRun.h
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

#ifndef TRIPLETRUN_H
#define	TRIPLETRUN_H

#include <cassert>
#include <string>
#include "Run.h"

using namespace std;

class TripletRun : public Run {
public:
    explicit TripletRun(int argc, char** argv);

    virtual ~TripletRun();

    virtual const bool isLogFileSupported(void) const {
        return true;
    };

    virtual const Mode getMode() const {
        return TRIPLET;
    };
    
    virtual void parseCmdLine();
    
    virtual void perform();


private:
    /** Number of sequences in triplet run */
    static const int SEQ_NUM = 3;
    
    TripletRun(const TripletRun& orig);

    TripletRun& operator=(const TripletRun& rhs);

    /** 
     * List of accession numbers of the sequences in TripletRun mode.<br>
     * 0. First parent<br>
     * 1. Second parent<br>
     * 2. Child
     */
    string accessionList[SEQ_NUM];

};

#endif	/* TRIPLETRUN_H */


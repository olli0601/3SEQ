/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   SingleBpRun.h
 * CREATED ON:  31 August 2011, 16:14
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date            Description
 *              1.0         2011-08-31      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   31 August 2011
 */

#ifndef SINGLEBPRUN_H
#define	SINGLEBPRUN_H

#include <cassert>
#include "Run.h"

class SingleBpRun : public Run {
public:
    explicit SingleBpRun(int argc, char** argv);

    virtual ~SingleBpRun();

    virtual const bool isLogFileSupported(void) const {
        return true;
    };

    virtual const Mode getMode() const {
        return SINGLE_BREAK_POINT;
    };

    virtual void parseCmdLine();

    virtual void perform();


private:
    SingleBpRun(const SingleBpRun& orig);

    SingleBpRun& operator=(const SingleBpRun& rhs);

    void prepare(void);

    void progress(void);
    
    void displayResult(void);

    double minPVal;

    /** Number of triplets where P-values are computed */
    double numComputed;
};

#endif	/* SINGLEBPRUN_H */


/* 
 * File:   WriteRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 19 August 2011, 16:26
 */

#ifndef WRITERUN_H
#define	WRITERUN_H

#include <cassert>
#include <string>
#include <iostream>
#include "Run.h"
#include "SequenceFile.h"

class WriteRun : public Run {
public:
    explicit WriteRun(int argc, char** argv);

    virtual ~WriteRun();

    virtual const bool isLogFileSupported(void) const {
        return false;
    };

    virtual const Mode getMode() const {
        return WRITE;
    };

    virtual void parseCmdLine();

    virtual void perform();


private:
    WriteRun(const WriteRun& orig);

    WriteRun& operator=(const WriteRun& rhs);
    
    SequenceFile* outputFile;
};

#endif	/* WRITERUN_H */


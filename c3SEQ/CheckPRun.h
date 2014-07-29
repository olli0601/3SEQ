/* 
 * File:   TestPRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 21 August 2011, 17:39
 */

#ifndef TESTPRUN_H
#define	TESTPRUN_H

#include <stdint.h> // for detecting the OS architecture through __WORDSIZE
#include <cassert>
#include "Run.h"
#include "PTableFile.h"

class CheckPRun : public Run {
public:
    struct MNK {
        int m;
        int n;
        int k;
    };
    
    explicit CheckPRun(int argc, char** argv);

    virtual ~CheckPRun();
    
    virtual const bool isLogFileSupported(void) const {
        return false;
    };

    virtual const Mode getMode() const {
        return TEST_P_VALUES;
    };
    
    virtual void parseCmdLine();
    
    virtual void perform();


private:

#if __WORDSIZE == 64   // 64-bit
    /** The maximum number of digits that M, N, K can have. */
    static const unsigned int MAX_NUM_LENGTH = 9;

#else   // 32-bit
    /** The maximum number of digits that M, N, K can have. */
    static const unsigned int MAX_NUM_LENGTH = 4;
#endif
    
    CheckPRun(const CheckPRun& orig);

    CheckPRun& operator=(const CheckPRun& rhs);
    
    /** Path to P-value table file */
    PTableFile* pTableFile;
    
};

#endif	/* TESTPRUN_H */


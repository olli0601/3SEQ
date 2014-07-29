/* 
 * File:   GenPRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 19 August 2011, 17:54
 */

#ifndef GENPRUN_H
#define	GENPRUN_H

#include <cassert>
#include "Run.h"

class GenPRun : public Run {
public:
    explicit GenPRun(int argc, char** argv);

    virtual ~GenPRun();

    virtual const bool isLogFileSupported(void) const {
        return false;
    };

    virtual const Mode getMode() const {
        return GEN_P_TABLE;
    };
    
    virtual void parseCmdLine();
    
    virtual void perform();

    
private:
    GenPRun(const GenPRun& orig);

    GenPRun& operator=(const GenPRun& rhs);

    string pTableFilePath;
    
    int pTableSize;
};

#endif	/* GENPRUN_H */


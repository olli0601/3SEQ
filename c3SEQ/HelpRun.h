/* 
 * File:   HelpRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 12 August 2011, 17:08
 */

#ifndef HELPRUN_H
#define	HELPRUN_H

#include <cassert>
#include "Run.h"

//class Run;

class HelpRun : public Run {
public:

    explicit HelpRun(int argc, char** argv);

    virtual ~HelpRun();

    virtual const bool isLogFileSupported(void) const {
        return false;
    };

    virtual const Mode getMode() const {
        return HELP;
    };

    virtual void parseCmdLine();
    
    virtual void perform();


private:

    HelpRun(const HelpRun& orig);

    HelpRun& operator=(const HelpRun& rhs);

};

#endif	/* HELPRUN_H */


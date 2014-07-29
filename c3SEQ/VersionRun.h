/* 
 * File:   VersionRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 17 August 2011, 11:57
 */

#ifndef VERSIONRUN_H
#define	VERSIONRUN_H

#include <cassert>
#include "Run.h"
#include "Interface.h"

class VersionRun : public Run {
public:
    explicit VersionRun(int argc, char** argv);

    virtual ~VersionRun();

    virtual const bool isLogFileSupported(void) const {
        return false;
    };

    virtual const Mode getMode() const {
        return VERSION;
    };
    
    virtual void parseCmdLine();

    virtual void perform();


private:
    VersionRun(const VersionRun& orig);

    VersionRun& operator=(const VersionRun& rhs);


};

#endif	/* VERSIONRUN_H */


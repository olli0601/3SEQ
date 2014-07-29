/* 
 * File:   InfoRun.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 18 August 2011, 17:50
 */

#ifndef INFORUN_H
#define	INFORUN_H

#include <cassert>
#include "Run.h"
#include "SequenceFile.h"

class InfoRun : public Run {
public:
    explicit InfoRun(int argc, char** argv);

    virtual ~InfoRun();

    virtual const bool isLogFileSupported(void) const {
        return true;
    };

    virtual const Mode getMode() const {
        return INFO;
    };

    virtual void parseCmdLine();

    virtual void perform();


private:
    InfoRun(const InfoRun& orig);

    InfoRun& operator=(const InfoRun& rhs);

};

#endif	/* INFORUN_H */


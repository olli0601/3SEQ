/* 
 * File:   Segment.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 28 July 2011, 11:28
 */

#ifndef SEGMENT_H
#define	SEGMENT_H

#include "Sequence.h"

class Segment : public Sequence {
public:

    explicit Segment(const string& segmentStr, const string& newName)
    : Sequence(segmentStr, newName) {
    }

    ~Segment() {
        // Do nothing
    }


private:
    Segment(const Segment& orig);

    Segment& operator=(const Segment& rhs);

};

#endif	/* SEGMENT_H */


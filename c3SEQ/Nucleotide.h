/* 
 * File:   Nucleotide.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 16 June 2011, 17:38
 */

#ifndef NUCLEOTIDE_H
#define	NUCLEOTIDE_H

#include <ctype.h>
//#include "Sequence.h"

using namespace std;


class Sequence; // forward declaration to get rid of recursive include error.

class Nucleotide {
public:
    static const char GAP = '-';
    static const char ADENINE = 'A';
    static const char CYTOSINE = 'C';
    static const char GUANINE = 'G';
    static const char THYMINE = 'T';

    explicit Nucleotide(const char& newCode, unsigned long newPosition,
            Sequence * newContainingSequence)
    : position(newPosition), containingSequence(newContainingSequence) {

        code = toupper(newCode);
        replaceUnclearByGap();
    }

    ~Nucleotide() {
        // Do nothing
    }

    /**
     * Indicates if this is a gap.
     * @return  <b>true</b> if this is a gap; otherwise, <b>false</b>.
     */
    bool const isGap(void) const {
        return (code == GAP);
    }

    const char getCode() const {
        return code;
    }

    /**
     * Get the sequence that contain this nucleotide.
     * @return The sequence that contain this nucleotide.
     */
    Sequence * const getContainingSequence() const {
        return containingSequence;
    }

    const unsigned long getBasePairPos() const {
        return position;
    }


private:
    /* Disable copy constructor and assignment */
    Nucleotide(const Nucleotide& orig);
    Nucleotide& operator=(const Nucleotide& rhs);

    /**
     * Indicates if this nucleotide is coded by A, T, G, C or not.
     * @return  <b>true</b> if this is coded by any one of A, T, G, C; 
     *          otherwise, <b>false</b>
     */
    bool const isClear(void) const {
        for (int i = 0; i < CLEAR_CODE_NUM; i++) {
            if (CLEAR_CODE[i] == toupper(code)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Replace this nucleotide by gap if it is ambiguous or invalid.
     */
    void replaceUnclearByGap(void) {
        if (!isClear()) {
            code = GAP;
        }
    }

    //    /**
    //     * Indicates if the code of this nucleotide is valid (clear, ambiguous, gap) 
    //     * or not.
    //     * @return  <b>true</b> if the code of this nucleotide is valid; 
    //     *          otherwise, <b>false</b>.
    //     */
    //    bool const isValid(void) const {
    //        return (isClear() || isAmbiguous() || isGap());
    //    };
    //
    //    /**
    //     * Indicates if this is an ambiguous nucleotide.
    //     * @return  <b>true</b> if ambiguous; otherwise, <b>false</b>.
    //     */
    //    bool const isAmbiguous(void) const {
    //        for (int i = 0; i < AMBIGUOUS_CODE_NUM; i++) {
    //            if (AMBIGUOUS_CODE[i] == toupper(code)) {
    //                return true;
    //            }
    //        }
    //        return false;
    //    };

    char code;
    unsigned long position;
    Sequence * containingSequence;

    /* Predefine acceptable nucleotide codes */
    static const int CLEAR_CODE_NUM = 4;
    static const char CLEAR_CODE[CLEAR_CODE_NUM];
    //    static const int AMBIGUOUS_CODE_NUM = 11;
    //    static const char AMBIGUOUS_CODE[AMBIGUOUS_CODE_NUM];
};


#endif	/* NUCLEOTIDE_H */


/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   Nucleotide.cpp
 * CREATED ON:  16 June 2011, 17:38
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: 
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-06-16  Created   
 * 
 * VERSION:     1.0
 * LAST EDIT:   16 June 2011
 */

#include "Nucleotide.h"


/* Initialise static members */
const char Nucleotide::CLEAR_CODE[CLEAR_CODE_NUM] = {
    ADENINE, GUANINE, THYMINE, CYTOSINE
};
//const char Nucleotide::AMBIGUOUS_CODE[AMBIGUOUS_CODE_NUM] = {
//    'R' /* G or A */,
//    'Y' /* C or T */,
//    'M' /* A or C */,
//    'K' /* G or T */,
//    'S' /* G or C */,
//    'W' /* G or A */,
//    'H' /* A or C or T */,
//    'B' /* C or G or T */,
//    'V' /* A or C or G */,
//    'D' /* A or G or T */,
//    'N' /* any */
//};

/* Finish initialising static members */


//Nucleotide::Nucleotide(const char& newCode, unsigned long newPosition,
//        Sequence * newContainingSequence)
//: position(newPosition), containingSequence(newContainingSequence) {
//
//    code = toupper(newCode);
//    replaceUnclearByGap();
//}


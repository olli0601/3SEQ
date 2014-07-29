/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   Sequence.cpp
 * CREATED ON:  27 July 2011, 16:46
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: 
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-07-27  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   27 July 2011
 */

#include "Sequence.h"


//Sequence::Sequence(const string& dnaStr, const string& newName)
//: name(newName) {
//
//    try {
//        nucleotides.reserve(dnaStr.length());
//
//        for (unsigned long i = 0; i < dnaStr.length(); i++) {
//            Nucleotide * nucleotidePtr = new Nucleotide(dnaStr[i], i, this);
//            nucleotides.push_back(nucleotidePtr);
//        }
//
//    } catch (bad_alloc) {
//        Interface::instance() << "Not enough memory to create new nucleotide object.\n";
//        Interface::instance().showError(true, true);
//    } catch (length_error) {
//        Interface::instance() << "Not enough memory to store the sequence " << name << "\n";
//        Interface::instance().showError(true, true);
//    }
//}

//Sequence::Sequence(const string& newName)
//: name(newName) {
//}

//Sequence::~Sequence() {
//    for (unsigned long i = 0; i < nucleotides.size(); i++) {
//        if (nucleotides[i]) {
//            delete nucleotides[i];
//        }
//    }
//}

const string Sequence::toString(void) const {
    string sequenceStr = "";

    for (unsigned long i = 0; i < nucleotides.size(); i++) {
        sequenceStr += nucleotides[i]->getCode();
    }

    return sequenceStr;
}

const unsigned long Sequence::ntDistanceTo(Sequence* another, bool ignoreGap) const {
    assert(another && this->getLength() == another->getLength());

    unsigned long distance = 0;

    for (unsigned long i = 0; i < this->getLength(); i++) {
        if (!ignoreGap
                || (!this->getNuAt(i)->isGap() && !another->getNuAt(i)->isGap())) {

            if (this->getNuAt(i)->getCode() != another->getNuAt(i)->getCode()) {
                distance++;
            }
        }
    }
    return distance;
}

const bool Sequence::isIdenticalTo(Sequence* another, bool ignoreGap) const {
    return (this->ntDistanceTo(another, ignoreGap) == 0);
}
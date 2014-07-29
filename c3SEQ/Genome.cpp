/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   Genome.cpp
 * CREATED ON:  27 July 2011, 16:33
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: This is a kind of sequence that can contain multiple segments.
 *              The name of this class is for distinguishing it among other 
 *              types of sequence.
 *              This is not a real genome since it may not contain a complete 
 *              set of genes.
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-07-27  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   27 July 2011
 */

#include "Genome.h"
#include "Triplet.h"
#include "Segment.h"

//Genome::Genome(const string& accession) : Sequence(accession) {
//    bestRecombinantTriplet = NULL;
//    recombinantType = Na_REC;
//}

Genome::~Genome() {
    for (unsigned long i = 0; i < segments.size(); i++) {
        if (segments[i]) {
            delete segments[i];
        }
    }

    if (bestRecombinantTriplet) delete bestRecombinantTriplet;

    /* All nucleotides in this genome have been deleted through deleting its 
     * segments. So, just clear the nucleotide vector to make sure that they 
     * will not be delete again when the destructor of the super class is 
     * executed. */
    nucleotides.clear();
}

Genome::Cursor * const Genome::newCursor(void) const {
    return new Cursor(this);
}

void Genome::extractByTemplate(const ExtractTemplate& extractTemplate) {
    unsigned long templateSize = extractTemplate.getSize();

    nucleotides.clear();
    nucleotides.reserve(templateSize);

    for (unsigned long i = 0; i < templateSize; i++) {
        nucleotides.push_back(
                segments[extractTemplate.segmentIndexAt(i)]->getNuAt(extractTemplate.nuIndexOnSegmentAt(i))
                );
    }

    assert(nucleotides.size() == templateSize);
}

unsigned long Genome::getBasePairPosOf(Nucleotide* nucleotide) const {
    Sequence* seq = nucleotide->getContainingSequence();
    unsigned long segIndex = getSegmentIndex(seq);
    
    unsigned long result = 0;
    
    for (unsigned long i = 0; i < segIndex; i++) {
        result += segments[i]->getLength();
    }
    
    result += nucleotide->getBasePairPos();
    
    return result;
}

double Genome::bpDistanceBetween(Nucleotide* firstNu, Nucleotide* secondNu) const {
    double distance = 0.0;
    bool isNegativeDistance = false;

    Sequence* firstSeg = firstNu->getContainingSequence();
    Sequence* secondSeg = secondNu->getContainingSequence();

    unsigned long firstBpPos = firstNu->getBasePairPos();
    unsigned long secondBpPos = secondNu->getBasePairPos();

    if (firstSeg == secondSeg) {
        /* 2 nucleotide are on the same segment */
        getSegmentIndex(firstSeg); // just to make sure the segment is in this genome

        unsigned long min = firstBpPos;
        unsigned long max = secondBpPos;

        if (firstBpPos > secondBpPos) {
            min = secondBpPos;
            max = firstBpPos;
            isNegativeDistance = true;
        }

        distance = static_cast<double> (max - min);


    } else {
        /* 2 nucleotide are on different segments */
        unsigned long firstSegIndex = getSegmentIndex(firstSeg);
        unsigned long secondSegIndex = getSegmentIndex(secondSeg);

        if (firstSegIndex > secondSegIndex) {
            unsigned long tmpSegIndex = firstSegIndex;
            firstSegIndex = secondSegIndex;
            secondSegIndex = tmpSegIndex;

            unsigned long tmpBpPos = firstBpPos;
            firstBpPos = secondBpPos;
            secondBpPos = tmpBpPos;

            isNegativeDistance = true;
        }

        distance = static_cast<double> (segments[firstSegIndex]->getLength() - firstBpPos);
        for (unsigned long i = firstSegIndex + 1; i < secondSegIndex; i++) {
            distance += static_cast<double> (segments[i]->getLength());
        }
        distance += static_cast<double> (secondBpPos);
    }

    if (isNegativeDistance) {
        distance = -distance;
    }

    return distance;
}

void Genome::setBestRecombinantTriplet(Triplet* triplet) {
    assert(triplet == NULL || this == triplet->getChild());
    if (bestRecombinantTriplet && bestRecombinantTriplet != triplet) {
        delete bestRecombinantTriplet;
    }
    bestRecombinantTriplet = triplet;
}




////////////////////////////////////////////////////////////////////////////////
//      CURSOR CLASS
////////////////////////////////////////////////////////////////////////////////

Genome::Cursor::Cursor(Genome const * const newGenomePtr) : genomePtr(newGenomePtr) {
    putAtFirstNu();
};

void Genome::Cursor::setPosition(unsigned long segIndex, unsigned long nuPosInSeg) {
    curSegPtr = genomePtr->getSegment(segIndex);
    curNuPtr = curSegPtr->getNuAt(nuPosInSeg);
};

void Genome::Cursor::setPosition(Nucleotide * const nuPtr) {    
    curNuPtr = nuPtr;  
    
    if (curNuPtr == NULL)
        curSegPtr = NULL;
    else
        curSegPtr = curNuPtr->getContainingSequence();
}

void Genome::Cursor::putAtFirstNu(void) {
    curSegPtr = NULL;
    curNuPtr = NULL;

    if (genomePtr->getSegmentNum() > 0) {
        curSegPtr = genomePtr->getSegment(0);
        if (curSegPtr->getLength() > 0) {
            curNuPtr = curSegPtr->getNuAt(0);
        }
    }
}

void Genome::Cursor::putAtLastNu(void) {
    curSegPtr = NULL;
    curNuPtr = NULL;

    unsigned long segNum = genomePtr->getSegmentNum();
    if (segNum > 0) {
        curSegPtr = genomePtr->getSegment(segNum - 1);

        unsigned long segLength = curSegPtr->getLength();
        if (segLength > 0) {
            curNuPtr = curSegPtr->getNuAt(segLength - 1);
        }
    }
}

Nucleotide * const Genome::Cursor::current(void) const {
    return curNuPtr;
}

Nucleotide * const Genome::Cursor::next(void) {
    if (curNuPtr == NULL) return NULL;
    
    unsigned long lastNuPos = curNuPtr->getBasePairPos();

    if (lastNuPos < curSegPtr->getLength() - 1) {
        curNuPtr = curSegPtr->getNuAt(lastNuPos + 1);

    } else {
        unsigned long lastSegIndex = genomePtr->getSegmentIndex(curSegPtr);
        if (lastSegIndex < genomePtr->getSegmentNum() - 1) {
            curSegPtr = genomePtr->getSegment(lastSegIndex + 1);

            if (curSegPtr->getLength() > 0) {
                curNuPtr = curSegPtr->getNuAt(0);
            } else {
                curNuPtr = NULL;
                curSegPtr = NULL;
            }

        } else {
            curSegPtr = NULL;
            curNuPtr = NULL;
        }
    }
    
    return curNuPtr;
}

Nucleotide * const Genome::Cursor::previous(void) {
    if (curNuPtr == NULL) return NULL;
    
    unsigned long lastNuPos = curNuPtr->getBasePairPos();

    if (lastNuPos > 0) {
        curNuPtr = curSegPtr->getNuAt(lastNuPos - 1);

    } else {
        unsigned long lastSegIndex = genomePtr->getSegmentIndex(curSegPtr);
        if (lastSegIndex > 0) {
            curSegPtr = genomePtr->getSegment(lastSegIndex - 1);

            unsigned long newSegLength = curSegPtr->getLength();
            if (newSegLength > 0) {
                curNuPtr = curSegPtr->getNuAt(newSegLength - 1);
            } else {
                curNuPtr = NULL;
                curSegPtr = NULL;
            }

        } else {
            curSegPtr = NULL;
            curNuPtr = NULL;
        }
    }
    
    return curNuPtr;
}




////////////////////////////////////////////////////////////////////////////////
//      EXTRACT TEMPLATE CLASS
////////////////////////////////////////////////////////////////////////////////

//Genome::ExtractTemplate::ExtractTemplate(void) {
//}

//Genome::ExtractTemplate::ExtractTemplate(const ExtractTemplate& orig) {
//    nuPositions = orig.nuPositions;
//}

//Genome::ExtractTemplate& Genome::ExtractTemplate::operator=(const ExtractTemplate& rhs) {
//    if (this != &rhs) {
//        nuPositions = rhs.nuPositions;
//    }
//    return *this;
//}

//Genome::ExtractTemplate::~ExtractTemplate() {
//}

void Genome::ExtractTemplate::initialize(Genome* modelGenome) {
    unsigned long templateSize = modelGenome->getFullLength();
    unsigned long indexOnGenome = 0;

    nuPositions.clear();
    nuPositions.reserve(templateSize);

    for (unsigned long segIndex = 0; segIndex < modelGenome->getSegmentNum(); segIndex++) {
        unsigned long segmentLength = modelGenome->getSegment(segIndex)->getLength();

        for (unsigned long nuIndex = 0; nuIndex < segmentLength; nuIndex++) {
            NucleotidePos nucleotidePos;
            nucleotidePos.segmentIndex = segIndex;
            nucleotidePos.indexOnSegment = nuIndex;
            nucleotidePos.indexOnGenome = indexOnGenome;

            nuPositions.push_back(nucleotidePos);
            indexOnGenome++;
        }
    }
    assert(nuPositions.size() == templateSize);
}

void Genome::ExtractTemplate::shrink(unsigned long beginIndex, unsigned long endIndex) {
    assert(beginIndex >= 0 && beginIndex < nuPositions.size());
    assert(endIndex > beginIndex && endIndex <= nuPositions.size());

    unsigned long newSize = 0;
    for (unsigned long i = beginIndex; i < endIndex; i++) {
        nuPositions[newSize] = nuPositions[i];
        newSize++;
    }
    assert(newSize <= nuPositions.size());
    nuPositions.resize(newSize);
}

void Genome::ExtractTemplate::cutOut(unsigned long beginIndex, unsigned long endIndex) {
    assert(beginIndex >= 0L && beginIndex < nuPositions.size());
    assert(endIndex > beginIndex && endIndex <= nuPositions.size());

    vector<NucleotidePos>::iterator begin = nuPositions.begin() + beginIndex;
    vector<NucleotidePos>::iterator end = nuPositions.begin() + endIndex;

    nuPositions.erase(begin, end);
}

void Genome::ExtractTemplate::deleteCodon1n2(void) {
    const unsigned long cKeptCodon = 2; // 0-base position - keep only position 3
    unsigned long newSize = 0;

    for (unsigned long i = 0; i < nuPositions.size(); i++) {
        if (nuPositions[i].indexOnGenome % 3 == cKeptCodon) {
            nuPositions[newSize] = nuPositions[i];
            newSize++;
        }
    }

    assert(newSize <= nuPositions.size());
    nuPositions.resize(newSize);
}

void Genome::ExtractTemplate::deleteCodon3(void) {
    unsigned long cDeletedCodon = 2; // 0-base position
    unsigned long newSize = 0;

    for (unsigned long i = 0; i < nuPositions.size(); i++) {
        if (nuPositions[i].indexOnGenome % 3 == cDeletedCodon) {
            continue;
        }
        nuPositions[newSize] = nuPositions[i];
        newSize++;
    }

    assert(newSize <= nuPositions.size());
    nuPositions.resize(newSize);
}

void Genome::ExtractTemplate::keepPositions(const vector<bool>& keptMarkers) {
    if (nuPositions.size() == 0 || keptMarkers.size() == 0) {
        return;
    }

    assert(nuPositions.size() == keptMarkers.size());
    unsigned long newSize = 0;

    for (unsigned long i = 0; i < nuPositions.size(); i++) {
        if (keptMarkers[i]) {
            if (newSize != i) {
                nuPositions[newSize] = nuPositions[i];
            }
            newSize++;
        }
    }

    assert(newSize <= nuPositions.size());
    nuPositions.resize(newSize);
}
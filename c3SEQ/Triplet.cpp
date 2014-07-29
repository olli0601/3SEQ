/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   Triplet.cpp
 * CREATED ON:  24 June 2011, 15:57
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: The triplet of Genome objects.
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-06-24  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   24 June 2011
 */

#include "Triplet.h"
#include "config.h"

//Triplet::Triplet(Genome* newDad, Genome* newMum, Genome* newChild)
//: dad(newDad), mum(newMum), child(newChild) {
//
//    updateSteps();
//}

//Triplet::Triplet(const Triplet& orig)
//: dad(orig.getDad()), mum(orig.getMum()), child(orig.getChild()) {
//    updateSteps();
//}

//Triplet::~Triplet() {
//    clearBreakPoints();
//    // Do not delete other objects because this is just a temporary storage.
//}

//void Triplet::clearBreakPoints(void) {
//    if (breakPoints.size() > 0) {
//        for (unsigned long i = 0; i < breakPoints.size(); i++) {
//            if (breakPoints[i] != NULL) {
//                delete breakPoints[i];
//            }
//        }
//        breakPoints.clear();
//    }
//}

void Triplet::updateSteps(void) {
    assert(dad != NULL && mum != NULL && child != NULL);
    assert(child->getExtractedLength() == dad->getExtractedLength()
            && child->getExtractedLength() == mum->getExtractedLength());

    maxHeight = 0.0;
    upStep = 0;
    downStep = 0;
    maxDecent = 0;

    /* The following variable is declared in "double" type 
     * because "unsigned long" cannot store negative values */
    double currentHeight = 0.0;
    unsigned long length = child->getExtractedLength();

    Nucleotide* dadNu = NULL;
    Nucleotide* mumNu = NULL;
    Nucleotide* childNu = NULL;
    
    unsigned long maxDecentToThisPoint;

    for (unsigned long i = 0; i < length; i++) {
        dadNu = dad->getExtractedNu(i);
        mumNu = mum->getExtractedNu(i);
        childNu = child->getExtractedNu(i);

        if (!dadNu->isGap() && !mumNu->isGap() && !childNu->isGap()) {
            if (dadNu->getCode() == childNu->getCode()) {
                if (mumNu->getCode() != childNu->getCode()) {
                    /* Up step */
                    currentHeight += 1.0;
                    upStep++;
                }

            } else {
                if (mumNu->getCode() == childNu->getCode()) {
                    /* Down step */
                    currentHeight -= 1.0;
                    downStep++;
                }
            }

            if (currentHeight > maxHeight) {
                maxHeight = currentHeight;
            }

            maxDecentToThisPoint = static_cast<unsigned long> (maxHeight - currentHeight);

            if (maxDecent < maxDecentToThisPoint) {
                maxDecent = maxDecentToThisPoint;
                firstMaxDecentHeight = currentHeight;
            }
        }
    }

    isBreakPointsCalculated = false;

    pValue = PTable::NaPVAL;
    isApproxPVal = false;
}

void Triplet::deleteLastType1BreakPoints(void) {
    if (breakPoints.size() <= 0) return;

    unsigned long newSize = breakPoints.size();
    BreakPoint* lastBp = breakPoints[newSize - 1];

    while (newSize > 0 && lastBp->getType() == BreakPoint::FIRST) {
        delete lastBp;
        newSize--;
        lastBp = breakPoints[newSize - 1];
    }

    breakPoints.resize(newSize, NULL);
}

void Triplet::seekBreakPoints(void) {
    assert(maxDecent > 0);

    clearBreakPoints();

    double currentHeight = 0.0;
    double lastTypeOneBpHeight = firstMaxDecentHeight + static_cast<double> (maxDecent);
    BreakPoint* lastBp = NULL;

    if (lastTypeOneBpHeight <= 0) {
        lastBp = new BreakPoint(child, BreakPoint::FIRST, NULL, 0.0);
        breakPoints.push_back(lastBp);
    }

    Genome::Cursor* childCursor = child->newCursor();

    for (unsigned long i = 0; i < child->getExtractedLength(); i++) {
        Nucleotide* dadNu = dad->getExtractedNu(i);
        Nucleotide* mumNu = mum->getExtractedNu(i);
        Nucleotide* childNu = child->getExtractedNu(i);

        if (!dadNu->isGap() && !mumNu->isGap() && !childNu->isGap()) {

            if (dadNu->getCode() == childNu->getCode()
                    && mumNu->getCode() != childNu->getCode()) {
                /* Up step */

                if (lastBp != NULL && lastBp->getType() == BreakPoint::SECOND && lastBp->isOpen()) {
                    /* Still not modify the currentHeight */
                    assert(currentHeight == lastBp->getHeight());

                    childCursor->setPosition(childNu);
                    lastBp->setRightBound(childCursor->previous());
                }

                /* Increase currentHeight */
                currentHeight += 1.0;

                /* If a new highest point detected */
                if (lastTypeOneBpHeight < currentHeight) {
                    lastTypeOneBpHeight = currentHeight;
                    deleteLastType1BreakPoints();
                }

                /* If this is one of the highest point (up to this place), save it */
                if (lastTypeOneBpHeight == currentHeight) {
                    lastBp = new BreakPoint(child, BreakPoint::FIRST, childNu, currentHeight);
                    breakPoints.push_back(lastBp);
                }

            } else if (dadNu->getCode() != childNu->getCode()
                    && mumNu->getCode() == childNu->getCode()) {
                /* Down step */

                if (lastBp != NULL && lastBp->isOpen()) {
                    /* Still not modify the currentHeight */
                    assert(lastBp->getType() == BreakPoint::FIRST && currentHeight == lastBp->getHeight());

                    childCursor->setPosition(childNu);
                    lastBp->setRightBound(childCursor->previous());
                }

                /* Decrease currentHeight */
                currentHeight -= 1.0;

                /* If maximum decent reached, save the point */
                if (breakPoints.size() > 0 && lastTypeOneBpHeight - currentHeight == maxDecent) {
                    lastBp = new BreakPoint(child, BreakPoint::SECOND, childNu, currentHeight);
                    breakPoints.push_back(lastBp);
                }
            }
        }
    }

    deleteLastType1BreakPoints();

    if (breakPoints.size() > 0) {
        lastBp = breakPoints[breakPoints.size() - 1];
        assert(lastBp->getType() == BreakPoint::SECOND);

        /* If the last breakpoint has not had the right bound, find it */
        if (lastBp->isOpen()) {
            childCursor->putAtLastNu();
            lastBp->setRightBound(childCursor->current());
        }
    }

    delete childCursor;
}

void Triplet::seekBreakPointPairs(void) {
    unsigned long newRecLength;

    Genome::Cursor* childCursor = child->newCursor();
    Nucleotide* tmpLeftNu;
    Nucleotide* tmpRightNu;

    clearBreakPoints();
    breakPointPairs.clear();

    minRecLength = child->getFullLength();
    //    maxRecLength = 0;

    seekBreakPoints();

    for (unsigned long i = 0; i < breakPoints.size(); i++) {
        if (breakPoints[i]->getType() == BreakPoint::SECOND) {
            double acceptedHeight = breakPoints[i]->getHeight() + maxDecent;
            unsigned long currentIndex = i;

            while (currentIndex > 0) {
                currentIndex--;
                if (breakPoints[currentIndex]->getType() == BreakPoint::SECOND) {
                    continue;
                } else {
                    if (breakPoints[currentIndex]->getHeight() < acceptedHeight) {
                        break; // finish while-loop
                    } else {
                        assert(breakPoints[currentIndex]->getHeight() == acceptedHeight);

                        /* Create breakpoint pair */
                        BreakPointPair bpPair;
                        bpPair.firstBp = breakPoints[currentIndex];
                        bpPair.secondBp = breakPoints[i];
                        breakPointPairs.push_back(bpPair);


                        /* Update minRecLength */
                        newRecLength = countNonGappedSiteInFull(
                                bpPair.firstBp->getRightBound(),
                                bpPair.secondBp->getLeftBound());
                        if (minRecLength > newRecLength) minRecLength = newRecLength;

                        childCursor->setPosition(bpPair.firstBp->getLeftBound());
                        tmpLeftNu = childCursor->previous();

                        childCursor->setPosition(bpPair.secondBp->getRightBound());
                        tmpRightNu = childCursor->next();

                        newRecLength = countNonGappedSiteInFull(NULL, tmpLeftNu)
                                + countNonGappedSiteInFull(tmpRightNu, NULL);
                        if (minRecLength > newRecLength) minRecLength = newRecLength;


                        /* Update maxRecLength */
                        //newRecLength = countNonGappedSiteInFull(
                        //        bpPair.firstBp->getLeftBound(),
                        //        bpPair.secondBp->getRightBound());
                        //if (maxRecLength < newRecLength) maxRecLength = newRecLength;
                        //
                        //childCursor->setPosition(bpPair.firstBp->getRightBound());
                        //tmpLeftNu = childCursor->previous();
                        //
                        //childCursor->setPosition(bpPair.secondBp->getLeftBound());
                        //tmpRightNu = childCursor->next();
                        //
                        //newRecLength =countNonGappedSiteInFull(NULL, tmpLeftNu)
                        //        + countNonGappedSiteInFull(tmpRightNu, NULL);
                        //if (maxRecLength < newRecLength) maxRecLength = newRecLength;                        
                    }
                }
            }
        }
    }

    delete childCursor;

    isBreakPointsCalculated = true;
}

unsigned long Triplet::countNonGappedSiteInFull(Nucleotide* leftChildNu,
        Nucleotide* rightChildNu) const {

    if (leftChildNu == NULL && rightChildNu == NULL)
        return 0;


    Genome::Cursor* childCursor = child->newCursor();
    Genome::Cursor* dadCursor = dad->newCursor();
    Genome::Cursor* mumCursor = mum->newCursor();

    if (leftChildNu == NULL) {
        childCursor->putAtFirstNu();
        leftChildNu = childCursor->current();
    }

    if (rightChildNu == NULL) {
        childCursor->putAtLastNu();
        rightChildNu = childCursor->current();
    }

    assert(leftChildNu != NULL && rightChildNu != NULL);

    unsigned long leftSegIndex = child->getSegmentIndex(leftChildNu->getContainingSequence());
    unsigned long leftNuIndex = leftChildNu->getBasePairPos();

    childCursor->setPosition(leftChildNu);
    dadCursor->setPosition(leftSegIndex, leftNuIndex);
    mumCursor->setPosition(leftSegIndex, leftNuIndex);

    unsigned long numNonGapped = 0;

    if (!dadCursor->current()->isGap() && !mumCursor->current()->isGap()
            && !childCursor->current()->isGap()) {
        numNonGapped++;
    }

    while (childCursor->current() != rightChildNu && childCursor->next() != NULL) {
        dadCursor->next();
        mumCursor->next();
        if (!dadCursor->current()->isGap() && !mumCursor->current()->isGap()
                && !childCursor->current()->isGap()) {
            numNonGapped++;
        }
    }

    delete childCursor;
    delete dadCursor;
    delete mumCursor;

    return numNonGapped;
}

const bool Triplet::calculatePVal(bool acceptApprox) {
    if (PTable::instance().canCalculateExact(upStep, downStep, maxDecent)) {
        pValue = PTable::instance().getExactPValue(upStep, downStep, maxDecent);
        isApproxPVal = false;

    } else if (acceptApprox) {
        pValue = PTable::instance().approxPValue(upStep, downStep, maxDecent);
        isApproxPVal = true;

    } else {
        pValue = PTable::NaPVAL;
        isApproxPVal = false;
        return false;
    }

    return true;
}

const double Triplet::getPValSingleBp(void) const {
    double closedPVal = 1.0;
    double dUpStep = static_cast<double> (upStep);
    double dDownStep = static_cast<double> (downStep);

    for (double loopCount = 1.0; loopCount <= maxHeight; loopCount += 1.0) {
        closedPVal = closedPVal * (dUpStep / (dDownStep + maxHeight));
        dUpStep = dUpStep - 1.0;
        dDownStep = dDownStep - 1.0;
    }

    return closedPVal;
}

string Triplet::toString(void) const {
    stringstream sStream;
    sStream << dad->getAccession() << "\t"
            << mum->getAccession() << "\t"
            << child->getAccession() << "\t"
            << upStep << "\t"
            << downStep << "\t"
            << maxDecent;
    return sStream.str();
}

void Triplet::generatePSFile(string fileName) const {
    /** Width of each nucleotides space in the PS file. */
    static const unsigned long WIDTH = 2;

    /** Height of the image. */
    static const unsigned long HEIGHT = 240;

    TextFile psFile(fileName);

    if (psFile.exists()) {
        Interface::instance()
                << "The file \"" << psFile.getPath() << "\" already exists.\n"
                << "Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* If the user say "No" */
            return;

        } else {
            psFile.removeFile();
        }
    }


    psFile.openToWrite();

    char line[200];
    psFile.writeLine("%%!PS-Adobe-3.0 EPSF-3.0");
    sprintf(line, "%%%%BoundingBox: 0 0 %lu %lu", (child->getFullLength() + 1) * WIDTH, HEIGHT);
    psFile.writeLine(line);
    psFile.writeLine("%%%%EndComments\n");

    psFile.writeLine("%%%%%%%%%%%%");
    psFile.writeLine("/Times-Roman findfont");
    psFile.writeLine("7 scalefont");
    psFile.writeLine("setfont\n");

    psFile.writeLine("0.001 setlinewidth");
    psFile.writeLine();

    sprintf(line, "/w {%lu} def    %% This is the width of a small box on the plot\n", WIDTH);
    psFile.writeLine(line);

    sprintf(line, "/h {%lu} def    %% This is the height of a small box on the plot\n", HEIGHT);
    psFile.writeLine(line);


    unsigned long segmentNum = child->getSegmentNum();
    double numPreNu = 0.0; // the total number of nucleotides of the previous segment

    for (unsigned long segIndex = 0; segIndex < segmentNum; segIndex++) {
        for (unsigned long nuIndex = 0; nuIndex < child->getSegment(segIndex)->getLength(); nuIndex++) {

            Nucleotide* dadNu = dad->getSegment(segIndex)->getNuAt(nuIndex);
            Nucleotide* mumNu = mum->getSegment(segIndex)->getNuAt(nuIndex);
            Nucleotide* childNu = child->getSegment(segIndex)->getNuAt(nuIndex);

            double indexInFull = numPreNu + static_cast<double> (nuIndex);

            if (!dadNu->isGap() || !mumNu->isGap() || !childNu->isGap()) {
                if (dadNu->getCode() == childNu->getCode() && mumNu->getCode() != childNu->getCode()) {
                    sprintf(line, "1 0 0 setrgbcolor %0.0lf 3 moveto", indexInFull * static_cast<double> (WIDTH));
                    psFile.writeLine(line);
                    psFile.writeLine("0 h rlineto w 0 rlineto 0 h neg rlineto closepath fill stroke\n");

                } else if (dadNu->getCode() != childNu->getCode() && mumNu->getCode() == childNu->getCode()) {
                    sprintf(line, "0 0 1 setrgbcolor %0.0lf 3 moveto", indexInFull * static_cast<double> (WIDTH));
                    psFile.writeLine(line);
                    psFile.writeLine("0 h rlineto w 0 rlineto 0 h neg rlineto closepath fill stroke\n");

                } else if (dadNu->getCode() != childNu->getCode() || mumNu->getCode() != childNu->getCode()) {
                    /* Draw a gray line if there is polymorphism but it is not an informative site */
                    sprintf(line, "0.4 0.4 0.4 setrgbcolor %0.0lf 3 moveto", indexInFull * static_cast<double> (WIDTH));
                    psFile.writeLine(line);

                    sprintf(line, "0 %lu rlineto w 0 rlineto 0 h neg rlineto closepath fill stroke\n", HEIGHT / 3);
                    psFile.writeLine(line);
                }
            }
        }

        numPreNu += static_cast<double> (child->getSegment(segIndex)->getLength());
    }

    psFile.writeLine("showpage");
    psFile.writeLine("%%EOF");

    psFile.close();
}




////////////////////////////////////////////////////////////////////////////////
//  BREAKPOINT PAIR
////////////////////////////////////////////////////////////////////////////////

const string Triplet::BreakPointPair::toString() const {
    stringstream tmpStream;
    tmpStream << firstBp->getLeftBoundPos() << "-"
            << firstBp->getRightBoundPos() << " & "
            << secondBp->getLeftBoundPos() << "-"
            << secondBp->getRightBoundPos();
    return tmpStream.str();
}
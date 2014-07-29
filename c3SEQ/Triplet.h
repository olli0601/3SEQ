/* 
 * File:   Triplet.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 24 June 2011, 15:57
 */

#ifndef TRIPLET_H
#define	TRIPLET_H

#include <sstream>
#include <string>
#include <vector>
#include "Genome.h"
#include "Nucleotide.h"
#include "PTable.h"

using namespace std;


class Triplet {
public:

    class BreakPoint;

    struct BreakPointPair {
        BreakPoint* firstBp;
        BreakPoint* secondBp;

        const string toString(void) const;
    };

    explicit Triplet(Genome* newDad, Genome* newMum, Genome* newChild)
    : dad(newDad), mum(newMum), child(newChild) {
        updateSteps();
    }

    ~Triplet() {
        clearBreakPoints();
        // Do not delete other objects because this is just a temporary storage.
    }

    Genome * const getDad(void) const {
        return dad;
    }

    Genome * const getMum(void) const {
        return mum;
    }

    Genome * const getChild(void) const {
        return child;
    }

    /**
     * Get the number of up-step.
     * @return The number of up-step.
     */
    unsigned long const getUpStep(void) const {
        return upStep;
    }

    /**
     * Get the number of down-step.
     * @return The number of down-step.
     */
    unsigned long const getDownStep(void) const {
        return downStep;
    }

    /**
     * Get the maximum decent.
     * @return The maximum decent.
     */
    unsigned long const getMaxDecent(void) const {
        return maxDecent;
    }

    /**
     * Get the P-value of this triplet.
     * @return  The P-value of this triplet.
     * @note    Be careful to make sure that the P-value can be calculated or
     *          approximated before calling this method.
     */
    const double getPVal(void) const {
        assert(pValue != PTable::NaPVAL);
        return pValue;
    }

    /**
     * Check to see if the P-value is approximated or exactly calculated.
     * @return <b>true</b> if the P-value is approximated.
     */
    const bool isPValApproximated(void) const {
        return isApproxPVal;
    }

    /**
     * Get the minimum recombinant segment lenght.
     * @return  The minimum recombinant segment lenght.
     */
    const unsigned long getMinRecLength(void) {
        if (!isBreakPointsCalculated) {
            seekBreakPointPairs();
        }
        return minRecLength;
    }

    /**
     * Get the maximum recombinant segment lenght.
     * @return  The maximum recombinant segment lenght.
     */
    //    const unsigned long getMaxRecLength(void) {
    //        if (!isBreakPointsCalculated) {
    //            seekBreakPointPairs();
    //        }
    //        return maxRecLength;
    //    }

    vector<BreakPointPair> getBreakPointPairs(void) {
        if (!isBreakPointsCalculated) {
            seekBreakPointPairs();
        }
        return breakPointPairs;
    }

    /**
     * Calculate the P-value for this triplet
     * @param acceptApprox  indicates that approximation can be used if the 
     *                      P-value cannot be calculated exactly.
     * @return  <b>true</b> if the P-value can be calculated or approximated;
     *          <b>false</b>, otherwise.
     */
    const bool calculatePVal(bool acceptApprox);

    /**
     * Get the P-value for the single break-point (the place where the height
     * of the random walk reaches the highest value).
     * @return The P-value for the single break-point.
     */
    const double getPValSingleBp(void) const;

    /**
     * Print out all the triplet info in 1 string.
     * @return A string that represents all the important info of the triplet.
     */
    string toString(void) const;

    /**
     * Generate a postscript file that represents the triplet.
     * @param fileName
     */
    void generatePSFile(string fileName) const;


private:

    Triplet(const Triplet& orig);
    //    Triplet(const Triplet& orig)
    //    : dad(orig.getDad()), mum(orig.getMum()), child(orig.getChild()) {
    //        updateSteps();
    //    }

    Triplet& operator=(const Triplet& rhs);
    //    Triplet& operator=(const Triplet& rhs) {
    //        assert(false);
    //    }

    /**
     * Delete all the element of the <b>breakPoints</b> vector, then clear it.
     */
    void clearBreakPoints(void);
    //    {
    //        if (breakPoints.size() > 0) {
    //            for (unsigned long i = 0; i < breakPoints.size(); i++) {
    //                if (breakPoints[i] != NULL) {
    //                    delete breakPoints[i];
    //                }
    //            }
    //            breakPoints.clear();
    //        }
    //    }

    void updateSteps(void);

    void seekBreakPoints(void);

    void seekBreakPointPairs(void);

    /**
     * Delete last type-one breakpoints in the breakpoint list to make sure that
     * the last breakpoint in the list is type-two.
     */
    void deleteLastType1BreakPoints(void);

    /**
     * Count the total number of non-gapped sites ranged between the 
     * <b>leftChildNu</b> and the <b>rightChildNu</b> (inclusive).
     * @param leftChildNu   The left bound of the range. If this value is
     *                      <b>NULL</b> then the 1st nucleotide of the child
     *                      sequence is considered as the left bound.
     * @param rightChildNu  The right bound of the range. If this value is
     *                      <b>NULL</b> then the last nucleotide of the child
     *                      sequence is considered as the right bound.
     * @return  The total number of non-gapped sites ranged between the 
     *          <b>leftChildNu</b> and the <b>rightChildNu</b> (inclusive).
     * @note    If both left bound and right bound are <b>NULL</b>, the returned
     *          value will be 0;
     *          We take into account all nucleotides in the range, not only the
     *          extracted ones.<br>
     *          A site is considered as non-gapped if and only if all the
     *          nucleotides at that site in the 3 sequences are not gap.
     */
    unsigned long countNonGappedSiteInFull(Nucleotide* leftChildNu, Nucleotide* rightChildNu) const;

    Genome * const dad;
    Genome * const mum;
    Genome * const child;

    unsigned long upStep;
    unsigned long downStep;
    unsigned long maxDecent;

    /**
     * Maximum height of the random walk. Use "double" instead of 
     * "unsigned long" since this value can be negative.
     */
    double maxHeight;

    /**
     * The height of the 1st maximum decent. This is useful to improve the
     * speed of breakpoint searching method.
     */
    double firstMaxDecentHeight;

    unsigned long minRecLength;
    //    unsigned long maxRecLength;

    double pValue;
    bool isApproxPVal;

    vector<BreakPoint*> breakPoints;
    vector<BreakPointPair> breakPointPairs;
    bool isBreakPointsCalculated;
};




////////////////////////////////////////////////////////////////////////////////
//  BREAKPOINT
////////////////////////////////////////////////////////////////////////////////

class Triplet::BreakPoint {
public:

    enum Type {
        FIRST, SECOND
    };

    explicit BreakPoint(Genome * const newChildGenome, const Type newType,
            Nucleotide * const newLeftBound, double newHeight)
    : childGenome(newChildGenome), type(newType) {

        leftBound = newLeftBound;
        rightBound = newLeftBound;
        height = newHeight;

        isClosed = false;
        isBpPosCalculated = false;
    }

    BreakPoint(const BreakPoint& orig)
    : childGenome(orig.childGenome), type(orig.type) {

        leftBound = orig.leftBound;
        rightBound = orig.rightBound;
        height = orig.height;
        isClosed = orig.isClosed;
    }

    BreakPoint& operator=(const BreakPoint& rhs) {
        if (this != &rhs) {
            assert(type == rhs.type);
            leftBound = rhs.leftBound;
            rightBound = rhs.rightBound;
            height = rhs.height;
            isClosed = rhs.isClosed;
        }
        return (*this);
    }

    ~BreakPoint() {
        // Do nothing
    }

    const Type getType(void) const {
        return type;
    }

    Nucleotide * const getLeftBound(void) const {
        return leftBound;
    }

    Nucleotide * const getRightBound(void) const {
        assert(isClosed);
        return rightBound;
    }

    const double getHeight(void) const {
        return height;
    }

    const bool isOpen(void) const {
        return !isClosed;
    }

    void setRightBound(Nucleotide * const newRightBound) {
        assert(!isClosed);
        rightBound = newRightBound;
        isClosed = true;
    }

    const unsigned long getLeftBoundPos(void) {
        if (!isBpPosCalculated) {
            calculateBpPos();
        }
        return leftBoundPos;
    }

    const unsigned long getRightBoundPos(void) {
        if (!isBpPosCalculated) {
            calculateBpPos();
        }
        return rightBoundPos;
    }


private:

    void calculateBpPos(void) {
        if (leftBound == NULL)
            leftBoundPos = 0;
        else
            leftBoundPos = childGenome->getBasePairPosOf(leftBound) + 1;

        if (rightBound == NULL)
            rightBoundPos = 0;
        else
            rightBoundPos = childGenome->getBasePairPosOf(rightBound) + 1;

        isBpPosCalculated = true;
    }

    Genome * const childGenome;
    const Type type;
    Nucleotide* leftBound;
    Nucleotide* rightBound;

    unsigned long leftBoundPos;
    unsigned long rightBoundPos;
    bool isBpPosCalculated;

    /**
     * The height of the breakpoint.<br>
     * This is also an integer. Use "double" to extend range.
     */
    double height;

    /** Indicates that the right bound of this breakpoint has been set. */
    bool isClosed;

};




////////////////////////////////////////////////////////////////////////////////
//  INLINE FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

inline void Triplet::clearBreakPoints(void) {
    if (breakPoints.size() > 0) {
        for (unsigned long i = 0; i < breakPoints.size(); i++) {
            if (breakPoints[i] != NULL) {
                delete breakPoints[i];
            }
        }
        breakPoints.clear();
    }
}

#endif	/* TRIPLET_H */


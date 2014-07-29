/* 
 * File:   Genome.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 27 July 2011, 16:33
 */

#ifndef GENOME_H
#define	GENOME_H

#include <cassert>
#include <vector>
#include <string>
#include "Sequence.h"
#include "Nucleotide.h"

class Triplet;

/**
 * For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class Genome : public Sequence {
public:

    enum RecombinantType {
        Na_REC /** Not a recombinant */,
        SHORT_REC /** Short recombinant */,
        LONG_REC /** Long recombinant */
    };

    class ExtractTemplate;

    class Cursor;

    explicit Genome(const string& accession) : Sequence(accession) {
        bestRecombinantTriplet = NULL;
        recombinantType = Na_REC;
        isActivated = true;
    }

    ~Genome();

    /**
     * Indicates if this genome is used for the analysis.
     * @return <b>true</b> if this genome is used.
     */
    bool isActive(void) {
        return isActivated;
    }
    
    /**
     * Indicates if this genome is used for the analysis.
     * @param active
     */
    void setActive(bool active) {
        isActivated = active;
    }
    
    /**
     * Get the accession of this genome.
     * @return the accession of this genome.
     * @note This method is an alias of <b>getName()</b>.
     */
    string const getAccession(void) const {
        return getName();
    }

    void setRecombinantType(RecombinantType newType) {
        recombinantType = newType;
    }

    const RecombinantType getRecombinantType(void) const {
        return recombinantType;
    }

    /**
     * Add a segment into the end of this genome.
     * @param segment
     */
    void addSegment(Sequence * segment) {
        segments.push_back(segment);
    }

    /**
     * Get the number of segments contained in this genome.
     * @return  The number of segments contained in this genome.
     */
    unsigned long const getSegmentNum(void) const {
        return segments.size();
    }

    /**
     * Access to the segment at the given index.
     * @param index
     * @return The segment pointer.
     */
    Sequence * const getSegment(const unsigned long& index) const {
        assert(index >= 0 && index < segments.size());
        return segments[index];
    }

    /**
     * Search for a segment by name.
     * @param segmentName
     * @return  A pointer that points to the 1st segment found. If no segment 
     *          is found, a NULL pointer will be returned.
     */
    Sequence * const getSegment(const string& segmentName) const {
        for (unsigned long i = 0; i < segments.size(); i++) {
            if (segments[i]->hasName(segmentName)) {
                return segments[i];
            }
        }
        return NULL;
    }

    /**
     * Get the nucleotide at the given index from the extracted part (vector).
     * @param index The index of the nucleotide in the vector (not the real position).
     * @return  The nucleotide.
     * @note    This is an alias of <b>getNuAt(index)</b> method.
     */
    Nucleotide * const getExtractedNu(unsigned long index) const {
        //        return Sequence::getNuAt(index);
        assert(index >= 0 && index < nucleotides.size());
        return nucleotides[index];
    }

    /**
     * Get the nucleotide at the given index from the extracted part (vector).
     * @param index The index of the nucleotide in the vector (not the real position).
     * @return  The nucleotide.
     * @note    This is an alias of <b>getExtractedNu(index)</b> method.
     */
    Nucleotide * const getNuAt(unsigned long index) const {
        return Sequence::getNuAt(index);
    }

    /**
     * Get the length of the extracted part of this genome.
     * @return  The length of the extracted part of this genome.
     * @note    This is an alias of <b>getLength()</b> method.
     */
    unsigned long const getExtractedLength(void) const {
        return Sequence::getLength();
    }

    /**
     * Get the length of the extracted part of this genome.
     * @return The length of the extracted part of this genome.
     * @note    This is an alias of <b>getExtractedLength()</b> method.
     */
    unsigned long const getLength(void) const {
        return Sequence::getLength();
    }

    /**
     * Get the full length (i.e the sum of lengths of all the segments) of this genome.
     * @return The full length of this genome.
     */
    unsigned long const getFullLength(void) const {
        unsigned long fullLength = 0;
        for (unsigned long i = 0; i < segments.size(); i++) {
            fullLength += segments[i]->getLength();
        }
        return fullLength;
    }

    /**
     * Get the index of the given segment.
     * @param segment
     * @return  The index of the given segment.
     * @note    If this genome does not contain the segment, the program will
     *          be terminated with an assertion fail.
     */
    unsigned long getSegmentIndex(Sequence* segment) const {
        assert(segment != NULL);

        for (unsigned long i = 0; i < segments.size(); i++) {
            if (segment == segments[i]) {
                return i;
            }
        }
        /* Just stop the program if cannot find the segment because this 
         * should never happen. */
        assert(false);
    }

    void setBestRecombinantTriplet(Triplet* triplet);

    Triplet* getBestRecombinantTriplet(void) {
        return bestRecombinantTriplet;
    }

    /**
     * Create a new cursor for this genome.
     * @return  The new cursor.
     * @note    This cursor must be manually deleted after use.
     */
    Cursor * const newCursor(void) const;

    /**
     * Based on a template (vector) which contains the positions of the 
     * nucleotides that are needed to analyse, extract the nucleotides 
     * (i.e put them into the extracted vector).
     * @param nuPositions   The vector which contains the positions of all 
     *                      nucleotides that are needed to be extracted.
     */
    void extractByTemplate(const ExtractTemplate& extractTemplate);

    /**
     * Calculate the base-pair position of the given nucleotide.
     * @param nucleotide
     * @return  The base-pair position of the given nucleotide.
     * @note    This position is 0-based.
     */
    unsigned long getBasePairPosOf(Nucleotide* nucleotide) const;

    /**
     * Calculate the base-pair distance between 2 given nucleotides.
     * @param firstNu
     * @param secondNu
     * @return  The base-pair distance between 2 given nucleotides. If the
     *          <b>firstNu</b> is on the left side of the <b>secondNu</b>,
     *          the return distance will be positive; otherwise, a negative
     *          distance will be returned.
     * @note    The distance depends on the order of the segments.<br>
     *          If any of the 2 given nucleotides is not in this genome, 
     *          the program will be suspended with an error.<br>
     *          The return value is "double" type since "unsigned long"
     *          cannot store negative values.<br>
     *          The distance is calculated based on the full genome (not only
     *          the extracted part).
     */
    double bpDistanceBetween(Nucleotide* firstNu, Nucleotide* secondNu) const;


private:
    Genome(const Genome& orig);

    Genome& operator=(const Genome& rhs);

    /** 
     * The storage of all the segments in this genome. 
     * @note    Not all parts of all segments will be taken into the analysis.
     *          Just the nucleotides that are in the <b>nucleotides</b> vector 
     *          will be used. We call that as extracted part.
     */
    vector<Sequence*> segments;

    /**
     * The best triplet where this genome is the child sequence. If this genome
     * is not a recombinant, this pointer will be NULL.
     */
    Triplet* bestRecombinantTriplet;

    
    RecombinantType recombinantType;
    
    /** Indicates if this genome is used for the analysis. */
    bool isActivated;
};




////////////////////////////////////////////////////////////////////////////////
//      CURSOR CLASS
////////////////////////////////////////////////////////////////////////////////

/**
 * The cursor that runs through all the nucleotides in the genome (not only the
 * extracted part).
 * @note    For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class Genome::Cursor {
public:
    /**
     * Initialise the cursor and put it at the first nucleotide.
     * @param newGenomePtr
     */
    explicit Cursor(Genome const * const newGenomePtr);

    void setPosition(unsigned long segIndex, unsigned long nuPosInSeg);

    void setPosition(Nucleotide * const nuPtr);

    /**
     * Put the cursor at the first nucleotide in the genome.
     */
    void putAtFirstNu(void);

    /**
     * Put the cursor at the last nucleotide in the genome.
     */
    void putAtLastNu(void);

    /**
     * Get the nucleotide which is currently pointed by the cursor.
     * @return  The currently pointed nucleotide.
     */
    Nucleotide * const current(void) const;

    /**
     * Jump to the next nucleotide.
     * @return  The pointer to the next nucleotide. If the current
     *          nucleotide is already at the last position, a <b>NULL</b>
     *          pointer will be returned.
     */
    Nucleotide * const next(void);

    /**
     * Jump to the previous nucleotide.
     * @return  The pointer to the previous nucleotide. If the current
     *          nucleotide is already at the first position, a <b>NULL</b>
     *          pointer will be returned.
     */
    Nucleotide * const previous(void);

private:
    Genome const * const genomePtr;
    Nucleotide * curNuPtr;
    Sequence* curSegPtr;
};




////////////////////////////////////////////////////////////////////////////////
//      EXTRACT TEMPLATE CLASS
////////////////////////////////////////////////////////////////////////////////

/**
 * For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class Genome::ExtractTemplate {
public:

    explicit ExtractTemplate(void) {
        // Do nothing
    }

    ExtractTemplate(const ExtractTemplate& orig) {
        nuPositions = orig.nuPositions;
    }

    ExtractTemplate& operator=(const ExtractTemplate& rhs) {
        if (this != &rhs) {
            nuPositions = rhs.nuPositions;
        }
        return *this;
    };

    ~ExtractTemplate() {
        // Do nothing
    }

    unsigned long getSize(void) const {
        return nuPositions.size();
    }

    unsigned long const segmentIndexAt(unsigned long templateIndex) const {
        assert(templateIndex >= 0 && templateIndex < nuPositions.size());
        return nuPositions[templateIndex].segmentIndex;
    }

    unsigned long const nuIndexOnSegmentAt(unsigned long templateIndex) const {
        assert(templateIndex >= 0 && templateIndex < nuPositions.size());
        return nuPositions[templateIndex].indexOnSegment;
    }

    /**
     * Initialise the template.
     */
    void initialize(Genome* modelGenome);

    /**
     * Delete all BUT the position (vector-position) from <b>beginIndex</b> to 
     * <b>endIndex - 1</b> in the template.
     * @param beginIndex
     * @param endIndex
     * @note    Make sure that the <b>extractTemplate</b> vector was initialised 
     *          before calling this method.
     */
    void shrink(unsigned long beginIndex, unsigned long endIndex);

    /**
     * Delete all elements that are in the position (vector-position) range from 
     * <b>beginIndex</b> to <b>endIndex - 1</b> in the template.
     * @param beginIndex
     * @param endIndex
     * @note    Make sure that the <b>extractTemplate</b> vector was initialised 
     *          before calling this method.
     */
    void cutOut(unsigned long beginIndex, unsigned long endIndex);

    /**
     * Delete all codon at position 1 and 2.
     */
    void deleteCodon1n2(void);

    /**
     * Delete all codon at position 3.
     */
    void deleteCodon3(void);

    /**
     * Remove all positions that are not marked as kept.
     * @param isPositionKeeped  The vector contains boolean values that indicate
     *                          which position should be kept.
     */
    void keepPositions(const vector<bool>& keptMarkers);

private:

    struct NucleotidePos {
        unsigned long segmentIndex;
        unsigned long indexOnSegment;
        unsigned long indexOnGenome;
    };

    vector<NucleotidePos> nuPositions;       
};

#endif	/* GENOME_H */


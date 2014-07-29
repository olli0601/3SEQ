/**
 * FILE NAME:   GenomeSet.cpp
 * CREATED ON:  17 June 2011, 11:15
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: See "GenomeSet.h"
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-06-17  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   17 June 2011
 */

#include "GenomeSet.h"

//GenomeSet::GenomeSet(void) {
//    nGap = 0;
//    nMonoallelic = 0;
//    nBiallelic = 0;
//    nTriallelic = 0;
//    nTetrallelic = 0;
//    nNoGap = 0;
//}

//GenomeSet::GenomeSet(const GenomeSet& orig) : extractTemplate(orig.extractTemplate) {
//    polymorphicMarkers = orig.polymorphicMarkers;
//    genomes = orig.genomes; // shallow copying is enough.
//
//    nGap = orig.nGap;
//    nMonoallelic = orig.nMonoallelic;
//    nBiallelic = orig.nBiallelic;
//    nTriallelic = orig.nTriallelic;
//    nTetrallelic = orig.nTetrallelic;
//    nNoGap = orig.nNoGap;
//}

//GenomeSet::~GenomeSet() {
//    for (unsigned i = 0; i < genomes.size(); i++) {
//        if (genomes[i]) {
//            delete genomes[i];
//        }
//    }
//}

unsigned long const GenomeSet::getDistinctGenomes(const bool autoRemoveDuplicate) {
    return getNonNeighborGenomes(0, autoRemoveDuplicate);
}

unsigned long const GenomeSet::getNonNeighborGenomes(
        const unsigned long maxNeighborDistance,
        const bool autoRemoveNeighbor) {

    unsigned long nonNeighborNum = 0;

    /* The array contain flags to indicate if the genomes should be delete or not. 
     * We use this to reduce the number of comparisons among genomes. */
    vector<bool> markedForDelete;

    /* Populate all the elements with false value. */
    markedForDelete.resize(getSize(), false);

    /* The outer loop must run to getSize()-1 , cannot be getSize()-2 */
    for (unsigned long i = 0; i < getSize(); i++) {
        if (!markedForDelete[i]) {
            if (autoRemoveNeighbor && nonNeighborNum != i) {
                /* Shift the genomes vector */
                genomes[nonNeighborNum] = genomes[i];
            }
            nonNeighborNum++;

            for (unsigned long k = i + 1; k < getSize(); k++) {
                if (!markedForDelete[k]
                        && genomes[i]->ntDistanceTo(genomes[k], true) <= maxNeighborDistance) {
                    markedForDelete[k] = true;
                }
            }
        }
    }

    if (autoRemoveNeighbor && nonNeighborNum < genomes.size()) {
        /* Shrink the genomes vector */
        genomes.resize(nonNeighborNum, NULL);
    }

    return nonNeighborNum;
}

void GenomeSet::updateStats(void) {
    nGap = 0;
    nMonoallelic = 0;
    nBiallelic = 0;
    nTriallelic = 0;
    nTetrallelic = 0;
    nNoGap = 0;

    if (genomes.size() <= 0) return;

    /* We assume that the all the extracted genomes in this set 
     * have the same length. */
    unsigned long extractedLength = getExtractedLength();

    polymorphicMarkers.resize(extractedLength, false);
    nNoGap = extractedLength;

    /* Allocate space for the array that hold integers between 0 to 31 
     * showing which NT are present at each position (by using bit mask) 
     * to detect if the positions is polymorphic or not. */
    vector<char> siteType; // 8 bits is enough
    siteType.resize(extractedLength, 0);

    /* Bit mask: GAP = 10000; A = 01000; C = 00100; G = 00010; T = 00001 */
    const char GAP = 16; /* Bit mask for gap */
    const char A = 8; /* Bit mask for allele A */
    const char C = 4; /* Bit mask for allele C */
    const char G = 2; /* Bit mask for allele G */
    const char T = 1; /* Bit mask for allele T */

    char nuCode;

    for (unsigned long nuIndex = 0; nuIndex < extractedLength; nuIndex++) {
        siteType[nuIndex] = 0;

        for (unsigned long seqIndex = 0; seqIndex < genomes.size(); seqIndex++) {
            nuCode = genomes[seqIndex]->getExtractedNu(nuIndex)->getCode();

            switch (nuCode) {
                case Nucleotide::GAP:
                    siteType[nuIndex] |= GAP;
                    break;
                case Nucleotide::ADENINE:
                    siteType[nuIndex] |= A;
                    break;
                case Nucleotide::CYTOSINE:
                    siteType[nuIndex] |= C;
                    break;
                case Nucleotide::GUANINE:
                    siteType[nuIndex] |= G;
                    break;
                case Nucleotide::THYMINE:
                    siteType[nuIndex] |= T;
                    break;
            }
        }
    }

    int numAlleleTypes;
    for (unsigned long nuIndex = 0; nuIndex < extractedLength; nuIndex++) {
        if ((siteType[nuIndex] & GAP) != 0) nNoGap--;

        numAlleleTypes = 0;

        if ((siteType[nuIndex] & A) != 0) numAlleleTypes++;
        if ((siteType[nuIndex] & C) != 0) numAlleleTypes++;
        if ((siteType[nuIndex] & G) != 0) numAlleleTypes++;
        if ((siteType[nuIndex] & T) != 0) numAlleleTypes++;

        if (numAlleleTypes >= 2) {
            polymorphicMarkers[nuIndex] = true;
        } else {
            polymorphicMarkers[nuIndex] = false;
        }

        switch (numAlleleTypes) {
            case 0:
                nGap++;
                break;
            case 1:
                nMonoallelic++;
                break;
            case 2:
                nBiallelic++;
                break;
            case 3:
                nTriallelic++;
                break;
            case 4:
                nTetrallelic++;
                break;
        }
    }
}

const GenomeSet::DistanceStats GenomeSet::getDistanceStats(void) const {
    DistanceStats distanceStats;
    double totalPairsNum = static_cast<double> (getSize()) * static_cast<double> (getSize() - 1);

    distanceStats.min = numeric_limits<unsigned long>::max();
    distanceStats.max = 0;
    distanceStats.mean = 0.0;
    distanceStats.total = 0.0;

    for (unsigned long i = 0; i < getSize() - 1; i++) {
        Genome* modelGenome = getGenome(i);

        for (unsigned long k = i + 1; k < getSize(); k++) {
            unsigned long newDistance = modelGenome->ntDistanceTo(getGenome(k), true);

            distanceStats.total += static_cast<double> (newDistance);
            if (distanceStats.min > newDistance) {
                distanceStats.min = newDistance;
            }
            if (distanceStats.max < newDistance) {
                distanceStats.max = newDistance;
            }
        }
    }
    distanceStats.mean = distanceStats.total / totalPairsNum;

    return distanceStats;
}

void GenomeSet::removeGenomes(const vector<string>& blackList) {
    unsigned long newSize = 0;

    for (unsigned long genomeIndex = 0; genomeIndex < genomes.size(); genomeIndex++) {
        bool isInBlackList = false;

        for (unsigned long i = 0; i < blackList.size(); i++) {
            if (genomes[genomeIndex]->hasName(blackList[i])) {
                isInBlackList = true;
                break;
            }
        }

        if (isInBlackList) {
            delete genomes[genomeIndex];
            genomes[genomeIndex] = NULL;
        } else {
            genomes[newSize] = genomes[genomeIndex];
            newSize++;
        }
    }

    assert(newSize <= genomes.size());
    genomes.resize(newSize, NULL);
}

void GenomeSet::keepGenomes(const vector<string>& whiteList) {
    unsigned long newSize = 0;

    for (unsigned long genomeIndex = 0; genomeIndex < genomes.size(); genomeIndex++) {
        bool isInWhiteList = false;

        for (unsigned long i = 0; i < whiteList.size(); i++) {
            if (genomes[genomeIndex]->hasName(whiteList[i])) {
                isInWhiteList = true;
                break;
            }
        }

        if (isInWhiteList) {
            genomes[newSize] = genomes[genomeIndex];
            newSize++;
        } else {
            delete genomes[genomeIndex];
            genomes[genomeIndex] = NULL;
        }
    }

    assert(newSize <= genomes.size());
    genomes.resize(newSize, NULL);
}

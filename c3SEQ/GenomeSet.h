/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   GenomeSet.cpp
 * CREATED ON:  17 June 2011, 11:15
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: 
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-06-17  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   17 June 2011
 */

#ifndef GENOMESET_H
#define	GENOMESET_H

#include <cstdlib>
#include <ctime>

#include <vector>
#include <string>
#include <limits>
#include "Interface.h"
#include "Genome.h"
#include "Nucleotide.h"

/**
 * @note For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class GenomeSet {
public:

    explicit GenomeSet(void) {
        nGap = 0;
        nMonoallelic = 0;
        nBiallelic = 0;
        nTriallelic = 0;
        nTetrallelic = 0;
        nNoGap = 0;
    }

    GenomeSet(const GenomeSet& orig) : extractTemplate(orig.extractTemplate) {
        polymorphicMarkers = orig.polymorphicMarkers;
        genomes = orig.genomes; // shallow copying is enough.

        nGap = orig.nGap;
        nMonoallelic = orig.nMonoallelic;
        nBiallelic = orig.nBiallelic;
        nTriallelic = orig.nTriallelic;
        nTetrallelic = orig.nTetrallelic;
        nNoGap = orig.nNoGap;
    }

    ~GenomeSet() {
        for (unsigned i = 0; i < genomes.size(); i++) {
            if (genomes[i]) {
                delete genomes[i];
            }
        }
    }

    struct DistanceStats {
        unsigned long min;
        unsigned long max;
        double mean;
        double total;
    };

    Genome::ExtractTemplate* getExtractTemplate(void) {
        return &extractTemplate;
    };

    void initExtractTemplate(void) {
        assert(getSize() > 0);
        extractTemplate.initialize(getGenome(0));
    };

    /**
     * Get the number of active sequence in this dataset.
     * @return The number of active sequence in this dataset.
     */
    const unsigned long getNumActiveSeq(void) const {
        unsigned long numActive = 0;
        for (unsigned long i = 0; i < genomes.size(); i++) {
            if (genomes[i]->isActive()) {
                numActive++;
            }
        }
        return numActive;
    }

    /**
     * Activate all sequences in this dataset.
     */
    void activateAllGenomes(void) {
        for (unsigned long i = 0; i < genomes.size(); i++) {
            genomes[i]->setActive(true);
        }
    }

    /**
     * Randomly activate a number of sequences.
     * @param numActive The number of sequences that need to be activated.
     */
    void randomlyActivate(unsigned long numActive) {
        activateAllGenomes();
        if (numActive < genomes.size()) {
            unsigned long inactiveSize = genomes.size();
            unsigned long inactive[inactiveSize];
            for (unsigned long i = 0; i < inactiveSize; i++) {
                inactive[i] = i;
            }

            srand(time(NULL));
            for (unsigned long i = 0; i < numActive; i++) {
                unsigned long removed = rand() % inactiveSize;
                inactiveSize--;
                for (unsigned long k = removed; k < inactiveSize; k++) {
                    inactive[k] = inactive[k + 1];
                }
            }
            
            /* Deactivate all sequences in the inactive list */
            for (unsigned long i = 0; i < inactiveSize; i++) {
                unsigned long seqIndex = inactive[i];
                genomes[seqIndex]->setActive(false);
            }
        }
    }

    /**
     * Get the index of the given genome.
     * @param genome
     * @return  The index of the given genome.
     * @note    Make sure that the genome is contained by this dataset before
     *          calling this method; otherwise, an assertion false will occur.
     */
    const unsigned long getIndex(Genome* genomePtr) {
        for (unsigned long i = 0; i < genomes.size(); i++) {
            if (genomes[i] == genomePtr) {
                return i;
            }
        }
        assert(false); // should never reach here
    }

    /**
     * Access to the genome at the given index.
     * @param index
     * @return The pointer to the genome.
     */
    Genome * const getGenome(unsigned long index) const {
        assert(index >= 0 && index < genomes.size());
        return genomes[index];
    };

    /**
     * Search a genome by its accession. If no genome is found, the result 
     * will be a <b>NULL</b> pointer. If there are more than 1 genome that 
     * could be found, only return the 1st occurrence.
     * @param accession The accession number of the genome that must be found.
     * @return  The pointer to the genome. If no genome is found, the result 
     *          will be a <b>NULL</b> pointer.
     */
    Genome * const getGenome(string accession) const {
        for (unsigned long i = 0; i < getSize(); i++) {
            Genome* genome = getGenome(i);
            if (genome->hasName(accession)) {
                return genome;
            }
        }
        return NULL;
    };

    /**
     * Get the number of genomes contained in this dataset.
     * @return  The number of genomes contained in this dataset.
     */
    unsigned long const getSize(void) const {
        return genomes.size();
    }

    unsigned long const getGapNum(void) const {
        return nGap;
    };

    unsigned long const getMonoallelicNum(void) const {
        return nMonoallelic;
    };

    unsigned long const getBiallelicNum(void) const {
        return nBiallelic;
    };

    unsigned long const getTriallelicNum(void) const {
        return nTriallelic;
    };

    unsigned long const getTetrallelicNum(void) const {
        return nTetrallelic;
    };

    unsigned long const getNoGapSitesNum(void) const {
        return nNoGap;
    };

    /**
     * Add a genome into this genome set.
     * @param genomePtr
     */
    void addGenome(Genome* genomePtr) {
        assert(genomePtr);
        genomes.push_back(genomePtr);
    };

    /**
     * Remove the genome at the given index and return the removed genome.
     * @param index
     * @return  The pointer to the genome that has just been removed.
     */
    Genome* removeGenome(unsigned long index) {
        assert(index >= 0 && index < genomes.size());
        Genome* result = genomes[index];

        for (unsigned long i = index; i < genomes.size() - 1; i++) {
            genomes[i] = genomes[i + 1];
        }
        genomes.resize(genomes.size() - 1);

        return result;
    }

    /**
     * Remove all genomes in this set.
     * @note    All the genomes in this set will not be deleted by this method.
     */
    void clear(void) {
        genomes.clear();
        updateStats();
    }

    /**
     * Get the original length of the genomes (the length before cutting, 
     * splicing, etc) in this set.
     * @return  The original length of the genomes in this set.
     * @note    We assume that all the genomes in this set have the same original 
     *          length, so, the original length of the first genome is considered 
     *          as the common original length of all genomes.
     */
    unsigned long const getOriginalLength(void) const {
        if (genomes.size() <= 0) return 0;
        return genomes[0]->getFullLength();
    };

    /**
     * Get the length of the extracted genomes (length of the parts of the 
     * genomes that are used in the analysis) in this set.
     * @return  The extracted length of the genomes in this set.
     * @note    We assume that the extracted parts of all the genomes in this set 
     *          have the same length, so, the extracted length of the first genome 
     *          is considered as the common extracted length of all genomes.
     */
    unsigned long const getExtractedLength(void) const {
        if (genomes.size() <= 0) return 0;
        return genomes[0]->getExtractedLength();
    };

    /**
     * Count the number of distinct genomes in this data set. If the 
     * <b>autoRemoveDuplicate</b> flag is turned on, whenever 2 identical 
     * genomes are detected, one of them will be automatically removed,
     * so that the dataset contains distinct genomes only.
     * @param isGapsIgnored         Indicate if we should ignore gap sites when 
     *                              calculating the distances among genomes.
     * @param autoRemoveDuplicate   Turn this on when you want to remove 
     *                              identical genomes from the dataset.
     * @return  The number of distinct genomes in this dataset.
     * @note    Only the extracted parts of the genomes are taken into account.
     */
    unsigned long const getDistinctGenomes(const bool autoRemoveDuplicate);

    /**
     * Count the number of non-neighbour (have nucleotide distance greater than
     * the <b>maxNeighborDistance</b>) genomes in this data set. If the 
     * <b>autoRemoveNeighbor</b> flag is turned on, whenever 2 neighbour 
     * genomes are detected, the 2nd in the dataset will be automatically 
     * removed, so that the dataset contains non-neighbour genomes only.
     * @param maxNeighborDistance   If the distance between 2 genomes is 
     *                              equal or smaller than this value, these 2 
     *                              genomes will be consider as neighbour 
     *                              genomes.
     * @param autoRemoveNeighbor    Turn this on when you want to remove 
     *                              neighbour genomes from the dataset.
     * @return  The number of non-neighbour genomes in this dataset.
     * @note    Only the extracted parts of the genomes are taken into account.
     */
    unsigned long const getNonNeighborGenomes(
            const unsigned long maxNeighborDistance,
            const bool autoRemoveNeighbor);

    /**
     * Re-test to see if the sites in extracted genomes is polymorphic or not 
     * (the array <b>isPolymorphicSite</b> will be updated). Then, recount the 
     * number of gapped, monoallelic, biallelic, triallelic, and tetrallelic sites.<br>
     */
    void updateStats(void);

    /**
     * Get the statistics on the distances among all strain in this dataset.
     * @return  The statistic values.
     * @note    The distances are calculated base on the extracted parts of 
     *          the strains.
     */
    const DistanceStats getDistanceStats(void) const;

    /**
     * Apply the <b>extractTemplate</b> on all the genomes in this set to get 
     * the new extracted parts.
     */
    void applyTemplateForAllGenomes(void) {
        for (unsigned long i = 0; i < genomes.size(); i++) {
            genomes[i]->extractByTemplate(extractTemplate);
        }
    };

    /**
     * Remove all monomorphic sites from the <b>extractTemplate</b> vector.
     * @note    Make sure that <b>updateStats()</b> was called before executing 
     *          this method because this method makes use of 
     *          <b>isPolymorphicSite</b> vector, which is updated only by 
     *          <b>updateStats()</b> method.
     */
    void removeNonPolymorphicFromTemplate(void) {
        extractTemplate.keepPositions(polymorphicMarkers);
    };

    /**
     * Remove all genomes that have accession listed in the <b>blackList</b>
     * @param blackList The list which contains the accession numbers of 
     *                  the genomes that must be removed.
     * @note Accessions are case-insensitive.
     */
    void removeGenomes(const vector<string>& blackList);

    /**
     * Remove all genomes BUT the ones that is in the given <b>whiteList</b>
     * @param whiteList The list which contains the accession numbers of 
     *                  the genomes that must be kept. Any genome, 
     *                  of which the accession number is not in this list 
     *                  will be removed.
     * @note Accessions are case-insensitive.
     */
    void keepGenomes(const vector<string>& whiteList);


private:

    /* Disallow assignment. */
    GenomeSet& operator=(const GenomeSet& rhs);

    /** 
     * The vector that holds boolean values showing which nucleotides in the 
     * extracted genomes are polymorphic and which are not.
     */
    vector<bool> polymorphicMarkers;

    /**
     * Vector that store pointers to all genomes.
     */
    vector<Genome*> genomes;

    /** Number of gapped sites */
    unsigned long nGap;

    /* Number of monoallelic sites */
    unsigned long nMonoallelic;

    /* Number of biallelic sites */
    unsigned long nBiallelic;

    /* Number of triallelic sites */
    unsigned long nTriallelic;

    /* Number of tetrallelic sites */
    unsigned long nTetrallelic;

    /* Number of sites that have no gap */
    unsigned long nNoGap;

    /**
     * This template is a vector used to store all the position (segment index & 
     * base-pair position in segment) of the nucleotides that will be analyse.<br>
     * Normally, after initialising, this vector contains all the nucleotide 
     * positions. Then, it will be splice, cut out, etc due to the user requests. 
     * Finally, we go through all the genomes in this set and take out the 
     * nucleotides which have position stored in this vector to put them into 
     * the extracted parts (parts that will be analysed) of the genomes.<br>
     * The reason to use this template is to reduce the cost-of-shifting in 
     * vectors (instead of modifying many genomes many times, we just maintain 
     * this template vector and then apply it on all the genomes).
     * @note    Make sure that the method <b>initialiseExtractTemplate()</b> 
     *          is called before using this vector.
     */
    Genome::ExtractTemplate extractTemplate;

};

#endif	/* GENOMESET_H */


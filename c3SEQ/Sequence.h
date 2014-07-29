/* 
 * File:   Sequence.h
 * Author: Maciej F. Boni, Ha Minh Lam
 *
 * Created on 27 July 2011, 16:46
 */

#ifndef SEQUENCE_H
#define	SEQUENCE_H

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include "Nucleotide.h"
#include "Interface.h"

using namespace std;

class Sequence {
public:

    virtual ~Sequence() {
        for (unsigned long i = 0; i < nucleotides.size(); i++) {
            if (nucleotides[i]) {
                delete nucleotides[i];
            }
        }
    }

    virtual const string getName(void) const {
        return name;
    }

    /**
     * Get the number of nucleotides contained in this sequence.
     * @return  The number of nucleotides contained in this sequence.
     */
    virtual const unsigned long getLength(void) const {
        return nucleotides.size();
    }

    /**
     * Indicates if the name of this sequence is the same as the given name.
     * @param anotherName
     * @return  <b>true</b> if the name of this sequence is the same as the 
     *          given name.
     * @note    The comparison between names is case insensitive.
     */
    virtual const bool hasName(string anotherName) const {
        if (name.length() != anotherName.length()) {
            return false;
        }
        for (unsigned long i = 0; i < name.length(); i++) {
            if (toupper(name[i]) != toupper(anotherName[i])) {
                return false;
            }
        }

        return true;
    }

    /**
     * Get the nucleotide at the given index.
     * @param index The index (in the vector) of the nucleotide. This might 
     *              not the real position of the nucleotide.
     * @return The nucleotide at the given index.
     */
    virtual Nucleotide * const getNuAt(unsigned long index) const {
        assert(index >= 0 && index < nucleotides.size());
        return nucleotides[index];
    }

    virtual const string toString(void) const;

    /**
     * Count the nucleotide distance between this sequence and the given sequence.
     * @param another   The sequence that you want to find the distance to.
     * @return  The nucleotide distance.
     * @note    Only non-gap sites are taken into account.
     */
    virtual const unsigned long ntDistanceTo(Sequence* another, bool ignoreGap) const;

    /**
     * Check to see if this sequence is identical with the given sequence.
     * @param another
     * @return  <b>true</b> if this and the given sequence are identical; 
     *          <b>false</b> otherwise.
     * @note    Only non-gap sites are taken into account.
     */
    virtual const bool isIdenticalTo(Sequence* another, bool ignoreGap) const;


protected:

    /**
     * Create the sequence with the given name and DNA string.
     * @param dnaStr
     * @param newName
     * @note    This constructor is protected to prevent of using this class 
     *          directly. Its sub classes must be used instead.
     */
    explicit Sequence(const string& dnaStr, const string& newName)
    : name(newName) {

        try {
            nucleotides.reserve(dnaStr.length());

            for (unsigned long i = 0; i < dnaStr.length(); i++) {
                Nucleotide * nucleotidePtr = new Nucleotide(dnaStr[i], i, this);
                nucleotides.push_back(nucleotidePtr);
            }

        } catch (bad_alloc) {
            Interface::instance() << "Not enough memory to create new nucleotide object.\n";
            Interface::instance().showError(true, true);
        } catch (length_error) {
            Interface::instance() << "Not enough memory to store the sequence " << name << "\n";
            Interface::instance().showError(true, true);
        }
    }

    /**
     * Create an empty sequence with the given name.
     * @param newName
     * @note    This constructor is protected to prevent of using this class 
     *          directly. Its sub classes must be used instead.
     */
    explicit Sequence(const string& newName)
    : name(newName) {
    }


    /** Sequence name */
    const string name;

    /** The container of all the nucleotides in this sequence. */
    vector<Nucleotide*> nucleotides;


private:
    Sequence(const Sequence& orig);

    Sequence& operator=(const Sequence& rhs);

};


#endif	/* SEQUENCE_H */


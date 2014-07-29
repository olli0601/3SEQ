/**
 * FILE NAME:   SingleBpRun.cpp
 * CREATED ON:  31 August 2011, 16:14
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "SingleBpRun.h"
 * 
 * HISTORY:     Version     Date            Description
 *              1.0         2011-08-31      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   31 August 2011
 */

#include "SingleBpRun.h"
#include "Triplet.h"

SingleBpRun::SingleBpRun(int argc, char** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Single Break-Point Run");
}

SingleBpRun::SingleBpRun(const SingleBpRun& orig) {
    assert(false); // should never reach here
}

SingleBpRun& SingleBpRun::operator=(const SingleBpRun& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach here
    }
    return *this;
}

SingleBpRun::~SingleBpRun() {
}

void SingleBpRun::parseCmdLine() {
    int argNum = getRunArgsNum();

    if (argNum < 1) {
        Interface::instance() << "Not enough parameter for single-break-point run.\n";
        Interface::instance().showError(true, true);
    }

    string filePath = argVector[2];
    SequenceFile parentFile(filePath, SequenceFile::UNKNOWN);
    if (!parentFile.exists()) {
        Interface::instance() << "File \"" << filePath << "\" not found.\n";
        Interface::instance().showError(true, true);
    }
    parentFile.detectFileType();
    parentDataset = new GenomeSet();
    parentFile.read(parentDataset, "");

    if (argNum >= 2) {
        argNum = 2;
        filePath = argVector[3];
        SequenceFile childFile(filePath, SequenceFile::UNKNOWN);
        if (!childFile.exists()) {
            Interface::instance() << "File \"" << filePath << "\" not found.\n";
            Interface::instance().showError(true, true);
        }
        childFile.detectFileType();
        childDataset = new GenomeSet();
        childFile.read(childDataset, "");
    }

    deleteCmdLineArgs(argNum + 1);
    Run::parseCmdLine(); // Get common run options
}

void SingleBpRun::perform() {
    Run::perform(); // process common data
    prepare();

    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);

    progress();

    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);

    displayResult();
}

void SingleBpRun::prepare(void) {
    /* Calculate the number of actual comparisons */
    double totalTripletNum = getNumTotalTriplets();

    if (totalTripletNum <= 0) {
        Interface::instance()
                << "The number of sequences is not enough (at least 3 sequences needed).\n";
        Interface::instance().showError(true, true);
    } else if (totalTripletNum > MAX_UL) {
        Interface::instance()
                << "The number of comparisons you are requesting is larger than LONG_MAX for an\n"
                << "\"unsigned long\" type on your system. The P-values below should be correct,\n"
                << "but the numbers of comparisons that were approximated/computed/skipped might\n"
                << "not be the exact values.";
        Interface::instance().showWarning();
    }

    /* Show some logs */
    Interface::instance()
            << "Using " << parentDataset->getSize() << " sequences as parents.\n"
            << "Using " << childDataset->getSize() << " sequences as children.\n";
    Interface::instance().showLog(true);

    Interface::instance() << "Need a p-value of "
            << cloRejectThreshold / getNumTripletsForCorrection()
            << " to survive multiple comparisons correction.\n";
    Interface::instance().showLog(true);
}

void SingleBpRun::progress(void) {
    bool isStopped = false;
    char minPValStr[20];
    minPVal = 1.0;

    Interface::instance() << "Progress    Time Elapsed    Minimum P-Value" << endl;
    Interface::instance().showLog(false);

    /* Initialise progress counter */
    double totalLoop = static_cast<double> (childDataset->getSize())
            * static_cast<double> (parentDataset->getSize());
    Interface::instance().initCounter("", 0, totalLoop);

    /* Progressing begins */
    for (unsigned long childIndex = 0; !isStopped && childIndex < childDataset->getSize(); childIndex++) {
        Genome* child = childDataset->getGenome(childIndex);

        for (unsigned long dadIndex = 0; !isStopped && dadIndex < parentDataset->getSize(); dadIndex++) {
            Genome* dad = parentDataset->getGenome(dadIndex);
            if (dad == child) continue;

            /* Show progress counter */
            double currentLoop = static_cast<double> (childIndex)
                    * static_cast<double> (parentDataset->getSize())
                    + static_cast<double> (dadIndex);
            sprintf(minPValStr, "%e", minPVal);
            Interface::instance()
                    << "     " << Interface::instance().getElapsedTime().c_str()
                    << "       " << minPValStr;
            Interface::instance().count(currentLoop, false);
            /* End progress counter */

            for (unsigned long mumIndex = 0; !isStopped && mumIndex < parentDataset->getSize(); mumIndex++) {
                Genome* mum = parentDataset->getGenome(mumIndex);
                if (mum == child || mum == dad) continue;

                numComputed += 1.0;
                Triplet triplet(dad, mum, child);
                double pVal = triplet.getPValSingleBp();

                if (pVal < minPVal) {
                    minPVal = pVal;
                    if (stats::correction::dunnSidak(minPVal, getNumTripletsForCorrection()) < cloRejectThreshold
                            && cloStopAtFirstRecombinant) {
                        /* Finish progressing early */
                        isStopped = true;
                    }
                }                
            }
        }
    }

    /* Progressing finished */
    sprintf(minPValStr, "%e", minPVal);
    Interface::instance()
            << "     " << Interface::instance().getElapsedTime().c_str()
            << "       " << minPValStr;
    Interface::instance().finishCounting();
}

void SingleBpRun::displayResult(void) {
    double numTotalTriplet = getNumTotalTriplets();
    
    Interface::instance()
            << "Number of triples tested :         " << numTotalTriplet << "\n"
            << "Number of p-values computed :      " << numComputed << "\n"
            << "Number of p-values not computed :  " << numTotalTriplet - numComputed << "\n";
    Interface::instance().showOutput(true);

    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);

    char formatedPVal[20];
    sprintf(formatedPVal,
            "%1.3e",
            static_cast<double> (stats::correction::dunnSidak(minPVal, getNumTripletsForCorrection()))
            );
    Interface::instance()
            << "Rejection of the null hypothesis of clonal evolution at p = "
            << stats::correction::dunnSidak(minPVal, getNumTripletsForCorrection()) << "\n"
            << "                                                        p = "
            << formatedPVal << "\n"
            << "                                            Uncorrected p = "
            << minPVal << "\n"
            << "                                            Bonferroni  p = "
            << stats::correction::bonferroni(minPVal, getNumTripletsForCorrection()) << "\n";
    Interface::instance().showOutput(true);    
}
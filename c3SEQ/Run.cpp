/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   Run.cpp
 * CREATED ON:  12 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "Run.h"
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-12  Created
 *
 * VERSION:     1.0
 * LAST EDIT:   12 August 2011
 */

#include "Run.h"

#include <cmath>
#include <cstdlib>

#include "Triplet.h"

#include "HelpRun.h"
#include "VersionRun.h"
#include "InfoRun.h"
#include "WriteRun.h"
#include "FullRun.h"
#include "GenPRun.h"
#include "CheckPRun.h"
#include "MatchRun.h"
#include "TripletRun.h"
#include "SingleBpRun.h"


////////////////////////////////////////////////////////////////////////////////
// Define static constants 
////////////////////////////////////////////////////////////////////////////////
const unsigned long Run::NaUL = numeric_limits<unsigned long>::max();
const unsigned long Run::MAX_UL = numeric_limits<unsigned long>::max() - 1;

////////////////////////////////////////////////////////////////////////////////

Run::Run() {
    /* This type of constructor should never be used */
    assert(false);
}

Run::Run(int argc, char** argv) : argCount(argc), argVector(argv) {
    Interface::instance().storeCommand(argc, argv);


    /* PRIVATE */
    isChildParentInDiffFiles = false;

    cloThirdPositionsOnly = false;
    cloFirstAndSecondPositionsOnly = false;

    cloBeginSequence = NaUL;
    cloEndSequence = NaUL;

    numTotalTriplets = -1.0; // indicates that it has not been calculated
    numTripletsForCorrection = -1.0; // indicates that it has not been calculated

    for (int i = 0; i < PVAL_HISTOGRAM_SIZE; i++) {
        pValsHistogram[i] = 0;
    }

    /* PROTECTED */
    id = "";

    cloIsSkippedTriplesRecorded = false;

    fileSubset = NULL;

    parentDataset = NULL;
    childDataset = NULL;

    cloUseAllSites = false;

    cloStopAtFirstRecombinant = false;
    cloIsStartedByPythonScript = false;

    cloRejectThreshold = DEFAULT_REJECT_THRESHOLD;

    cloFirstNucleotidePos = NaUL;
    cloLastNucleotidePos = NaUL;

    cloIsFirstToLastCutOut = false;

    cloIsAllTriplesUsed = true;

    cloIsSubsetRemoved = false;
    cloIsIdenticalSeqRemoved = false;
    cloRemoveDistance = 0;

    cloAllBpCalculated = false;

    cloNoBpCalculated = false;

    cloNo3sRecFile = false;

    cloUseSiegmundApprox = false;

    cloMinLongRecombinantLength = DEFAULT_MIN_RECOMBINANT_LENGTH;

    cloOutputFileType = SequenceFile::PHYLIP;
}

Run::Run(const Run& orig) : argCount(orig.argCount), argVector(orig.argVector) {
    assert(false); // should never reach to this line
}

Run& Run::operator=(const Run& rhs) {
    if (this != &rhs) {
        assert(false); // should never reach to this line
    }
    return *this;
}

Run::~Run() {
    if (childDataset != NULL) {
        if (!isChildParentInDiffFiles) {
            /* This means parentDataset and childDataset share pointers */
            childDataset->clear();
        }
        delete childDataset;
    }

    if (parentDataset != NULL) delete parentDataset;

    if (fileSubset != NULL) delete fileSubset;

}

Run* Run::getRun(int argc, char** argv) {
    if (argc < 2) {
        return NULL;
    }

    string mode = argv[1];

    if (mode == "-help" || mode == "-h") {
        return new HelpRun(argc, argv);

    } else if (mode == "-version" || mode == "-v") {
        return new VersionRun(argc, argv);

    } else if (mode == "-gen-p" || mode == "-g") {
        return new GenPRun(argc, argv);

    } else if (mode == "-check-p" || mode == "-c") {
        return new CheckPRun(argc, argv);

    } else if (mode == "-info" || mode == "-i") {
        return new InfoRun(argc, argv);

    } else if (mode == "-write" || mode == "-w") {
        return new WriteRun(argc, argv);

    } else if (mode == "-full" || mode == "-f") {
        return new FullRun(argc, argv);

    } else if (mode == "-match" || mode == "-m") {
        return new MatchRun(argc, argv);

    } else if (mode == "-triplet" || mode == "-t") {
        return new TripletRun(argc, argv);

    } else if (mode == "-single" || mode == "-s") {
        return new SingleBpRun(argc, argv);

    }

    return NULL;
}

void Run::deleteCmdLineArgs(int numDelete) {
    if (numDelete >= argCount) {
        numDelete = argCount - 1;
    }

    assert(numDelete >= 1);

    argCount -= numDelete;

    for (int i = 1; i < argCount; i++) {
        argVector[i] = argVector[i + numDelete];
    }
}

const int Run::getRunArgsNum(void) const {
    if (argCount < 2) return -1;

    int result = 0;
    for (int i = 2; i < argCount; i++) {
        if (argVector[i][0] != '-') {
            result++;
        } else {
            return result;
        }
    }
    return result;
}

void Run::generateAutoID(void) {
    time_t currentTime;
    struct tm *timeInfo;
    stringstream tmpStream;

    time(&currentTime);
    timeInfo = localtime(&currentTime);

    tmpStream
            << String::formatInt(timeInfo->tm_year % 100, 2)
            << String::formatInt(timeInfo->tm_mon + 1, 2)
            << String::formatInt(timeInfo->tm_mday, 2)
            << "_"
            << String::formatInt(timeInfo->tm_hour, 2) << "h"
            << String::formatInt(timeInfo->tm_min, 2) << "m"
            << String::formatInt(timeInfo->tm_sec, 2) << "s";

    /* Add millisecond into ID to make sure of getting rid of identicalness.*/
    struct timeval timeInSec;
    gettimeofday(&timeInSec, NULL);
    int milisec = timeInSec.tv_usec / 1000;
    tmpStream << "" << String::formatInt(milisec, 2);

    id = tmpStream.str();
}

void Run::parseCmdLine() {
    try {
        for (int i = 1; i < argCount; i++) {
            if (argVector[i][0] == (char) 0) continue;
            
            if (argVector[i][0] != '-') {
                throw (argVector[i]);
            }


            string option = argVector[i];

            switch (argVector[i][1]) {
                case '1':
                    if (option == "-12po") {
                        cloFirstAndSecondPositionsOnly = true;
                        cloThirdPositionsOnly = false;
                    } else {
                        throw (argVector[i]);
                    }
                    break;

                case '3':
                    if (option == "-3po") {
                        cloFirstAndSecondPositionsOnly = false;
                        cloThirdPositionsOnly = true;
                    } else {
                        throw (argVector[i]);
                    }
                    break;

                case 'a':
                    cloUseAllSites = true;
                    break;

                case 'b':
                    if (option == "-bp-all") {
                        cloAllBpCalculated = true;
                        cloNoBpCalculated = false;
                    } else if (option == "-bp-none") {
                        cloAllBpCalculated = false;
                        cloNoBpCalculated = true;
                    } else if (option == "-b") {
                        Interface::instance()
                                << "With -b option, you must give a starting position in the sequence array,\n"
                                << "e.g. -b0 for the first sequence or -b2 for the third sequence.\n";
                        Interface::instance().showError(true, true);
                    } else {
                        cloBeginSequence = strtoul(argVector[i] + 2, NULL, 0);
                    }
                    break;

                case 'd':
                    cloIsIdenticalSeqRemoved = true;
                    if (argVector[i] + 2 != '\0') {
                        cloRemoveDistance = strtoul(argVector[i] + 2, NULL, 0);
                    }
                    break;

                case 'e':
                    if (option == "-e") {
                        Interface::instance()
                                << "With -e option, you must give an ending position in the sequence array,\n"
                                << "e.g. -e0 for the first sequence or -e2 for the third sequence.\n";
                        Interface::instance().showError(true, true);
                    } else {
                        cloEndSequence = strtoul(argVector[i] + 2, NULL, 0) + 1;
                    }
                    break;

                case 'f':
                    if (option == "-fasta") {
                        cloOutputFileType = SequenceFile::FASTA;
                    } else {
                        cloFirstNucleotidePos = strtoul(argVector[i] + 2, NULL, 0) - 1;
                    }
                    break;

                case 'h':
                    if (option == "-hs") {
                        cloUseSiegmundApprox = true;
                    } else {
                        throw (argVector[i]);
                    }
                    break;

                case 'i':
                    if (option == "-id") {
                        if (i == argCount - 1) {
                            Interface::instance()
                                    << "No identifier string given after -id option.\n";
                            Interface::instance().showError(true, true);
                        }
                        if (argVector[i + 1][0] == '-') {
                            Interface::instance()
                                    << "Identifier string after -id option cannot start with \"-\" character.\n";
                            Interface::instance().showError(true, true);
                        }
                        i++;
                        if (id.length() <= 0) {
                            id = argVector[i];
                        }
                    } else if (option == "-id-auto") {
                        if (id.length() <= 0) {
                            generateAutoID();
                        }
                    } else {
                        throw (argVector[i]);
                    }
                    break;

                case 'l':
                    cloLastNucleotidePos = strtoul(argVector[i] + 2, NULL, 0);
                    break;

                case 'L':
                    cloMinLongRecombinantLength = strtoul(argVector[i] + 2, NULL, 0);
                    break;

                case 'n': option = argVector[i];
                    if (option == "-nexus") {
                        cloOutputFileType = SequenceFile::NEXUS;
                    } else if (option == "-nr") {
                        cloNo3sRecFile = true;
                    } else {
                        throw (argVector[i]);
                    }
                    break;

                case 'p':
                    cloIsStartedByPythonScript = true;
                    break;

                case 'r':
                    cloIsSkippedTriplesRecorded = true;
                    break;

                case 's': option = argVector[i];
                    if (i == argCount - 1) {
                        Interface::instance() << "No filename given after subset option.\n";
                        Interface::instance().showError(true, true);
                    }
                    if (option == "-subset") {
                        i++;
                        fileSubset = new TextFile(argVector[i]);
                    } else if (option == "-subset-remove") {
                        i++;
                        fileSubset = new TextFile(argVector[i]);
                        cloIsSubsetRemoved = true;
                    } else {
                        throw (argVector[i]);
                    }
                    if (!fileSubset->exists()) {
                        Interface::instance() << "File \"" << fileSubset->getPath()
                                << "\" does not exist.\n";
                        Interface::instance().showError(true, true);
                    }
                    break;

                case 't': option = argVector[i];
                    if (option == "-t") {
                        Interface::instance()
                                << "With -t option, you must give a value for the rejections threshold,\n"
                                << "e.g. -t0.20 or -t1e-6\n";
                        Interface::instance().showError(true, true);
                    } else {
                        // here we assume -t was used properly, i.e. "-t0.01"
                        cloRejectThreshold = atof(argVector[i] + 2);
                    }
                    break;

                case 'x':
                    cloIsFirstToLastCutOut = true;
                    break;

                case 'y':
                    cloStopAtFirstRecombinant = true;
                    break;

                case '#':
                    cloIsAllTriplesUsed = false;
                    break;

                default:
                    throw (argVector[i]);
            }
        }

    } catch (char* unknownArg) {
        Interface::instance() << "Invalid option " << unknownArg << "\n";
        Interface::instance().showError(true, true);
    }
}

void Run::initLogFile(void) const {
    if (!isLogFileSupported()) return;

    string logFileName = DEFAULT_LOG_FILE_NAME;
    if (id.length() > 0) {
        logFileName = id + "." + logFileName;
    }
    Interface::instance().initLogFile(logFileName);
}

void Run::perform() {
    /* Prepare file names */
    pValHistogramFileName = DEFAULT_PVAL_HISTOGRAM_FILE_NAME;
    if (id.length() > 0) {
        pValHistogramFileName = id + "." + pValHistogramFileName;
    }

    /* Check the dataset */
    if (parentDataset == NULL || parentDataset->getSize() <= 0) {
        Interface::instance() << "No sequence has been read in.\n";
        Interface::instance().showError(true, true);
    }

    if (childDataset != NULL) {
        /* This means we have children and parents in separate files */
        isChildParentInDiffFiles = true;

        if (parentDataset->getOriginalLength() != childDataset->getOriginalLength()) {
            Interface::instance()
                    << "Parents sequences and children sequences must have the same length.\n"
                    << "Please, check your dataset.\n";
            Interface::instance().showError(true, true);
        }
        if (cloBeginSequence != NaUL || cloEndSequence != NaUL) {
            Interface::instance()
                    << "Parents and children sequences are in different files. Ignore -b and -e options.";
            Interface::instance().showWarning(true);
            cloBeginSequence = NaUL;
            cloEndSequence = NaUL;
        }
        if (cloIsAllTriplesUsed) {
            Interface::instance()
                    << "Auto-enable -# option; actual number of comparisons will be used in corrections.";
            Interface::instance().showWarning(true);
            cloIsAllTriplesUsed = false;
        }

    } else {
        isChildParentInDiffFiles = false;
    }

    /* Check if the user wants to analyse a particular subset of sequences */
    if (fileSubset != NULL) {
        applySubset();
    }

    /* Verify -f and -l options */
    verifyFirstAndLastNuPos();

    /* Cut, splice genome, strip out certain codon positions due to user's request */
    extractGenomes();

    /* Remove non-polymorphic sites if needed */
    removeNonPolymorphic();

    /* Remove neighbour and identical sequences if needed */
    removeNeighborAndIdenticalGenomes();

    if (!isChildParentInDiffFiles) {
        /* Verify -b and -e options */
        verifyBeginAndEndSeqIndex();

        if (cloBeginSequence > 0 || cloEndSequence < parentDataset->getSize()) {
            childDataset = new GenomeSet();
            for (unsigned long i = cloBeginSequence; i < cloEndSequence; i++) {
                childDataset->addGenome(parentDataset->getGenome(i));
            }
            childDataset->updateStats();

        } else {
            /* Just copy the whole parent set */
            childDataset = new GenomeSet(*parentDataset);
        }
    }
}

void Run::applySubset(void) {
    vector<string> accessionList = fileSubset->readAllLines();
    vector<GenomeSet*> datasetList; // the reason to use this vector is just to make the code more concise.
    datasetList.push_back(parentDataset);
    datasetList.push_back(childDataset);

    for (unsigned long i = 0; i < datasetList.size(); i++) {
        GenomeSet* dataset = datasetList[i];
        if (dataset == NULL || dataset->getSize() <= 0) continue;

        if (cloIsSubsetRemoved) {
            dataset->removeGenomes(accessionList);
        } else {
            dataset->keepGenomes(accessionList);
        }
    }

    if (cloIsIdenticalSeqRemoved == true
            || cloBeginSequence != NaUL || cloEndSequence != NaUL) {
        Interface::instance() << "-d, -b and -e options will be ignore when using subset option.";
        Interface::instance().showWarning(true);

        /* Ignore -d, -b and -e options */
        cloIsIdenticalSeqRemoved = false;
        cloBeginSequence = NaUL;
        cloEndSequence = NaUL;
    }
}

void Run::extractGenomes(void) {
    vector<GenomeSet*> datasetList; // the reason to use this vector is just to make the code more concise.
    datasetList.push_back(parentDataset);
    datasetList.push_back(childDataset);

    for (unsigned long i = 0; i < datasetList.size(); i++) {
        GenomeSet* dataset = datasetList[i];
        if (dataset == NULL || dataset->getSize() <= 0) continue;

        dataset->initExtractTemplate(); // initialise template vector for extraction

        /* Cut the template */
        if (cloIsFirstToLastCutOut) {
            dataset->getExtractTemplate()->cutOut(cloFirstNucleotidePos, cloLastNucleotidePos);
        } else {
            dataset->getExtractTemplate()->shrink(cloFirstNucleotidePos, cloLastNucleotidePos);
        }

        /* Stripping out certain codon positions */
        if (cloFirstAndSecondPositionsOnly) {
            dataset->getExtractTemplate()->deleteCodon3();
        } else if (cloThirdPositionsOnly) {
            dataset->getExtractTemplate()->deleteCodon1n2();
        }

        /* Actually, the sequences still have not been modified until here */
        dataset->applyTemplateForAllGenomes();
        dataset->updateStats();
    }
}

void Run::removeNonPolymorphic(void) {
    if (cloUseAllSites) return;
    assert(parentDataset);

    GenomeSet* tmpDataset = new GenomeSet(*parentDataset);

    if (isChildParentInDiffFiles && childDataset && childDataset->getSize() > 0) {
        /* Merge 2 dataset */
        for (unsigned long i = 0; i < childDataset->getSize(); i++) {
            tmpDataset->addGenome(childDataset->getGenome(i));
        }
        tmpDataset->updateStats();
    }

    /* Remove non-polymorphic sites from template */
    tmpDataset->removeNonPolymorphicFromTemplate();

    /* Modify all sequences */
    tmpDataset->applyTemplateForAllGenomes();

    if (parentDataset) {
        parentDataset->updateStats();
        *(parentDataset->getExtractTemplate()) = *(tmpDataset->getExtractTemplate());
    }
    if (childDataset) {
        childDataset->updateStats();
        *(childDataset->getExtractTemplate()) = *(tmpDataset->getExtractTemplate());
    }

    tmpDataset->clear();
    delete tmpDataset;
}

void Run::removeNeighborAndIdenticalGenomes(void) {
    if (!cloIsIdenticalSeqRemoved) return;

    vector<GenomeSet*> datasetList; // the reason to use this vector is just to make the code more concise.
    datasetList.push_back(parentDataset);
    datasetList.push_back(childDataset);

    if (cloRemoveDistance > 0) {
        Interface::instance() << "Removing neighbour sequences (nt distance <= " << cloRemoveDistance << ").\n";
    } else {
        Interface::instance() << "Removing identical sequences.\n";
    }
    Interface::instance().showLog(true);

    for (unsigned long i = 0; i < datasetList.size(); i++) {
        GenomeSet* dataset = datasetList[i];
        if (dataset == NULL || dataset->getSize() <= 0) continue;

        if (cloRemoveDistance > 0) {
            dataset->getNonNeighborGenomes(cloRemoveDistance, true);
        } else {
            dataset->getDistinctGenomes(true);
        }
        dataset->updateStats();
    }
}

void Run::verifyFirstAndLastNuPos(void) {
    assert(parentDataset && parentDataset->getSize() > 0);

    if (cloFirstNucleotidePos == NaUL) {
        cloFirstNucleotidePos = 0;
    }
    if (cloLastNucleotidePos == NaUL) {
        cloLastNucleotidePos = parentDataset->getOriginalLength();
    }

    /* Check to make sure -f and -l options are within bounds */
    if (cloFirstNucleotidePos < 0 || cloFirstNucleotidePos > parentDataset->getOriginalLength() - 1) {
        Interface::instance()
                << "Cannot use " << cloFirstNucleotidePos + 1
                << " as first position. Setting first position to 1.";
        Interface::instance().showWarning(true);
        cloFirstNucleotidePos = 0;
    }

    if (cloLastNucleotidePos < 1 || cloLastNucleotidePos > parentDataset->getOriginalLength()) {
        Interface::instance()
                << "Cannot use " << cloLastNucleotidePos
                << " as last position. Setting last position to "
                << parentDataset->getOriginalLength() << ".";
        Interface::instance().showWarning(true);
        cloLastNucleotidePos = parentDataset->getOriginalLength();
    }

    if (cloFirstNucleotidePos >= cloLastNucleotidePos) {
        Interface::instance()
                << "First position comes after last position. Ignoring -f, -l, and -x options.\n"
                << "Using whole genome.";
        Interface::instance().showWarning(true);
        cloFirstNucleotidePos = 0;
        cloLastNucleotidePos = parentDataset->getOriginalLength();
    }

    if (cloFirstNucleotidePos == 0 && cloLastNucleotidePos == parentDataset->getOriginalLength()) {
        cloIsFirstToLastCutOut = false;
    }

    if (!cloIsStartedByPythonScript && cloIsFirstToLastCutOut && cloThirdPositionsOnly) {
        Interface::instance()
                << "You are cutting out a middle segment of the sequence and requesting that" << endl
                << "only third positions be used. Make sure that your starting position is the" << endl
                << "first position of a codon and that the segment you are cutting out has a" << endl
                << "length that is divisible by three.";
        Interface::instance().showWarning(true);
    }
}

void Run::verifyBeginAndEndSeqIndex(void) {
    assert(parentDataset && parentDataset->getSize() > 0);

    if (cloBeginSequence == NaUL) {
        cloBeginSequence = 0;
    }
    if (cloEndSequence == NaUL) {
        cloEndSequence = parentDataset->getSize();
    }

    if (cloBeginSequence < 0 || cloBeginSequence >= parentDataset->getSize()) {
        Interface::instance()
                << "-b option: beginning sequence index " << cloBeginSequence
                << " must be between 0 and " << parentDataset->getSize() - 1 << ".\n"
                << "If you used the \"-d\" option to remove duplicate sequences, your\n"
                << "requested start point with the \"-b\" option may now be out of range.\n";
        Interface::instance().showError(true, true);
    }

    if (cloEndSequence == 0) {
        cloEndSequence = parentDataset->getSize();
    } else if (cloEndSequence < 0 || cloEndSequence > parentDataset->getSize()) {
        Interface::instance()
                << "-e option: end sequence index " << cloEndSequence - 1
                << " must be between 0 and " << parentDataset->getSize() - 1 << ".\n"
                << "If you used the \"-d\" option to remove duplicate sequences, your "
                << "requested end point with the \"-e\" option may now be out of range.\n";
        Interface::instance().showError(true, true);
    }

    if (cloBeginSequence >= cloEndSequence) {
        Interface::instance()
                << "-b & -e options: beginning sequence has index after end sequence.\n";
        Interface::instance().showError(true, true);
    }
}

void Run::calculateTripletNum(void) {
    /* Calculate the number of actual comparisons */
    if (parentDataset->getNumActiveSeq() < 2 || childDataset->getNumActiveSeq() < 1) {
        numTotalTriplets = 0.0;

    } else {
        numTotalTriplets = static_cast<double> (childDataset->getNumActiveSeq())
                * static_cast<double> (parentDataset->getNumActiveSeq())
                * static_cast<double> (parentDataset->getNumActiveSeq() - 1);
        if (!isChildParentInDiffFiles) {
            numTotalTriplets = static_cast<double> (childDataset->getNumActiveSeq())
                    * static_cast<double> (parentDataset->getNumActiveSeq() - 1)
                    * static_cast<double> (parentDataset->getNumActiveSeq() - 2);
        }
    }

    /* Calculate the number of comparisons used in Dunn-Sidak corrections */
    numTripletsForCorrection = numTotalTriplets;
    if (!isChildParentInDiffFiles && cloIsAllTriplesUsed) {
        numTripletsForCorrection = static_cast<double> (parentDataset->getNumActiveSeq())
                * static_cast<double> (parentDataset->getNumActiveSeq() - 1)
                * static_cast<double> (parentDataset->getNumActiveSeq() - 2);
    }
}

void Run::addPValIntoHistogram(double pValue) {
    int index = static_cast<int> (-1.0 * log10(pValue));
    if (index >= PVAL_HISTOGRAM_SIZE) {
        index = PVAL_HISTOGRAM_SIZE - 1;
    }
    if (index < 0) {
        index = 0;
    }
    pValsHistogram[index]++;
}

void Run::savePValHistogram(void) {
    TextFile pHistFile(pValHistogramFileName);

    if (pHistFile.exists()) {
        Interface::instance()
                << "The histogram of P-values will be recorded into the file \""
                << pValHistogramFileName << "\".\n"
                << "This file already exists. Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            Interface::instance() << "The histogram of P-values is not saved.\n";
            Interface::instance().showLog(true);
            return;

        } else {
            pHistFile.removeFile();
        }
    }

    try {
        pHistFile.openToWrite();
        double totalTripletsNum = getNumTotalTriplets();

        char bufStr[200];
        for (int i = 0; i < PVAL_HISTOGRAM_SIZE; i++) {
            double histVal = static_cast<double> (pValsHistogram[i]);
            double histPerTriplet = (pValsHistogram[i] == 0) ? 0.0 : histVal / totalTripletsNum;

            if (histPerTriplet != 0.0) {
                sprintf(bufStr,
                        "%d\t%20lu\t%1.8f\t%1.3f",
                        i, pValsHistogram[i], histPerTriplet, log10(histPerTriplet)
                        );
            } else {
                sprintf(bufStr,
                        "%d\t%20lu\t%1.8f\tN/A",
                        i, pValsHistogram[i], histPerTriplet
                        );
            }

            pHistFile.writeLine(bufStr);
        }

        pHistFile.close();

        Interface::instance()
                << "The histogram of P-values has been saved in the file \""
                << pValHistogramFileName << "\".\n";
        Interface::instance().showLog(true);

    } catch (...) {
        Interface::instance() << "An error occurs during saving progress. Cannot save histogram of P-values.\n";
        Interface::instance().showError(true, false); // show error but not exit
    }
}

void Run::loadPTable(PTableFile* pTableFile) {
    if (pTableFile != NULL) {
        
        PTableFile::ReadResult readResult = pTableFile->tryLoadInto(PTable::instance());        
        switch (readResult) {
            case PTableFile::INVALID_FILE:
                Interface::instance() << "Invalid P-value table file.\n";
                Interface::instance().showError(true, true);
                break;

            case PTableFile::WRONG_ARCH:
                Interface::instance() << "The P-value table file is not compatible with your system architecture.\n";
                Interface::instance().showError(true, true);
                break;

            case PTableFile::CANCELLED:
                Interface::instance() << "Loading is cancelled.\n";
                Interface::instance().showLog(true);
                Interface::instance().throwExitSignal(0);
                break;

            case PTableFile::FILE_CORRUPT:
                Interface::instance() << "Loading fail! The P-value table file may be corrupt.\n";
                Interface::instance().showError(true, true);
                break;
                
            case PTableFile::SUCCESS:
//                Interface::instance() << "Loading successful.\n";
//                Interface::instance().showLog(true);
                break;
        }

    } else {
        Interface::instance() << "Searching for P-value table file...\n";
        Interface::instance().showLog(true);

        PTableFile::ReadResult readResult = PTableFile::searchAndLoad(PTable::instance());
        switch (readResult) {
            case PTableFile::WRONG_ARCH:    // fall through
            case PTableFile::INVALID_FILE:
                Interface::instance() << "Cannot find any valid P-value table file.\n";
                Interface::instance().showError(true, true);
                break;

            case PTableFile::CANCELLED:
                Interface::instance() << "Loading is cancelled.\n";
                Interface::instance().showLog(true);
                Interface::instance().throwExitSignal(0);
                break;

            case PTableFile::FILE_CORRUPT:
                Interface::instance() << "Loading fail! The P-value table file may be corrupt.\n";
                Interface::instance().showError(true, true);
                break;
                
            case PTableFile::SUCCESS:
//                Interface::instance() << "Loading successful.\n";
//                Interface::instance().showLog(true);
                break;
        }
    }
}

/**
 * FILE NAME:   PTableFile.cpp
 * CREATED ON:  07 October 2011, 15:48
 * AUTHOR:      Ha Minh Lam
 *
 * DESCRIPTION: See "PTableFile.h"
 *
 * HISTORY:     Version     Date            Description
 *              1.0         2011-10-07      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   07 October 2011
 */

#include "PTableFile.h"
#include "Interface.h"
#include "Util.h"

#include <dirent.h>
#include <string>


////////////////////////////////////////////////////////////////////////////////
//  STATIC CONSTANTS
////////////////////////////////////////////////////////////////////////////////

const char PTableFile::FILE_MARKER[] = "P-table";

const string PTableFile::CONF_KEY = "pTableFile";

const string PTableFile::getConfigKey() {
    return CONF_KEY;
}

////////////////////////////////////////////////////////////////////////////////

const bool PTableFile::exists(void) {
    if (filePath.length() <= 0) {
        return false;
    }

    fstream stream;
    stream.open(filePath.c_str(), ios::in);

    if (stream.is_open() && !stream.fail() && !stream.eof()) {
        stream.close();
        return true;
    } else {
        return false;
    }
}

const bool PTableFile::save(const PTable& pTable) {
    try {
        fstream binaryFile(filePath.c_str(), ios::out | ios::binary);

        /* Write the marker */
        binaryFile.write((char*) &FILE_MARKER, sizeof (FILE_MARKER));

        /* Save the sizeof(int) in the current system (this may help to indicate
         * the system architecture where the P-value table file is created) */
        int intSize = sizeof (int);
        binaryFile.write((char*) &intSize, sizeof (int));

        int mSize = pTable.getMSize();
        int nSize = pTable.getNSize();
        int kSize = pTable.getKSize();
        binaryFile.write(reinterpret_cast<char*> (&mSize), sizeof (int));
        binaryFile.write(reinterpret_cast<char*> (&nSize), sizeof (int));
        binaryFile.write(reinterpret_cast<char*> (&kSize), sizeof (int));

        float* dataPtr = pTable.getDataPtr();
        binaryFile.write(reinterpret_cast<char*> (dataPtr), pTable.getNumStoredVals() * sizeof (float));

        binaryFile.close();
        return true;

    } catch (...) {
        return false;
    }
}

const PTableFile::ReadResult PTableFile::tryLoadInto(PTable& pTable) {
    int mSize, nSize, kSize;

    fstream binaryFile(filePath.c_str(), ios::in | ios::binary);
    if (!binaryFile.is_open() || binaryFile.fail() || binaryFile.eof()) {
        return INVALID_FILE;
    }

    try {
        /* Test marker */
        char* marker = new char[sizeof (FILE_MARKER)];
        binaryFile.read(marker, sizeof (FILE_MARKER));
        string strMarker(marker);
        string strFileMarker(FILE_MARKER);
        delete marker;

        if (strMarker.compare(strFileMarker) != 0)
            return INVALID_FILE;


        /* Test the system architecture */
        int intSize;
        binaryFile.read((char*) &intSize, sizeof (int));
        if (intSize != sizeof (int))
            return WRONG_ARCH;


        /* Read table size */
        binaryFile.read((char*) &mSize, sizeof (int));
        binaryFile.read((char*) &nSize, sizeof (int));
        binaryFile.read((char*) &kSize, sizeof (int));

        char lineBreak = ' ';
        if (filePath.length() > 40) lineBreak = '\n';
        Interface::instance()
                << "Loading P-value table from file" << lineBreak
                << "\"" << filePath << "\"\n"
                << Interface::DEFAULT_INDENT << "Table size : "
                << mSize << " * " << nSize << " * " << kSize << "\n"
                << endl
                << "This will take " << PTable::estimateMemNeededInMB(mSize, nSize, kSize)
                << " MB of memory (RAM). Proceed?";

        if (!Interface::instance().yesNoQuestion(Interface::YES)) {
            /* If user says NO */
            return CANCELLED;
        }

    } catch (...) {
        return INVALID_FILE;
    }


    try {
        /* Show a warning */
        Interface::instance()
                << "Please be patient while the table loads into memory.\n"
                << "This will typically take less than twenty seconds.\n";
        Interface::instance().showLog(true);

        /* Load file into RAM */
        pTable.initialize(mSize, nSize, kSize);
        float* dataPtr = pTable.getDataPtr();
        binaryFile.read(reinterpret_cast<char*> (dataPtr),
                pTable.getNumStoredVals() * sizeof (float));

        /* Finalise */
        if (!binaryFile) {
            return FILE_CORRUPT;
        }
        binaryFile.close();
    } catch (...) {
        return FILE_CORRUPT;
    }

    saveFilePath();
    return SUCCESS;
}

bool PTableFile::saveFilePath(void) {
    try {
        string absolutePath = filePath;

        /* Convert relative path to absolute path */
        if (absolutePath.length() > 1 && absolutePath[0] != '/') {
            if (absolutePath.length() > 2 && absolutePath.substr(0, 2) == "./") {
                absolutePath = absolutePath.substr(2);
            }

            string currentDir = Util::getCurrentDir();
            absolutePath = currentDir + "/" + absolutePath;
        }

        /* Save the path */
        return Util::saveConfig(CONF_KEY, absolutePath);

    } catch (...) {
        return false;
    }
}

PTableFile::ReadResult PTableFile::searchAndLoad(PTable& pTable) {
    ReadResult loadResult = INVALID_FILE;

    /* Try to reuse the file of last run */
    loadResult = tryLoadLastUsedFile(pTable);
    if (loadResult == SUCCESS || loadResult == CANCELLED || loadResult == FILE_CORRUPT) {
        return loadResult;
    }

    /* Create a list of search directories */
    vector<string> searchDirList;
    string searchDir = Interface::instance().getProgramDir();
    searchDirList.push_back(searchDir);
    if (searchDir != "./") {
        searchDirList.push_back("./");
    }

    /* Search */
    for (unsigned long i = 0; i < searchDirList.size(); i++) {
        loadResult = searchAndLoad(searchDirList[i], pTable);
        if (loadResult == SUCCESS || loadResult == CANCELLED || loadResult == FILE_CORRUPT) {
            return loadResult;
        }
    }

    return INVALID_FILE;
}

PTableFile::ReadResult PTableFile::searchAndLoad(string searchDir, PTable& pTable) {
    DIR* localDir = opendir(searchDir.c_str());
    if (localDir == NULL) {
        return INVALID_FILE;
    }

    ReadResult loadResult = INVALID_FILE;
    struct dirent* nextDirent;

    while ((nextDirent = readdir(localDir)) != NULL) {
        string fileName(nextDirent->d_name);
        if (searchDir != "./") {
            fileName = searchDir + fileName;
        }

        PTableFile file(fileName);
        loadResult = file.tryLoadInto(pTable);

        if (loadResult == SUCCESS || loadResult == CANCELLED || loadResult == FILE_CORRUPT) {
            closedir(localDir);
            return loadResult;
        }
    }

    closedir(localDir);
    return INVALID_FILE;
}

const PTableFile::ReadResult PTableFile::tryLoadLastUsedFile(PTable& pTable) {
    try {
        string lastUsedFilePath = Util::getConfig(CONF_KEY);
        PTableFile lastUsedFile(lastUsedFilePath);
        return lastUsedFile.tryLoadInto(pTable);

    } catch (...) {
        return FILE_CORRUPT;
    }
}


/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   PTableFile.h
 * CREATED ON:  07 October 2011, 15:48
 * AUTHOR:      Ha Minh Lam
 *
 * DESCRIPTION: This class represent for P-value table files.
 *
 * HISTORY:     Version     Date            Description
 *              1.0         2011-10-07      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   07 October 2011
 */

#ifndef PTABLEFILE_H
#define	PTABLEFILE_H

#include <cassert>
#include <string>
#include "PTable.h"

using namespace std;

class PTableFile {
public:

    /** The values that could be return when reading a P-value table file */
    enum ReadResult {
        INVALID_FILE /** The file is invalid or does not exist. */,
        WRONG_ARCH /** The P-value table file cannot be used on the current system architecture. */,
        CANCELLED /** The reading process is cancelled by the user. */,
        FILE_CORRUPT /** The file is corrupt. */,
        SUCCESS /** The P-value table has been loaded into RAM successfully. */
    };

    explicit PTableFile(string newFilePath) : filePath(newFilePath) {
    }

    ~PTableFile() {
    }

    /**
     * Get the path to the file from which the P-table is loaded.
     * @return The path to the file from which the P-table is loaded.
     */
    const string getFilePath(void) const {
        return filePath;
    }

    /**
     * Indicates whether the file exists already or not.
     * @return <b>true</b> if the file already exists.
     */
    const bool exists(void);

    /**
     * Save the given P-value table into this file.
     * @param pTable    The pointer to the memory location where the table should
     *                  be load onto.
     * @return  <b>true</b> if saving succeed, <b>false</b> otherwise.
     */
    const bool save(const PTable& pTable);

    /**
     * Try to load this file into the given P-value table.
     * @param pTable
     * @return  INVALID_FILE    The file is invalid or does not exist.
     *          WRONG_ARCH      The P-value table file cannot be used on the
     *                          current system architecture.
     *          CANCELLED       The reading process is cancelled by the user.
     *          FILE_CORRUPT    The file is corrupt.
     *          SUCCESS         The P-value table has been loaded into RAM
     *                          successfully.
     */
    const ReadResult tryLoadInto(PTable& pTable);

    /**
     * Search in current directory for a binary file that contains a P-value table.
     * When a P-value table file is found, it will be automatically loaded into RAM
     * by this method.
     * @param pTable
     * @return  INVALID_FILE    No valid file found.
     *          CANCELLED       A valid file is found but the loading process is
     *                          interrupted by the user.
     *          FILE_CORRUPT    A valid file is found (correct marker and
     *                          architecture) but may be corrupt, so the reading
     *                          process cannot be done.
     *          SUCCESS         A valid file has been found and loaded into RAM.
     */
    static ReadResult searchAndLoad(PTable& pTable);

    static const string getConfigKey();


private:

    /**
     * This string will be written at the beginning of every P-value table files.
     * It is helpful to test if a binary file is a P-value table quickly.
     */
    static const char FILE_MARKER[];

    /** The config key */
    static const string CONF_KEY;

    PTableFile(const PTableFile& orig) : filePath(orig.filePath) {
        assert(false); // should never reach here
    }

    PTableFile& operator=(const PTableFile& rhs) {
        if (this != &rhs) {
            assert(false); // should never reach here
        }
        return *this;
    }

    /**
     * Save the path of the file into the config file so that this path can be
     * used for the next run.
     * @return <b>true</b> if the path is successfully saved.
     */
    bool saveFilePath(void);

    /**
     * Try to load the P-value table from the last used file. Normally, the path
     * to this file is stored in a config file by last run.
     * @param pTable
     * @return  INVALID_FILE    The file is invalid or does not exist.
     *          WRONG_ARCH      The P-value table file cannot be used on the
     *                          current system architecture.
     *          CANCELLED       The reading process is cancelled by the user.
     *          FILE_CORRUPT    The file is corrupt.
     *          SUCCESS         The P-value table has been loaded into RAM
     *                          successfully.
     */
    static const ReadResult tryLoadLastUsedFile(PTable& pTable);

    /**
     * Search in the given directory for a binary file that contains a P-value table.
     * When a P-value table file is found, it will be automatically loaded into RAM
     * by this method.
     * @param searchDir The directory to start searching.
     * @param pTable
     * @return  INVALID_FILE    No valid file found.
     *          CANCELLED       A valid file is found but the loading process is
     *                          interrupted by the user.
     *          FILE_CORRUPT    A valid file is found (correct marker and
     *                          architecture) but may be corrupt, so the reading
     *                          process cannot be done.
     *          SUCCESS         A valid file has been found and loaded into RAM.
     */
    static ReadResult searchAndLoad(string searchDir, PTable& pTable);

    /**
     * The path to the file from which this P-value table is loaded.
     */
    string filePath;


};

#endif	/* PTABLEFILE_H */


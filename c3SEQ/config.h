/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   config.h
 * CREATED ON:  15 June 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: This file contains the declaration of all default values of
 *              variables in the 3SEQ program.
 *
 * HISTORY:     Version    Date          Description
 *              1.0        2011-06-15    Created.
 *
 * VERSION:     1.0
 * LAST EDIT:   15 June 2011
 */

#ifndef CONFIG_H
#define	CONFIG_H

#include <string>


/* Default program info */
const std::string DEFAULT_PROGRAM_NAME = "3SEQ";
const std::string PROGRAM_VERSION = "1.2 beta 2 build 140417";
const std::string DEFAULT_PROGRAM_DESC = "Software For Identifying Recombination In Sequence Data";
const std::string AUTHOR_EMAIL = "maciej.boni@ndm.ox.ac.uk, mboni@oucru.org";


/* Default file names */
const std::string DEFAULT_LOG_FILE_NAME = "3s.log";
const std::string DEFAULT_RECOMBINANTS_FILE_NAME = "3s.rec";
const std::string DEFAULT_LONG_RECOMBINANTS_FILE_NAME = "3s.longRec";
const std::string DEFAULT_PVAL_HISTOGRAM_FILE_NAME = "3s.pvalHist";
const std::string DEFAULT_SKIPPED_TRIPLETS_FILE_NAME = "3s.skipped";

const std::string DEFAULT_MATCH_FILE_NAME = "3s.match"; // used in match-run
const std::string DEFAULT_TRIPLET_PS_FILE_EXT = "3s.triplet.eps"; // used in triplet-run


/* Default analysis settings */
const unsigned long DEFAULT_MIN_RECOMBINANT_LENGTH = 100;
const double DEFAULT_REJECT_THRESHOLD = 0.05;


/** The size of the histogram of P-values */
const int PVAL_HISTOGRAM_SIZE = 41; // 0 to 40


/** The directory where all the program data are stored */
const std::string DEFAULT_CONF_DIR = ".3seq";

/** The name of the file where all the program configurations are stored */
const std::string DEFAULT_CONF_FILE = "3seq.conf";


#endif	/* CONFIG_H */

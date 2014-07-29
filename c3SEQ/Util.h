/**
 * Copyright (c) ___. All Rights Reserved.
 *
 * FILE NAME:   util.h
 * CREATED ON:  26 September 2011, 15:28
 * AUTHOR:      lamhm
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date            Description
 *              1.0         26 September 2011      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   26 September 2011
 */

#ifndef UTIL_H
#define	UTIL_H

#include <string>
#include "config.h"

class Util {
public:

    /**
     * Get a key from standard input without showing anything on the output stream.
     * @return  The code of the key.
     */
    static int getKey(void);

    static const std::string getCurrentDir(void);

    static bool saveConfig(std::string key, std::string value);

    static std::string getConfig(std::string key);


private:

    static const std::string getConfigDir(void);
};

#endif	/* UTIL_H */


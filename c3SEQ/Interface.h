/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   Interface.h
 * CREATED ON:  04 August 2011, 11:33
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: This class handles all the output/log/error messages of the 3Seq
 *              program. It is also the interface between the 3Seq-core (console)
 *              and the 3Seq-gui.
 *
 * HISTORY:     Version    Date         Description
 *              1.0        2011-08-04   Created
 *
 * VERSION:     1.0
 * LAST EDIT:   04 August 2011
 */

#ifndef INTERFACE_H
#define	INTERFACE_H

#include <iostream>
#include <cstdio>
#include <cassert>
#include <sstream>
#include <string>

#include "config.h"
#include "TextFile.h"

using namespace std;

/**
 * This class handles all the output/log/error messages of the 3Seq program.
 * It is also the interface between the 3Seq-CORE (console) and the 3Seq-GUI.
 * @note For performance reason, all the methods of this class are NOT VIRTUAL.
 */
class Interface : public stringstream {
public:

    enum ExitSignal {
        SUCCESS, FAIL
    };

    /** The separator between parts of the output */
    static const string SEPARATOR;
    static const string BOLD_SEPARATOR;

    /** Additional indentation in console output */
    static const string DEFAULT_INDENT;

    enum Mode {
        CONSOLE, GUI, SILENT
    };

    enum ThreeStateOpt {
        YES, NO, UNKNOWN
    };

    ~Interface() {
        if (logFile != NULL) delete logFile;
        cleanOutputStreamFormat();
        cleanErrorStreamFormat();
    }

    static Interface& instance(void) {
        static Interface instance;
        return instance;
    }

    /**
     * Clear the string in this string stream.
     */
    void clear(void) {
        str("");
        stringstream::clear();
    };

    /**
     * Set a new program name.
     * @param newProgramName
     */
    void setProgramName(char* newProgramName) {
        programName = newProgramName;
    };

    /**
     * Get the directory that contains the executable file.
     * @return  The directory that contains the executable file.
     */
    string getProgramDir(void) const;

    /**
     * Set new output mode.
     * @param newMode
     */
    void setMode(Mode newMode) {
        mode = newMode;
    };

    /**
     * Save the command line as a string for logging.
     * @param argc
     * @param argv
     */
    void storeCommand(int argc, char** argv) {
        const string OP_SEPARATOR = " ";

        runCommand = "";
        for (int i = 1; i < argc; i++) {
            if (runCommand.length() > 0) {
                runCommand += OP_SEPARATOR;
            }
            runCommand += argv[i];
        }
    };

    /**
     * Initialise the log file.
     * @param logFileName
     */
    void initLogFile(string logFileName);

    /**
     * Show starting message.
     * @param startMessage
     */
    void startProgram(string startMessage);

    /**
     * Show the path to the current associated p-table file.
     */
    void showPTableLink();

    /**
     * Exit the program.
     * @param status    The error code of the program when exiting (0 means no error).
     */
    void throwExitSignal(int status);

    /**
     * Show output message
     * @param autoEndl  indicates if this method should add one more 'endl'
     *                  character at the end of the message.
     */
    void showOutput(bool autoEndl);

    /**
     * Show log message (this method will show nothing in SILENT mode).
     * @param autoEndl  indicates if this method should add one more 'endl'
     *                  character at the end of the message.
     */
    void showLog(bool autoEndl);

    /**
     * Show warning message and return the user's option (if necessary).
     * @param ignoreAndCont If this flag is <b>true</b>, the warning will be
     *                      shown as a normal message. If it is <b>false</b>,
     *                      the message will appear as a yes/no question to see
     *                      what the user want to do. Default value is <b>true</b>.
     * @return  The user's option. If the warning is just shown as a normal
     *          message (the user don't have to choose), this method will
     *          return <b>true</b>.
     */
    bool const showWarning(bool ignoreAndCont = true);

    /**
     * Show error (and exit if necessary).
     * @param autoEndl  indicates if this method should add one more 'endl'
     *                  character at the end of the message.
     * @param exitAfterShowing  Indicates if the program should exit after
     *                          showing the error message; default value is
     *                          <b>true</b>.
     */
    void showError(bool autoEndl, bool exitAfterShowing);

    /**
     * Show the yes/no question and get the answer.
     * @param defaultAnswer In SILENT mode, this will be taken as the user's
     *                      answer if it is YES or NO. The question only shown
     *                      when the <b>defaultAnswer</b> is UNKNOWN.
     * @return  <b>true</b> if the user answers YES, <b>false</b> otherwise.
     */
    bool yesNoQuestion(ThreeStateOpt defaultAnswer);

    /**
     * Initialise and activate the process counter.
     * @param processingText
     * @param minValue
     * @param maxValue
     * @param percentCounter
     * @note    After activating the process counter, trying to show any message
     *          (except ERROR MESSAGE) will throw a false assertion. Before
     *          continuing showing any message, the <b>finishCounting()</b>
     *          or <b>count(currentValue, true)</b> method must be called to
     *          deactivate the process counter.
     */
    void initCounter(string processingText, double minValue,
            double maxValue, bool percentCounter = true);

    /**
     * Show the counting value and deactivate the counter if needed.
     * @param currentValue  The current value that the counter should count.
     * @param isLastValue   The counter will be deactivate if this is <b>true</b>.
     *                      This parameter is useful when breaking the counting
     *                      is necessary. The default value is <b>false</b>.
     */
    void count(double currentValue, bool isLastValue = false);

    /**
     * Show the current value of the counter. This is useful to show the counter
     * for the first time.
     */
    void count(void) {
        count(counter.currentValue, false);
    };

    /**
     * Show the last (maximum) value of the counter and deactivate it.
     */
    void finishCounting();

    /**
     * Get the elapsed time as a string which is formated as "hhh:mm:ss".
     * @return The elapsed time as a string which is formated as "hhh:mm:ss".
     */
    const string getElapsedTime(void) const {
        char timeFormat[] = "%03d:%02d:%02d";
        char cTime[20];

        long elapsedTime = counter.getElapsedSeconds();
        int hour = static_cast<int> (elapsedTime / 3600);
        int minute = static_cast<int> ((elapsedTime / 60) % 60);
        int second = static_cast<int> (elapsedTime % 60);

        sprintf(cTime, timeFormat, hour, minute, second);

        string result(cTime);
        return result;
    }

    /**
     * Get elapsed time (in second) since the counter was initialised.
     * @return Elapsed time (in second) since the counter was initialised.
     */
    const long getElapsedSeconds(void) const {
        return counter.getElapsedSeconds();
    }

    /**
     * Catch the usage of available modes of this 3SEQ program into this stream.
     */
    void catchUsageModes(void);

    /**
     * Catch version info into this stream.
     */
    void catchVersion(void);

    void catchOptions(void);


private:

    struct Counter {
    public:
        static const int DECIMAL_PLACES;

        void initialize(string processingText, double minValue, double maxValue, bool percentCounter);

        const string showValue() const;

        const long getElapsedSeconds(void) const;

        string processingText;
        double minValue;
        double maxValue;
        double currentValue;
        bool isPercentCounter;

    private:
        time_t startTime;
    };


    static const int PRECISION = 15;

    static const string ERROR_FORMAT;
    static const string WARNING_FORMAT;
    static const string OUTPUT_FORMAT;

    static const string YES_STR;
    static const string NO_STR;
    static const string WARNING_STR;
    static const string ERROR_STR;

    static void setOutputStreamFormat(string format) {
        cout << "\e[" << format << "m";
    };

    static void setErrorStreamFormat(string format) {
        cerr << "\e[" << format << "m";
    };

    static void cleanOutputStreamFormat(void) {
        cout << "\e[m";
    };

    static void cleanErrorStreamFormat(void) {
        cerr << "\e[m";
    };

    explicit Interface(void) {
        programName = DEFAULT_PROGRAM_NAME;
        programVersion = PROGRAM_VERSION;
        mode = CONSOLE;
        isCounting = false;

        runDescription = DEFAULT_PROGRAM_NAME;
        runCommand = "";

        logFile = NULL;

        cleanOutputStreamFormat();
        cleanErrorStreamFormat();
    }

    Interface(const Interface& orig) {
        assert(false);
    }

    Interface& operator=(const Interface& rhs) {
        if (this != &rhs) {
            assert(false);
        }

        return *this;
    }

    const string align(string firstIndent);

    /**
     * The output mode.
     * @default CONSOLE
     */
    Mode mode;

    /** The process counter */
    Counter counter;

    /** This flag indicates if the process counter is working or not. */
    bool isCounting;

    string programName;
    string programVersion;

    string runDescription;
    string runCommand;

    TextFile* logFile;
};

#endif	/* INTERFACE_H */


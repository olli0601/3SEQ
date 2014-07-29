/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   stats.h
 * CREATED ON:  23 August 2011, 16:45
 * AUTHOR:      Maciej F. Boni
 * 
 * DESCRIPTION: Implement statistical methods.
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-08-23   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   23 August 2011
 * 
 * NOTE:        All the methods here are implemented with
 *              "long double" type for more precision.
 */

#ifndef STATS_H
#define	STATS_H

namespace stats {

    namespace siegmund {

        const long double PI = 3.14159265358979323846264338327950288419716939937510L;

        long double continuousApprox(int m, int n, int k);

        long double discreteApprox(int m, int n, int k);

        /**
         * Calculate 1 - exp(x) where x is negative and very close to 0;
         * @param x
         * @return 
         */
        long double oneSubtractExpX(long double x);

        long double normPdf(long double x);

        long double nu(long double x);

        long double normCdf(long double x);
    };


    
    namespace correction {

        /**
         * Dunn-Sidak statistical correction for multiple comparisons
         * @param pVal              The P-value of the current sample.
         * @param numTotalSamples   Total number of samples.
         * @return  The statistically corrected P-value.
         * @note    If the given P-value is too small, the result will be
         *          the same as using Bonferroni correction.
         */
        long double dunnSidak(long double pVal, long double numTotalSamples);

        /**
         * Bonferroni statistical correction for multiple comparisons
         * @param pVal              The P-value of the current sample.
         * @param numTotalSamples   Total number of samples.
         * @return  The statistically corrected P-value.
         */
        long double bonferroni(long double pVal, long double numTotalSamples);
    };

}

#endif	/* STATS_H */


/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   stats.cpp
 * CREATED ON:  23 August 2011, 16:45
 * AUTHOR:      Maciej F. Boni
 *
 * DESCRIPTION: See "stats.h"
 *
 * HISTORY:     Version    Date         Description
 *              1.0        2011-08-23   created
 *
 * VERSION:     1.0
 * LAST EDIT:   23 August 2011
 */

#include "stats.h"
#include <cmath>
#include <cassert>


namespace stats {

    ////////////////////////////////////////////////////////////////////////////////
    //  SIEGMUND
    ////////////////////////////////////////////////////////////////////////////////

    namespace siegmund {

        long double continuousApprox(int m, int n, int k) {
            long double mPlusN = static_cast<long double> (m) + static_cast<long double> (n);
            long double b = (static_cast<long double> (k)) - 0.5L;
            long double xi = static_cast<long double> (n - m);

            long double p1 = exp(-2.0L * b * (b - xi) / mPlusN);

            long double p2 = p1 * (2.0L * (2.0L * b - xi) * (b - xi) / mPlusN + 1.0L);

            //return (1.0L - exp(-p2));
            return oneSubtractExpX(-p2);
        }

        long double discreteApprox(int m, int n, int k) {
            long double mPlusN = static_cast<long double> (m) + static_cast<long double> (n);
            long double b = (static_cast<long double> (k)) - 0.5L;
            long double xi = static_cast<long double> (n - m);

            long double p1 = exp(-2.0L * b * (b - xi) / mPlusN);

            long double p2 = p1 * (2.0L * (2.0L * b - xi)*(b - xi) / mPlusN + 1.0L);

            long double NU = nu(2.0L * (2.0L * b - xi) / mPlusN);

            long double p3 = NU * NU * p2;

            //return (1.0L - exp(-p3));
            return oneSubtractExpX(-p3);
        }

        /**
         * Calculate <i> 1 - exp(x) </i>.
         * @param x A number
         * @return <i>1 - exp(x)</i>.
         * @Note
         * If <i> x >= 0 </i>, the ordinary exponential function of C will be used.
         * Otherwise, the result will be calculated by using the Taylor series: <br>
         * <i> 1 - exp(x) = - x - (x^2)/(2!) - (x^3)/(3!) - (x^4)/(4!) - ... </i><br>
         */
        long double oneSubtractExpX(long double x) {
            /* The number of terms in the Taylor series.
             * The number is chosen based on the fact that floating point arithmetic
             * is considered accurate up to 15 digits. */
            static const int MAX_TERMS_COUNT = 50;

            if (x >= 0.0L) {
                // If x is non-negative, the ordinary exp() function of C will be used.
                return 1.0L - exp(x);

            } else if (x >= -1.0L) {
                long double seriesTerms[MAX_TERMS_COUNT];
                long double xPower = 1.0L; // = x^0
                long double factorial = 1.0L; //  = 0!

                for (int i = 0; i < MAX_TERMS_COUNT; i++) {
                    xPower *= x; //  = x^(i+1)
                    factorial *= static_cast<long double> (i + 1); //  = (i+1)!

                    seriesTerms[i] = xPower / factorial;
                }

                long double result = 0.0L;

                for (int i = MAX_TERMS_COUNT - 1; i >= 0; i--) {
                    /* The order of adding the terms is important to reduce precision loss */
                    result -= seriesTerms[i];
                }

                return result;

            } else {
                // x < 0 && |x| > 1
                long double absX = -x;
                long double expAbsX = exp(absX);
                return (expAbsX - 1.0L) / expAbsX;
            }
        }

        long double normPdf(long double x) {
            return ( 1.0L / sqrt(2.0L * PI)) * exp(-0.5L * x * x);
        }

        long double nu(long double x) {
            return ((normCdf(x / 2.0L) - 0.5L) * 2.0L) / x / (normPdf(x / 2.0L) + (x * normCdf(x / 2.0L)) / 2.0L);
        }

        long double normCdf(long double x) {
            static const long double B1 = 0.319381530L;
            static const long double B2 = -0.356563782L;
            static const long double B3 = 1.781477937L;
            static const long double B4 = -1.821255978L;
            static const long double B5 = 1.330274429L;
            static const long double P = 0.2316419L;
            static const long double C = 0.39894228L;

            if (x >= 0.0L) {
                long double t = 1.0L / (1.0L + P * x);
                return (
                        1.0L - C * exp(-x * x / 2.0L) * t
                        * (t * (t * (t * (t * B5 + B4) + B3) + B2) + B1)
                        );

            } else {
                long double t = 1.0L / (1.0L - P * x);
                return (
                        C * exp(-x * x / 2.0L) * t
                        * (t * (t * (t * (t * B5 + B4) + B3) + B2) + B1)
                        );
            }
        }
    };




    ////////////////////////////////////////////////////////////////////////////////
    //  CORRECTION
    ////////////////////////////////////////////////////////////////////////////////

    namespace correction {

        long double dunnSidak(long double pVal, long double numTotalSamples) {
            /* The minimum value of pVal that can be applied Dunn-Sidak correction */
            static const long double DS_THRESHOLD = 1e-15L;

            if (pVal >= 1.0L) return 1.0L;

            if (pVal < DS_THRESHOLD) {
                /* Simply return the Bonferroni correction because if the P-value is
                 * too small, then 1 - PValue rounds off to one and you cannot do a
                 * Dunn-Sidak correction. */
                return bonferroni(pVal, numTotalSamples);

            } else {
                return 1.0L - pow(1.0L - pVal, numTotalSamples);
            }
        }

        long double bonferroni(long double pVal, long double numTotalSamples) {
            return pVal * numTotalSamples;
        }
    }

}
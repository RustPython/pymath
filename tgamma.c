/*
 * Implementation of the gamma function based on CPython's implementation.
 * This is a standalone version without dependencies on CPython.
 */

#include <math.h>
#include <errno.h>
#include <float.h>
#include <assert.h>
#include <stdbool.h>

/* Constants and definitions needed for gamma calculation */
static const double pi = 3.141592653589793238462643383279502884197;
static const double logpi = 1.144729885849400174143427351353058711647;

/* Check if a value is finite */
#ifndef Py_IS_FINITE
#define Py_IS_FINITE(x) isfinite(x)
#endif

/* Check if a value is infinity */
#ifndef Py_IS_INFINITY
#define Py_IS_INFINITY(x) isinf(x)
#endif

/* Check if a value is NaN */
#ifndef Py_IS_NAN
#define Py_IS_NAN(x) isnan(x)
#endif

/* Define HUGE_VAL if not available */
#ifndef Py_HUGE_VAL
#define Py_HUGE_VAL HUGE_VAL
#endif

/* Define mathematical constants required for Lanczos approximation */
#define LANCZOS_N 13
static const double lanczos_g = 6.024680040776729583740234375;
static const double lanczos_g_minus_half = 5.524680040776729583740234375;
static const double lanczos_num_coeffs[LANCZOS_N] = {
    23531376880.410759688572007674451636754734846804940,
    42919803642.649098768957899047001988850926355848959,
    35711959237.355668049440185451547166705960488635843,
    17921034426.037209699919755754458931112671403265390,
    6039542586.3520280050642916443072979210699388420708,
    1439720407.3117216736632230727949123939715485786772,
    248874557.86205415651146038641322942321632125127801,
    31426415.585400194380614231628318205362874684987640,
    2876370.6289353724412254090516208496135991145378768,
    186056.26539522349504029498971604569928220784236328,
    8071.6720023658162106380029022722506138218516325024,
    210.82427775157934587250973392071336271166969580291,
    2.5066282746310002701649081771338373386264310793408
};

/* denominator is x*(x+1)*...*(x+LANCZOS_N-2) */
static const double lanczos_den_coeffs[LANCZOS_N] = {
    0.0, 39916800.0, 120543840.0, 150917976.0, 105258076.0, 45995730.0,
    13339535.0, 2637558.0, 357423.0, 32670.0, 1925.0, 66.0, 1.0};

/* gamma values for small positive integers, 1 though NGAMMA_INTEGRAL */
#define NGAMMA_INTEGRAL 23
static const double gamma_integral[NGAMMA_INTEGRAL] = {
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
    1307674368000.0, 20922789888000.0, 355687428096000.0,
    6402373705728000.0, 121645100408832000.0, 2432902008176640000.0,
    51090942171709440000.0, 1124000727777607680000.0,
};

/* Helper function for sin(pi*x) used in gamma calculation */
static double m_sinpi(double x)
{
    double y, r;
    int n;
    /* this function should only ever be called for finite arguments */
    assert(Py_IS_FINITE(x));
    y = fmod(fabs(x), 2.0);
    n = (int)round(2.0*y);
    assert(0 <= n && n <= 4);
    switch (n) {
    case 0:
        r = sin(pi*y);
        break;
    case 1:
        r = cos(pi*(y-0.5));
        break;
    case 2:
        /* N.B. -sin(pi*(y-1.0)) is *not* equivalent: it would give
           -0.0 instead of 0.0 when y == 1.0. */
        r = sin(pi*(1.0-y));
        break;
    case 3:
        r = -cos(pi*(y-1.5));
        break;
    case 4:
        r = sin(pi*(y-2.0));
        break;
    default:
        /* Should never happen due to the assert above */
        r = 0.0;
    }
    return copysign(1.0, x)*r;
}

/* Implementation of Lanczos sum for gamma function */
static double lanczos_sum(double x)
{
    double num = 0.0, den = 0.0;
    int i;
    assert(x > 0.0);
    /* evaluate the rational function lanczos_sum(x). For large
       x, the obvious algorithm risks overflow, so we instead
       rescale the denominator and numerator of the rational
       function by x**(1-LANCZOS_N) and treat this as a
       rational function in 1/x. This also reduces the error for
       larger x values. The choice of cutoff point (5.0 below) is
       somewhat arbitrary; in tests, smaller cutoff values than
       this resulted in lower accuracy. */
    if (x < 5.0) {
        for (i = LANCZOS_N; --i >= 0; ) {
            num = num * x + lanczos_num_coeffs[i];
            den = den * x + lanczos_den_coeffs[i];
        }
    }
    else {
        for (i = 0; i < LANCZOS_N; i++) {
            num = num / x + lanczos_num_coeffs[i];
            den = den / x + lanczos_den_coeffs[i];
        }
    }
    return num/den;
}

/**
 * Computes the gamma function value for x.
 *
 * The implementation is based on the Lanczos approximation with parameters 
 * N=13 and g=6.024680040776729583740234375, which provides excellent accuracy
 * across the domain of the function.
 *
 * @param x The input value
 * @return The gamma function value
 *
 * Special cases:
 * - If x is NaN, returns NaN
 * - If x is +Inf, returns +Inf
 * - If x is -Inf, sets errno to EDOM and returns NaN
 * - If x is 0, sets errno to EDOM and returns +/-Inf (with the sign of x)
 * - If x is a negative integer, sets errno to EDOM and returns NaN
 */
double tgamma(double x) {
    double absx, r, y, z, sqrtpow;

    /* special cases */
    if (!Py_IS_FINITE(x)) {
        if (Py_IS_NAN(x) || x > 0.0)
            return x;  /* tgamma(nan) = nan, tgamma(inf) = inf */
        else {
            errno = EDOM;
            return NAN;  /* tgamma(-inf) = nan, invalid */
        }
    }
    if (x == 0.0) {
        errno = EDOM;
        /* tgamma(+-0.0) = +-inf, divide-by-zero */
        return copysign(INFINITY, x);
    }

    /* integer arguments */
    if (x == floor(x)) {
        if (x < 0.0) {
            errno = EDOM;  /* tgamma(n) = nan, invalid for */
            return NAN; /* negative integers n */
        }
        if (x <= NGAMMA_INTEGRAL)
            return gamma_integral[(int)x - 1];
    }
    absx = fabs(x);

    /* tiny arguments:  tgamma(x) ~ 1/x for x near 0 */
    if (absx < 1e-20) {
        r = 1.0/x;
        if (Py_IS_INFINITY(r))
            errno = ERANGE;
        return r;
    }

    /* large arguments: assuming IEEE 754 doubles, tgamma(x) overflows for
       x > 200, and underflows to +-0.0 for x < -200, not a negative
       integer. */
    if (absx > 200.0) {
        if (x < 0.0) {
            return 0.0/m_sinpi(x);
        }
        else {
            errno = ERANGE;
            return Py_HUGE_VAL;
        }
    }

    y = absx + lanczos_g_minus_half;
    /* compute error in sum */
    if (absx > lanczos_g_minus_half) {
        /* note: the correction can be foiled by an optimizing
           compiler that (incorrectly) thinks that an expression like
           a + b - a - b can be optimized to 0.0. This shouldn't
           happen in a standards-conforming compiler. */
        double q = y - absx;
        z = q - lanczos_g_minus_half;
    }
    else {
        double q = y - lanczos_g_minus_half;
        z = q - absx;
    }
    z = z * lanczos_g / y;
    if (x < 0.0) {
        r = -pi / m_sinpi(absx) / absx * exp(y) / lanczos_sum(absx);
        r -= z * r;
        if (absx < 140.0) {
            r /= pow(y, absx - 0.5);
        }
        else {
            sqrtpow = pow(y, absx / 2.0 - 0.25);
            r /= sqrtpow;
            r /= sqrtpow;
        }
    }
    else {
        r = lanczos_sum(absx) / exp(y);
        r += z * r;
        if (absx < 140.0) {
            r *= pow(y, absx - 0.5);
        }
        else {
            sqrtpow = pow(y, absx / 2.0 - 0.25);
            r *= sqrtpow;
            r *= sqrtpow;
        }
    }
    if (Py_IS_INFINITY(r))
        errno = ERANGE;
    return r;
}

/**
 * Computes the natural logarithm of the absolute value of the gamma function.
 *
 * @param x The input value
 * @return The log of the absolute gamma function value
 *
 * Special cases:
 * - If x is NaN, returns NaN
 * - If x is +/-Inf, returns +Inf
 * - If x is a non-positive integer, sets errno to EDOM and returns +Inf
 * - If x is 1 or 2, returns 0.0
 */
double lgamma(double x) {
    double r;
    double absx;

    /* special cases */
    if (!Py_IS_FINITE(x)) {
        if (Py_IS_NAN(x))
            return x;  /* lgamma(nan) = nan */
        else
            return HUGE_VAL; /* lgamma(+-inf) = +inf */
    }

    /* integer arguments */
    if (x == floor(x) && x <= 2.0) {
        if (x <= 0.0) {
            errno = EDOM;  /* lgamma(n) = inf, divide-by-zero for */
            return HUGE_VAL; /* integers n <= 0 */
        }
        else {
            return 0.0; /* lgamma(1) = lgamma(2) = 0.0 */
        }
    }

    absx = fabs(x);
    /* tiny arguments: lgamma(x) ~ -log(fabs(x)) for small x */
    if (absx < 1e-20)
        return -log(absx);

    /* Lanczos' formula. We could save a fraction of a ulp in accuracy by
       having a second set of numerator coefficients for lanczos_sum that
       absorbed the exp(-lanczos_g) term, and throwing out the lanczos_g
       subtraction below; it's probably not worth it. */
    r = log(lanczos_sum(absx)) - lanczos_g;
    r += (absx - 0.5) * (log(absx + lanczos_g - 0.5) - 1);
    if (x < 0.0)
        /* Use reflection formula to get value for negative x. */
        r = logpi - log(fabs(m_sinpi(absx))) - log(absx) - r;
    if (Py_IS_INFINITY(r))
        errno = ERANGE;
    return r;
}

#include <stdio.h>
#include <stdint.h>

union Result {
    double d;
    uint64_t u;
};

int main() {
    union Result result;
    result.d = tgamma(-3.8510064710745118);
    printf("The result of tgamma(-3.8510064710745118) is: %f\n", result.d);

    // Print the binary representation of a double
    printf("Bit representation of result: %llu\n", result.u);
    return 0;
}
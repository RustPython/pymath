#include <stdio.h>
#include <stdint.h>
#include <string.h> // For memcpy
#include <math.h> // For M_PI if needed, though PI is defined below

// Function to print the bit representation of a double
void print_double_bits(const char *name, double val) {
    uint64_t bits;
    // Use memcpy to safely copy the bits, avoiding potential strict-aliasing issues
    memcpy(&bits, &val, sizeof(double));
    printf("%s = %.17g (0x%016lx)\n", name, val, bits);
}

int main() {
    // Constants from gamma.rs
    const double PI = 3.141592653589793238462643383279502884197;
    const double LOG_PI = 1.144729885849400174143427351353058711647;
    const double LANCZOS_G = 6.024680040776729583740234375;
    const double LANCZOS_G_MINUS_HALF = 5.524680040776729583740234375;

    const int LANCZOS_N = 13;
    const double LANCZOS_NUM_COEFFS[LANCZOS_N] = {
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
        2.5066282746310002701649081771338373386264310793408,
    };
    const double LANCZOS_DEN_COEFFS[LANCZOS_N] = {
        0.0,
        39916800.0,
        120543840.0,
        150917976.0,
        105258076.0,
        45995730.0,
        13339535.0,
        2637558.0,
        357423.0,
        32670.0,
        1925.0,
        66.0,
        1.0,
    };

    const int NGAMMA_INTEGRAL = 23;
    const double GAMMA_INTEGRAL[NGAMMA_INTEGRAL] = {
        1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
        362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0,
        87178291200.0, 1307674368000.0, 20922789888000.0,
        355687428096000.0, 6402373705728000.0, 121645100408832000.0,
        2432902008176640000.0, 51090942171709440000.0,
        1124000727777607680000.0,
    };

    printf("--- Single Constants ---\n");
    print_double_bits("PI", PI); // Added PI for completeness
    print_double_bits("LOG_PI", LOG_PI);
    print_double_bits("LANCZOS_G", LANCZOS_G);
    print_double_bits("LANCZOS_G_MINUS_HALF", LANCZOS_G_MINUS_HALF);

    printf("\n--- LANCZOS_NUM_COEFFS ---\n");
    for (int i = 0; i < LANCZOS_N; ++i) {
        char name[32];
        snprintf(name, sizeof(name), "LANCZOS_NUM_COEFFS[%d]", i);
        print_double_bits(name, LANCZOS_NUM_COEFFS[i]);
    }

    printf("\n--- LANCZOS_DEN_COEFFS ---\n");
    for (int i = 0; i < LANCZOS_N; ++i) {
        char name[32];
        snprintf(name, sizeof(name), "LANCZOS_DEN_COEFFS[%d]", i);
        print_double_bits(name, LANCZOS_DEN_COEFFS[i]);
    }

    printf("\n--- GAMMA_INTEGRAL ---\n");
    for (int i = 0; i < NGAMMA_INTEGRAL; ++i) {
        char name[32];
        snprintf(name, sizeof(name), "GAMMA_INTEGRAL[%d]", i);
        print_double_bits(name, GAMMA_INTEGRAL[i]);
    }

    return 0;
}

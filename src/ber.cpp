#include "ber.hpp"

CArray ber(CArray &pulse, float fs, int snr, int nbits, int M, float omega)
{
    int N = nbits/2;
    double coefficient = pow(10, snr/20);

    RArray nbits_linspace (nbits), sumaT2 (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        if (i != 0)
        {
            nbits_linspace[i] = i-N;
        }
    }
    for (int c=0; c < OFFSETS.size(); c++)
    {
        double g0 = coefficient;
    }
}

// -- END OF FILE -- //
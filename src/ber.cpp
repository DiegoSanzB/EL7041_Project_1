#include "ber.hpp"

valarray<double> ber_isi(Pulse &pulse, float fs, int snr, int nbits, int M, float omega)
{
    int N = nbits/2;
    double snr_coefficient = pow(10, snr/20);

    valarray<double> nbits_linspace (nbits), p_isi (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        if (i != 0)
        {
            // [-N, -N+1, ... , -2, -1, 1, 2, ... , N-1, N]
            nbits_linspace[i] = i-N;
        }
    }
    for (int c = 0; c<OFFSETS.size(); c++)
    {
        double g0 = snr_coefficient*pulse.evaluate(OFFSETS[c]*time_unit_period);
        valarray<double> gk (nbits);

        // Calculate values for gk
        for (int i=0; i<nbits; i++)
        {

            gk[i] = snr_coefficient * ((int)rand() % 2 + 1) * pulse.evaluate((OFFSETS[c] - nbits_linspace[i])*time_unit_period);
        }
        // Calculate the sum and product
        double sum = 0;
        for (int m=1; m<=M; m += 2)
        {
            double mult = 1;
            // Mult
            for (int k=0; k<nbits; k++) // gk.size() == nbits
            {
                mult *= (double)(cos(m*omega*gk[k]));
            }
            sum += (double)((exp(-pow(m*omega, 2)/2) * sin(m*omega*g0)) / m) * mult;
        }
        p_isi[c] = (double)(0.5 - 2.0/PI * sum);
    }
    return p_isi;
}

// -- END OF FILE -- //
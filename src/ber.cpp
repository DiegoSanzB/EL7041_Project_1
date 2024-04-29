#include "ber.hpp"

valarray<double> ber_isi(Pulse pulse, double fs, double snr, int nbits, int M, double omega)
{
    int N = floor(nbits/2);
    double snr_coefficient = pow(10, snr/20.0);

    valarray<double> nbits_linspace (nbits), p_isi (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        if (i < N)
        {
            // [-N, -N+1, ... , -2, -1, 1, 2, ... , N-1, N]
            nbits_linspace[i] = i-N;
        } else
        {
            nbits_linspace[i] = i-N+1;
        }
    }
    for (int c = 0; c<OFFSETS.size(); c++)
    {
        double g0 = snr_coefficient*pulse.evaluate(OFFSETS[c]*time_unit_period);
        valarray<double> gk (nbits);

        // Calculate values for gk
        for (int i=0; i<nbits; i++)
        {
            double pm_one = ((int)rand() % 2)*2 - 1;
            gk[i] = pm_one * snr_coefficient * pulse.evaluate((OFFSETS[c] - nbits_linspace[i])*time_unit_period);
            // cout << "i = " << i << " -> " << nbits_linspace[i] << endl;
        }
        // Calculate the sum and product
        double sum = 0.0, mult = 1.0;
        for (double m=1.0; m<=M; m += 2.0)
        {
            for (int k=0; k<nbits; k++) // gk.size() == nbits
            {
                // if (abs(m*omega*gk[k]) > lim)
                mult = mult * cos(m*omega*gk[k]);
                // cout << mult << endl;
                // cout << "k = " << k << " -> " << cos(m*omega*gk[k]) << endl;
            }
            // cout << m << " : " << mult << endl;
            sum = sum + (exp(-0.5*pow(m*omega, 2)) * sin(m*omega*g0) / m) * mult;
            // cout << sum << endl;
            mult = 1;
        }
        p_isi[c] = (double)(0.5 - 2.0/PI * sum);
    }
    return p_isi;
}

// -- END OF FILE -- //
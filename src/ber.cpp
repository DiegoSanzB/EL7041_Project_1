#include "ber.hpp"

valarray<double> ber_isi(Pulse pulse, double fs, double snr, int nbits, int M, double truncation = numeric_limits<double>::infinity())
{
    int N = floor(nbits/2);
    double snr_coefficient = pow(10, snr/20.0);
    double omega = 1.0 / fs;

    valarray<double> nbits_linspace (nbits), p_isi (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        // [-N, -N+1, ... , -2, -1, 1, 2, ... , N-1, N]
        if (i < N)
        {
            nbits_linspace[i] = i-N;
        } else
        {
            nbits_linspace[i] = i-N+1;
        }
    }
    for (int c = 0; c<OFFSETS.size(); c++)
    {
        double g0 = snr_coefficient*pulse.evaluate(OFFSETS[c]*time_unit_period, truncation*fs);
        valarray<double> gk (nbits);

        // Calculate values for gk
        for (int i=0; i<nbits; i++)
        {
            double pm_one = ((int)rand() % 2)*2 - 1;
            gk[i] = pm_one * snr_coefficient * pulse.evaluate((OFFSETS[c] - nbits_linspace[i])*time_unit_period, truncation*fs);
        }
        // Calculate the sum and product
        double sum = 0.0, mult = 1.0;
        for (double m=1.0; m<=M; m += 2.0)
        {
            for (int k=0; k<nbits; k++) // gk.size() == nbits
            {
                mult *= cos(m*omega*gk[k]);
            }
            sum +=(exp(-0.5*pow(m*omega, 2)) * sin(m*omega*g0) / m) * mult;
            mult = 1;
        }
        p_isi[c] = (double)(0.5 - 2.0/PI * sum);
    }
    return p_isi;
}


valarray<double> ber_cci(Pulse pulse, double fs, double snr, double sir, int L, int nbits, int M, double truncation = numeric_limits<double>::infinity())
{
    int N = floor(nbits/2);
    double snr_coefficient = pow(10, snr/20.0);
    double sir_coefficient = pow(10, sir/20.0);
    double omega = 1 / fs;

    valarray<double> nbits_linspace (nbits), p_isi (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        // [-N, -N+1, ... , -2, -1, 1, 2, ... , N-1, N]
        if (i < N)
        {
            nbits_linspace[i] = i-N;
        } else
        {
            nbits_linspace[i] = i-N+1;
        }
    }
    for (int c = 0; c<OFFSETS.size(); c++)
    {
        double g0 = snr_coefficient*pulse.evaluate(OFFSETS[c]*time_unit_period, truncation*fs);

        // Calculate r_i assuming each co-channel has same power
        double r_i = sqrt(g0*g0/(L*sir_coefficient));   // TODO: coefficient seems off given P_CCI values in paper

        // Calculate the sum and product
        double sum = 0.0, mult = 1.0;
        for (double m=1.0; m<=M; m += 2.0)
        {
            for (int i=1; i<=L; i++) // L co-channels
            {
                mult *= cyl_bessel_j(0, (double)(m*omega*r_i));
            }
            sum += (exp(-0.5*pow(m*omega, 2)) * sin(m*omega*g0) / m) * mult;
            mult = 1;
        }
        p_isi[c] = (double)(0.5 - 2.0/PI * sum);
    }
    return p_isi;
}


valarray<double> ber_isi_cci(Pulse pulse, double fs, double snr, double sir, int L, int nbits, int M, double truncation = numeric_limits<double>::infinity())
{
    int N = floor(nbits/2);
    double snr_coefficient = pow(10, snr/20.0);
    double sir_coefficient = pow(10, sir/20.0);
    double omega = 1 / fs;

    valarray<double> nbits_linspace (nbits), p_isi (OFFSETS.size());

    for (int i=0; i <= nbits; i++)
    {
        // [-N, -N+1, ... , -2, -1, 1, 2, ... , N-1, N]
        if (i < N)
        {
            nbits_linspace[i] = i-N;
        } else
        {
            nbits_linspace[i] = i-N+1;
        }
    }
    for (int c = 0; c<OFFSETS.size(); c++)
    {
        double g0 = snr_coefficient*pulse.evaluate(OFFSETS[c]*time_unit_period, truncation*fs);
        valarray<double> gk (nbits);

        double r_i = sqrt(g0*g0/(L*sir_coefficient));

        // Calculate values for gk
        for (int i=0; i<nbits; i++)
        {
            double pm_one = ((int)rand() % 2)*2 - 1;
            gk[i] = pm_one * snr_coefficient * pulse.evaluate((OFFSETS[c] - nbits_linspace[i])*time_unit_period, truncation*fs);
        }
        // Calculate the sum and product
        double sum = 0.0, mult = 1.0;
        for (double m=1.0; m<=M; m += 2.0)
        {
            for (int k=0; k<nbits; k++) // isi
            {
                mult *= cos(m*omega*gk[k]);
            }
            for (int i=1; i<=L; i++) // L co-channels
            {
                mult *= cyl_bessel_j(0, (double)(m*omega*r_i));
            }

            sum +=(exp(-0.5*pow(m*omega, 2)) * sin(m*omega*g0) / m) * mult;
            mult = 1;
        }
        p_isi[c] = (double)(0.5 - 2.0/PI * sum);
    }
    return p_isi;
}

// -- END OF FILE -- //
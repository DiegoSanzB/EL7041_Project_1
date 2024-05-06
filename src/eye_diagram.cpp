#include "eye_diagram.hpp"

valarray<double> eye_diagram(size_t N, valarray<double> pulse, float fs)
{
    valarray<int> am = random_bpsk(N);
    valarray<int> amU = bpsk_upsample(am, N, fs);

    valarray<double> out, am_upsampled_complex (N*fs);

    for (int i = 0; i < N*fs; i++)
    {
        am_upsampled_complex[i] = amU[i];
    }

    out = conv(am_upsampled_complex, pulse);

    return out;
}


valarray<int> random_bpsk(size_t N)
{
    srand(time(NULL));

    valarray<int> out (N);

    for (size_t i = 0; i < N; i++)
    {
        int num = ((int)rand() % 2);
        out[i] = num;
    }
    out = 2*out - 1;

    return out;
}


valarray<int> bpsk_upsample(valarray<int> &bpsk, size_t N, float fs)
{
    valarray<int> out (N*fs);

    for (size_t i = 0; i < N*fs; i++)
    {
        if (i % (int)fs == 0)
        {
            out[i] = bpsk[(int)(i/fs)];
        } else
        {
            out[i] = 0;
        }
    }
    return out;
}


valarray<double> conv(valarray<double> &a, valarray<double> &b)
{
    size_t n = a.size(), m = b.size();
    valarray<double> r(n + m -1);

    for (size_t i=0; i<r.size(); i++)
    {
        r[i] = 0;
        size_t s = (i+2-m > 0) ? i+2-m : 0;
        for (size_t j=s; j < min(n, i+2); j++)
        {
            r[i] += a[j] * b[m - (i-j) - 1];
        }
    }
    return r;
}

// -- END OF FILE -- //
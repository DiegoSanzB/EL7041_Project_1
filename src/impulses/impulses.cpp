#include "impulses/impulses.hpp"

CArray raised_cosine(float alpha)
{
    // DEFAULT
    float fs = 10;
    int NFFT = 1024;
    return raised_cosine(alpha, fs, NFFT);    
}

CArray raised_cosine(float alpha, float fs, int NFFT)
{
    CArray SampledPulse (2*fs*fs + 1);
    for (int i = 0; i <= 2*fs*fs; i++)
    {
        SampledPulse[i] = -fs + i/fs;
    }
    // SINC
    CArray sincNum = sin(PI*SampledPulse);
    CArray sincDen = PI*SampledPulse;

    CArray sincOp (2*fs*fs + 1);
    
    for (int i = 0; i <= 2*fs*fs; i++)
    {
        if (abs(sincDen[i].real()) > pow(10, -20))
        {
            sincOp[i] = sincNum[i] / sincDen[i];
        } else
        {
            sincOp[i] = 1;
        }
    }
    // COSINE
    CArray cosNum = cos(alpha*PI*SampledPulse);
    CArray cosDen = (1 - pow(2*alpha*SampledPulse, 2.0));

    CArray cosOp (2*fs*fs + 1);

    for (int i = 0; i <= 2*fs*fs; i++)
    {
        if (abs(cosDen[i].real()) > pow(10, -20))
        {
            cosOp[i] = cosNum[i] / cosDen[i];
        } else
        {
            cosOp[i] = PI/4;
        }
    }
    // RAISED COSINE
    CArray RC_Filter = sincOp * cosOp;

    return RC_Filter;
}


CArray BTRC(float alpha)
{

}


CArray IPLCP(float alpha, float mu, float epsilon, float gamma)
{

}

CArray ELP(float alpha, float beta)
{

}

// -- END OF FILE -- //
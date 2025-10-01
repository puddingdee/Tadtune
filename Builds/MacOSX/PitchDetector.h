//
//  PitchDetector.h
//  Tadtune
//
//  Created by Merritt Hyman on 10/1/25.
//

#pragma once

#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <limits>


class PitchDetector
{
public:
    PitchDetector() = default;
    
    // initialization
    void prepare(double sampleRate, int maxPeriodSamples = 880)
    {
        fs = sampleRate;
        maxPeriod = maxPeriodSamples;
        
        //tracking mode arrays
        //once pitch is detected, maintain small arrays around the detected period
        const int N = 8;
        E.resize(N, 0.0); //energy function
        H.resize(N, 0.0); //correlation function
        
        //detection mode arrays
        //downsampled data to search wide frequency range more efficiently
        //44100/8 = 5512.5Hz sample rate, so we're checking L=2: 2756Hz through L=110: 50.1Hz
        Edown.resize(110, 0.0);
        Hdown.resize(110, 0.0);
        
        //input buffer, circular, needs to be large enough to hold multiple periods of the lowest possible frequency
        inputBuffer.resize(maxPeriod * 4, 0.0f);
        
        //downsampled buffer stores 8:1 downsampled data used during detection mode
        downsampledBuffer.resize(512, 0.0f);
        
        reset();
        
        
    }
    
    void reset()
    {
        detectionMode = true;
        inputWritePos = 0;
        downSampleCounter = 0;
        downsampleWritePos = 0;
        cyclePeriod = 100.0;
        EH_offset = 0;
        
        std::fill(E.begin(), E.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(Edown.begin(), Edown.end(), 0.0);
        std::fill(Hdown.begin(), Hdown.end(), 0.0);
    }
    
    //sample processing entry point
    //2 modes: detection mode when a pitch is unknown, search a wide range
    // tracking mode when a pitch is known, track more efficiently with small arrays
    bool processSample(float sample)
    {
        inputBuffer[inputWritePos] = sample;
        inputWritePos = (inputWritePos + 1) % inputBuffer.size();
        
        if (detectionMode)
        {
            //downsample by 8
            downSampleCounter++;
            if (downSampleCounter >= 8){
                downSampleCounter = 0;
                
                //simple average for downsampling
                //TODO implement FIR lowpass anti aliasing
                float avg = 0.0f;
                for (int i = 0; i < 8; i++)
                {
                    int idx = (inputWritePos - 8 + i + inputBuffer.size()) % inputBuffer.size();
                    avg += inputBuffer[idx];
                }
                avg /= 8.0f;
                
                downsampledBuffer[downSampleWritePos] = avg;
                downSampleWritePos = (downSampleWritePos + 1) % downsampledBuffer.size();
                
                //try to detect pitch from downsampled data
                return detectPitch();
            }
        }
        else{
            return trackPitch();
        }
        return false;
    }
    
    //auto correlation functions
    /**
         * Compute E(L) - the energy function
         *
         * E(L) = sum of [x_i^2 + x_(i-L)^2] for i from (i-2L) to i
         *
         * This represents the total energy in two consecutive periods of length L.
         * At the true period, this captures the "steady state" energy.
         *
         * @param L          The lag (prospective period in samples)
         * @param buffer     The audio buffer to analyze
         * @param writePos   Current write position in circular buffer
         * @return           The computed energy value
         */
    //ai implementation of E(L)
    
    double computeE(int L, const std::vector<float>& buffer, int writePos){
        double sum = 0.0;
        int bufSize = buffer.size();
        
        //sum over 2 periods (2L samples)
        for (int k = 0; k < 2 * L; k++)
        {
            // Get sample at position (writePos - 2L + k)
            int idx1 = (writePos - 2 * L + k + bufSize) % bufSize;
            // Get sample L samples earlier (one period back)
            int idx2 = (writePos - L + k + bufSize) % bufSize;
            
            float x1 = buffer[idx1];
            float x2 = buffer[idx2];
            
            // Accumulate x_i^2 + x_(i-L)^2
            sum += x1 * x1 + x2 * x2;
        }
        return sum;
    
            
    }
    double computeH(int L, const std::vector<float>& buffer, int writePos) {
            double sum = 0.0;
            int bufSize = buffer.size();
            
            for (int k = 0; k < 2 * L; k++) {
                int idx1 = (writePos - 2 * L + k + bufSize) % bufSize;
                int idx2 = (writePos - L + k + bufSize) % bufSize;
                
                float x1 = buffer[idx1];
                float x2 = buffer[idx2];
                
                // H(L) = sum of [x_i * x_(i-L)] - correlation function
                sum += x1 * x2;
            }
            return sum;
        }
    /*inital pitch detection
    from DOWNSAMPLED data
    this is equation 6 from the patent
    eps * E(L) > E(L) - 2*H(L)
    we want the smallest L that satisfies this and is also a local minimum
    */
    
    bool detectPitch()
    {
        for (int L = 2; L < 110; L++)
        {
            Edown[L] = computeE(L, downsampledBuffer, downSampleWritePos);
            Hdown[L] = computeH(L, downsampledBuffer, downSampleWritePos);
        }
        
        //periodicity test
        //find first L where (E(L) - 2*H(L) is small relative to E(L)
        const double eps = 0.1;
        int Lmin1 = -1;
        
        for (int L = 2; L < 110; L++) {
            double threshold = eps * Edown[L];
            double diff = Edown[L] - 2.0 * Hdown[L];
            
            if (diff < threshold && diff >= 0) {
                if (L > 2 && L < 109) {
                    double prevDiff = Edown[L-1] - 2.0 * Hdown[L-1];
                    double nextDiff = Edown[L+1] - 2.0 * Hdown[L+1];
                    
                    if (diff < prevDiff && diff < nextDiff) {
                        Lmin1 = L;
                        break;
                    }
                }
            }
        }

    if (Lmin1 < 0) return false; // no periodicity detected
    
    //TODO handle missing fundamental
    
    // scale from downsampled period to full rate period
    int Lmin = 8 * Lmin1;
    
    //init tracking arrays, we don't need to search the whole range anymore
    const int N = E.size();
    EH_offset = Lmin - N / 2 + 1; // center arrays around detected period
    
    //compute E and H at full sample rate around detected period
    for (int i = 0; i < N; i++)
    {
        int L = EH_offset + i;
        E[i] = computeE(L, inputBuffer, inputWritePos);
        H[i] = computeH(L, inputBuffer, inputWritePos);
    }
    
    cyclePeriod = static_cast<double>(Lmin);
    
    //switch to tracking mode
    detectionMode = false;
    
    return true;
    
    
    //continuous pitch tracking
    bool trackPitch()
    {
        //update E and H arrays around current period estimate
        const int N = E.size();
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
            
        }
        
        //find which L gives minimum E(L) - 2*H(L)
        double minVal = std::numeric_limits<double>::max();
        int minIdx = -1;
        
        for (int i = 0; i < N; i++)
        {
            double val = E[i] - 2.0 * H[i];
            if (val < minVal)
            {
                minVal = val;
                minIdx = i;
            }
        }
        if (minIdx < 0) return false;
        
        //tracking validation: check if minimum still satisfies periodicity criterion
        const double eps = 0.1;
        if (minVal > eps * E[minIdx])
        {
            //lost tracking, return to detection mode
            detectionMode = true;
            return false;
        }
        
        //adaptive window shifting
        if (minIdx < N/2)
        {
            //then period decreased, so shift array down to smaller L
            for (int i = N - 1; i > 0; i--)
            {
                E[i] = E[i - 1];
                H[i] = H[i - 1];

            }
            
            EH_offset--; //array represents valuse one smaller
            
            //compute new E[0] and H[0] for new smallest L
            int L = EH_offset + N - 1;
            E[0] = computeE(L, inputBuffer, inputWritePos);
            H[0] = computeH(L, inputBuffer, inputWritePos);
        }
        else if (minIdx > N / 2 + 1)
        {
            //period increased, shift array up
            for (int i = 0; i < N - 1; i++)
            {
                E[i] = E[i + 1];
                H[i] = H[i + 1];
            }
            EH_offset++; // one larger
            
            //compute new for largest L
            int L - EH_offset;
            E[N-1] = computeE(L, inputBuffer, inputWritePos);
            H[N-1] = computeH(L, inputBuffer, inputWritePos);
        }
        
        //integer L is not accurate enough, so we interpolate to get subsample recision using quadratic fit
        double Pmin = interpolateMinimum(minIdx);
        
        //final period in samples with fractional precision
        cyclePeriod = EH_offset + Pmin - 1.0;
        return true;
    }
    
    //subsample precision
    //from ai:
    /**
         * Interpolate the minimum of E(i) - 2*H(i) for fractional period.
         *
         * Uses quadratic interpolation of three points around the integer minimum.
         * This gives sub-sample accuracy needed for high-quality pitch correction.
         *
         * For a parabola through points (x-1, y0), (x, y1), (x+1, y2):
         * The minimum occurs at x + offset where:
         *    offset = (y0 - y2) / (2 * (y0 - 2*y1 + y2))
         */
    
    double interpolateMinimum(int idx)
    {
        const int N = E.size();
        
        //edge case: can't interpolate at boundaries
        if (idx <= 0 || idx >= N - 1)
            return static_cast<double>(idx);
        
        //get three points around minimum
        double y0 = E[idx - 1] - 2.0 * H[idx - 1];
        double y1 = E[idx] - 2.0 * H[idx];
        double y2 = E[idx + 1] - 2.0 * H[idx + 1];
        
        //quadratic fit formula
        double denom = 2.0 * (y0 - 2.0 * y1 + y2);
        
        //avoid division by 0
        if (std::abs(denom) < 1e-10)
            return static_cast<double>(idx);
        double offset = (y0 - y2) / denom;
        
        return idx + offset;
    }
    
    //public methods
    double getCyclePeriod() const { return cyclePeriod; }
    double getFrequency() const
    {
        return (cyclePeriod > 0.0) ? (fs / cyclePeriod) : 0.0;
    }
    bool isInDetectionMode() const { return detectionMode; }
    
private:
    //audio params
    double fs = 44100.0;
    int maxPeriod = 880;
    
    bool detectionMode = true;
    
    //circular buffers
    std::vector<float> inputBuffer;
    std::vector<float> downsampledBuffer;
    
    //auto-correlation arrays
    std::vector<double> E, H;
    std::vector<double> Edown, Hdown;
    
    //buffer management
    int inputWritePos = 0;
    int downSampleWritePos = 0;
    int downSampleCounter = 0;
    
    //tracking state
    double cyclePeriod = 100.0; // current detected period
    int EH_offset = 0; //L value corresponding to E(0) and H(0)
};


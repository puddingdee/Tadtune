//
//  PitchDetector.h
//  Tadtune
//
//

#ifndef PitchDetector_h
#define PitchDetector_h


#endif /* PitchDetector_h */

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
        fs = sampleReate;
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
        Hdown.resize(110.0.0);
        
        //input buffer, circular, needs to be large enough to hold multiple periods of the lowest possible frequency
        inputBuffer.resize(maxPeriod * 4, 0.0f);
        
        //downsampled buffer stores 8:1 downsampled data used during detection mode
        downsampledBuffer.resize(512, 0.0f);
        
        reset();
        
        
    }
    
    void reset()
    {
        detectionMode = truue;
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
            return detectPitch();
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
            Edown[L] = computeE(L, downsampledBuffer, downSampledWritePos);
            Hdown[L] = computeH(L, downsampledBuffer, downSampledWritePos);
        }
    }
    
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


#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <limits>

class PitchDetector
{
public:
    PitchDetector() = default;
    
    //initalization
    void prepare(double sampleRate, int maxPeriodSamples = 880)
    {
        fs = sampleRate;
        maxPeriod = maxPeriodSamples;
        //tracking mode arrays
        // once pitch is detected, small arrays are used to track pitch
        const int N = 8;
        E.resize(N, 0.0); // energy function for tracking
        H.resize(N, 0.0); // correlation function for tracking
        
        //detection mode arrays
        // 44100/8 = 5512.5 Hz sample rate gives us:
        // - L=2:   2756Hz
        // - L=110: 50.1Hz
        Edown.resize(110, 0.0);
        Hdown.resize(110, 0.0);
        
        //circular input buffer
        inputBuffer.resize(maxPeriod * 4, 0.0f);
        
        //downsampled buffer
        downsampledBuffer.resize(512, 0.0f);
        
        reset();
    }
    
    void reset()
    {
        detectionMode = true; //start with detection
        inputWritePos = 0;
        downSampleCounter = 0;
        downSampleWritePos = 0;
        cyclePeriod = 100.0; //init guess
        EH_offset = 0;
        
        std::fill(E.begin(), E.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(Edown.begin(), Edown.end(), 0.0);
        std::fill(Hdown.begin(), Hdown.end(), 0.0);
    }
    
    //sample processing entry point
    bool processSample(float sample)
    {
        //store sample in circular buffer
        inputBuffer[inputWritePos] = sample;
        inputWritePos = (inputWritePos + 1) % inputBuffer.size();
        
        if (detectionMode)
        {
            //downsample by 8 to reduce computation
            downSampleCounter++;
            if (downSampleCounter >= 8)
            {
                downSampleCounter = 0;
                
                //simple averaging for downsampling
                //TODO anti aliasing filter
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
        else
        {
            //change to tracking mode
            return trackPitch();
        }
        
        return false;
    }
    
    //autocorrelation functions, equations from the patent, claude ai implementation
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
    double computeE(int L, const std::vector<float>& buffer, int writePos)
    {
        double sum = 0.0;
        int bufSize = buffer.size();
        
        // Sum over two periods (2L samples)
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
    
    /**
     * Compute H(L) - the cross-correlation function
     *
     * H(L) = sum of [x_i * x_(i-L)] for i from (i-2L) to i
     *
     * This measures how well the waveform correlates with itself L samples ago.
     * At the true period L, consecutive cycles are nearly identical, so:
     * - x_i ≈ x_(i-L), meaning x_i * x_(i-L) ≈ x_i^2
     * - Therefore H(L) ≈ E(L)/2 at the true period
     *
     * @param L          The lag (prospective period in samples)
     * @param buffer     The audio buffer to analyze
     * @param writePos   Current write position in circular buffer
     * @return           The computed correlation value
     */
    double computeH(int L, const std::vector<float>& buffer, int writePos)
    {
        double sum = 0.0;
        int bufSize = buffer.size();
        
        // Sum over two periods (2L samples)
        for (int k = 0; k < 2 * L; k++)
        {
            int idx1 = (writePos - 2 * L + k + bufSize) % bufSize;
            int idx2 = (writePos - L + k + bufSize) % bufSize;
            
            // Accumulate x_i * x_(i-L)
            sum += buffer[idx1] * buffer[idx2];
        }
        return sum;
    }
    
    // initial pitch detection
    bool detectPitch()
    {
        // compute E and H for all positive periods
        for (int L = 2; L < 110; L++)
        {
            Edown[L] = computeE(L, downsampledBuffer, downSampleWritePos);
            Hdown[L] = computeH(L, downsampledBuffer, downSampleWritePos);
        }
        
        // PERIODICITY TEST:
        // find first L where E(L) - 2*H(L) is small relative to E(L)
        const double eps = 0.1;  // threshold where smaller = more strict matching
        int Lmin1 = -1;
        
        for (int L = 2; L < 110; L++)
        {
            double threshold = eps * Edown[L];
            double diff = Edown[L] - 2.0 * Hdown[L];
            
            // must be positive and small
            if (diff < threshold && diff >= 0)
            {
                // Verify it's a local minimum, not just noise
                if (L > 2 && L < 109)
                {
                    double prevDiff = Edown[L-1] - 2.0 * Hdown[L-1];
                    double nextDiff = Edown[L+1] - 2.0 * Hdown[L+1];
                    
                    if (diff < prevDiff && diff < nextDiff)
                    {
                        Lmin1 = L;
                        break;  // found first minimum, this is the period
                    }
                }
            }
        }
        
        if (Lmin1 < 0) return false;  // no periodicity detected
        
        //TODO handle missing fundamental on vowels with weak fundamentals
        
        // scale from downsampled period to full-rate period
        int Lmin = 8 * Lmin1;
        
        // init smaller tracking arrays now that a pitch is detected
        const int N = E.size();  // should be 8
        EH_offset = Lmin - N / 2 + 1;  // center arrays around detected period
        
        // compute E and H at full sample rate around the detected period
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
        }
        
        cyclePeriod = static_cast<double>(Lmin);
        
        // switch to tracking mode
        detectionMode = false;
        
        return true;
    }
    
    //continuous pitch tracking
    bool trackPitch()
    {
        // update E and H arrays around current period estimate
        const int N = E.size();
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
        }
        
        // find which L gives minimum E(L) - 2*H(L)
        // this is the most likely period
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
        
        // TRACKING VALIDATION:
        // check if minimum still satisfies periodicity criterion
        const double eps = 0.1;
        if (minVal > eps * E[minIdx])
        {
            // lost tracking, return to detection mode
            detectionMode = true;
            return false;
        }
        
        // if the minimum is near the edge of our array, shift the array
        if (minIdx < N / 2)
        {
            // period decreased, shift array down
            for (int i = N - 1; i > 0; i--)
            {
                E[i] = E[i - 1];
                H[i] = H[i - 1];
            }
            EH_offset--;  // array now represents L values one smaller
            
            // compute new E[0] and H[0] for the new smallest L
            int L = EH_offset;
            E[0] = computeE(L, inputBuffer, inputWritePos);
            H[0] = computeH(L, inputBuffer, inputWritePos);
        }
        else if (minIdx > N / 2 + 1)
        {
            // period increased, shift array up
            for (int i = 0; i < N - 1; i++)
            {
                E[i] = E[i + 1];
                H[i] = H[i + 1];
            }
            EH_offset++;  // array now represents L values one larger
            
            // compute new E[N-1] and H[N-1] for the new largest L
            int L = EH_offset + N - 1;
            E[N - 1] = computeE(L, inputBuffer, inputWritePos);
            H[N - 1] = computeH(L, inputBuffer, inputWritePos);
        }
        
        //subsample precision with quadratic interpolation
        double Pmin = interpolateMinimum(minIdx);
        
        //final period
        cyclePeriod = EH_offset + Pmin - 1.0;
        
        return true;
    }
    
//    subsample precision
    double interpolateMinimum(int idx)
    {
        const int N = E.size();
        
        //edge case at boundaries
        if (idx <= 0 || idx >= N - 1)
            return static_cast<double>(idx);
        
        //get three points around minimum
        double y0 = E[idx - 1] - 2.0 * H[idx - 1];
        double y1 = E[idx] - 2.0 * H[idx];
        double y2 = E[idx + 1] - 2.0 * H[idx + 1];
        
        // quadratic fit formula
        double denom = 2.0 * (y0 - 2.0 * y1 + y2);
        
        // avoid division by zero
        if (std::abs(denom) < 1e-10)
            return static_cast<double>(idx);
        
        double offset = (y0 - y2) / denom;
        
        // return fractional index
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
    
    // operating mode
    bool detectionMode = true;
    
    // circular buffers
    std::vector<float> inputBuffer;
    std::vector<float> downsampledBuffer;
    
    // auto-correlation arrays
    std::vector<double> E, H; // tracking arrays size 8
    std::vector<double> Edown, Hdown; // detection arrays size 110
    
    // buffer management
    int inputWritePos = 0;
    int downSampleWritePos = 0;
    int downSampleCounter = 0;
    
    // tracking state
    double cyclePeriod = 100.0;  // current detected period
    int EH_offset = 0; // L value corresponding to E[0], H[0]
};

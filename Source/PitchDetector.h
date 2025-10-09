#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <limits>

class PitchDetector
{
public:
    PitchDetector() = default;
    
    // initialize buffers
    
    void prepare(double sampleRate, int maxPeriodSamples = 1760)
    {
        fs = sampleRate;
        maxPeriod = maxPeriodSamples; // maximum period length to detect, affects low frequency limit
        
        // tracking arrays store E(L) and H(L) values for current period area
        const int N = 12;
        E.resize(N, 0.0);
        H.resize(N, 0.0);
        
        // downsampled detection arrays store values for init pitch detection
        Edown.resize(150, 0.0);// energy vals (downsampled)
        Hdown.resize(150, 0.0);// correlation
        
        // circular buffers store recent audio samples for analysis
        inputBuffer.resize(maxPeriod * 8, 0.0f);
        downsampledBuffer.resize(1024, 0.0f);
        
        downsampleAccumulator = 0.0f;
        
        reset();
    }
    
    void reset()
    {
        // vals to use when starting or when tracking is lost
        detectionMode = true;
        inputWritePos = 0;
        downSampleCounter = 0;
        downSampleWritePos = 0;
        cyclePeriod = 100.0;
        EH_offset = 0; // current center period for tracking window
        downsampleAccumulator = 0.0f;
        dcOffset = 0.0f; // for highpass filter
        goodTrackingCount = 0;
        consecutiveBad = 0;
        
        std::fill(E.begin(), E.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(Edown.begin(), Edown.end(), 0.0);
        std::fill(Hdown.begin(), Hdown.end(), 0.0);
        std::fill(inputBuffer.begin(), inputBuffer.end(), 0.0f);
        std::fill(downsampledBuffer.begin(), downsampledBuffer.end(), 0.0f);
    }
    
    bool processSample(float sample)
    {
        // apply high-pass filter to avoid false low frequency detecting
        float filteredSample = highPassFilter(sample);
        
        // store in circular buffer using modulo for wrap around
        inputBuffer[inputWritePos] = filteredSample;
        inputWritePos = (inputWritePos + 1) % inputBuffer.size();
        
        
        if (detectionMode)
        {
            // downsampling for detection
            // reduce computational load by searching wide range after downsampling
            downsampleAccumulator += filteredSample;
            downSampleCounter++;
            
            // downsample by a factor of 4
            if (downSampleCounter >= 4)
            {
                // average 4 samples to create one downsampled sample & simple antialiasing filter
                downsampledBuffer[downSampleWritePos] = downsampleAccumulator / 4.0f;
                downSampleWritePos = (downSampleWritePos + 1) % downsampledBuffer.size();
                
                // reset for next downsampling block
                downsampleAccumulator = 0.0f;
                downSampleCounter = 0;
                
                if (downSampleWritePos > 300)
                {
                    
                    // try detection periodically to balance cpu load vs responsiveness
                    if (downSampleWritePos % 4 == 0)
                    {
                        return detectPitch();
                    }
                }
            }
        }
        else
        {
            return trackPitch();
        }
        
        return false;
    }
    
    // NORMALIZED DIFFERENCE FUNCTION is the core pitch detection algorithm
    // which measures similarity between segmentds of signal separated by period L
    // returns 0 for perfect match and 1 for no correlation
    /*
     * Let x[n] be the signal, then we compare segments:
     segment1: x[t-2L] to x[t-2L+L-1]
     segment2: x[t-L] to x[t-L+L-1]
     the normalized correlation = sum(x1*x2) / sqrt(sum(x1^2) * sum(x2^2))
     
     
     */
    double computeDifference(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0) return 1.0; // invalid period
        
        double sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
        int bufSize = buffer.size();
        
        // compute for 2 consecutive segments of length L
        for (int k = 0; k < L; k++)
        {
            // calculate buffer indices with wrap around using modulo
            // the +10 ensures we stay positive before modulo
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            float x = buffer[idx1];
            float y = buffer[idx2];
            
            // accumulate stats for normalized coreelation
            sumXX += x * x; // energy of first segment
            sumYY += y * y; // energy of second segment
            sumXY += x * y; // cross correlation between segments
        }
        
        // Normalized Difference Function
        // by geometric mean
        double denominator = std::sqrt(sumXX * sumYY);
        if (denominator < 1e-10) return 1.0;
        
        // pearson correlation coefficient between the two segments
        double correlation = sumXY / denominator;
        return 1.0 - correlation;  // 0 = perfect correlation, 1 = no correlation
    }
    /*
     ENERGY COMPUTATION: calculate E(L) = average energy of two consecutive segments
     used in tracking mode with the E-2H minimization
     
     E(L) = (1/L) * Σ [x₁²(k) + x₂²(k)] for k=0 to L-1
     where x₁ = segment [t-2L, t-2L+L-1], x₂ = segment [t-L, t-1]
     
     L is period candidate length in samples
     returns normalized energy measure
     */
    double computeE(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0) return 0.0;
        
        double sum = 0.0;
        int bufSize = buffer.size();
        
        for (int k = 0; k < L; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            float x1 = buffer[idx1];
            float x2 = buffer[idx2];
            
            // sum of squares gives energy of both segments
            sum += x1 * x1 + x2 * x2;
        }
        return sum / L; //normalize by length to compare across different L
    }
    
    /*
     CORRELATION COMPUTATION:  calculate H(L) = average cross-correlation between segments
     used in tracking mode with the E-2H minimization approach
     
     H(L) = (1/L) * Σ [x₁(k) * x₂(k)] for k=0 to L-1
     measures how similar the two segments are
     
     L is period candidate length in samples
     returns normalized correlation measure
     */
    double computeH(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0) return 0.0;
        
        double sum = 0.0;
        int bufSize = buffer.size();
        
        for (int k = 0; k < L; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            // cross correlation measures similarity between segments
            sum += buffer[idx1] * buffer[idx2];
        }
        return sum / L;
    }
    
    /*
     PITCH DETECTION using normalized difference function
     searches for period L that minimizes the difference between consecutive segments
     includes octave correction to avoid harmonic errors
     
     steps:
     1. compute difference function for all candidate periods
     2. find period with minimum difference (best correlation)
     3. apply octave correction to find true fundamental
     4. validate frequency range and initialize tracking
     
     returns true if pitch is successfully detectd, false otherwise
     */
    
    bool detectPitch()
    {
        // storage for difference function vals
        std::vector<double> differences(110, 1.0);
        double minDifference = 1.0;
        int bestL = -1;
        
        //compute difference function for candidate periods
        for (int L = 20; L < 100; L++)  // Focus on reasonable range, covers 176-880Hz at 44.1
        {
            differences[L] = computeDifference(L, downsampledBuffer, downSampleWritePos);
            //get best
            if (differences[L] < minDifference)
            {
                minDifference = differences[L];
                bestL = L;
            }
        }
        
        // check if periodic signal is found
        // threshold at 0.3 means 70% correlation between segments is good
        if (bestL > 0 && minDifference < 0.3)
        {
            // octave error correction: the detected period might be a multiple of the true fundamental
            int fundamentalL = bestL;
            double fundamentalDiff = minDifference;
            
            // Test if this could be a harmonic of a shorter period
            for (int divisor = 2; divisor <= 6; divisor++)
            {
                int testL = bestL / divisor;
                if (testL >= 15)  // reasonable minimum period
                {
                    double testDiff = computeDifference(testL, downsampledBuffer, downSampleWritePos);
                    
                    // if divisor period has similar correlation, it might be fundamental
                    // the *1.2 allows slightly worse correlation for shorter periods
                    if (testDiff < fundamentalDiff * 1.2)
                    {
                        fundamentalL = testL;
                        fundamentalDiff = testDiff;
                        DBG("Octave correction: " << bestL << " -> " << testL << " (divisor " << divisor << ")");
                    }
                }
            }
            
            // convert downsampled period to full sample rate period
            // we downsampled by 4, so multiply by 4
            int fullRatePeriod = 4 * fundamentalL;
            double detectedFreq = fs / fullRatePeriod;
            
            DBG("Detection: L=" << fundamentalL << ", period=" << fullRatePeriod << ", freq=" << detectedFreq << "Hz, diff=" << fundamentalDiff);
            
            // validate that detected frequency is in reasonable vocal range
            if (detectedFreq < 100.0 || detectedFreq > 900.0)
            {
                DBG("Frequency out of range: " << detectedFreq << "Hz");
                return false;
            }
            
            // initialize tracking with the found period
            initializeTracking(fullRatePeriod);
            return true;
        }
        
        return false;
    }
    
    /*
     initialize tracking sets up tracing arrays around detected period
     creates a window of E(L) and H(L) values centered on the detected period for tracking small period variations
     */
    void initializeTracking(int period)
    {
        const int N = E.size();
        
        // center the tracking window around the detected period
        // EH_offset represents the smallest period in our tracking window
        
        EH_offset = period - N / 2;
        if (EH_offset < 20) EH_offset = 20;
        
        //pre compute E and H vals for all periods in tracking window
        // this gives us local landscape to find a precise minimum
        
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
        }
        
        // set init state for tracking mode
        cyclePeriod = static_cast<double>(period);
        detectionMode = false;
        goodTrackingCount = 10;  // start with some confidence
        consecutiveBad = 0;
        
        double freq = fs / period;
        DBG("Tracking initialized: " << freq << "Hz, period=" << period);
    }
    
    /*
     PITCH TRACKING: refine and follow pitch using E-2H minimization
     once pitch is detected, this method tracks small variations using local search
     
     E(L) measures total energy, H(L) measures correlation
     when segments are perfectly periodic: E(L) ≈ 2H(L), so E(L)-2H(L) ≈ 0
     */
    bool trackPitch()
    {
        const int N = E.size();
        
        //update E and H arrays with current audio data
        // this maintains a sliding window of period candidates
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
        }
        
        // Find minimum
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
        
        if (minIdx < 0)
        {
            DBG("Tracking: no minimum found");
            return switchToDetection();
        }
        
        // TRACKING VALIDATION: chack if we still have good periodicity
        // ratio = (E-2H)/E measures how close we are to perfect periodicity
        // perfect periodicity: E ≈ 2H, so ratio ≈ 0
        if (E[minIdx] > 1e-8)
        {
            double ratio = minVal / E[minIdx];
            
            if (ratio < 0.3)  // good periodicity (at least 70% correlation)
            {
                goodTrackingCount++;
                consecutiveBad = 0;
            }
            else // bad periodicity
            {
                consecutiveBad++;
                goodTrackingCount = std::max(0, goodTrackingCount - 1);
                
                //only switch back to detection after many consecutive failures
                //this provides hysteresis to prevent mode bouncing
                if (consecutiveBad > 35)
                {
                    DBG("Tracking lost after " << consecutiveBad << " bad frames");
                    return switchToDetection();
                }
            }
        }
        
        // adaptive window shifting
        if (minIdx <= 1 && EH_offset > 30)
        {
            shiftArraysLeft();
        }
        else if (minIdx >= N - 2 && EH_offset < 600)
        {
            shiftArraysRight();
        }
        
        // subsample interpolation
        // find precise peeriod between integer samples
        double Pmin = interpolateMinimum(minIdx);
        double newPeriod = EH_offset + Pmin;
        
        // validate period change, prevent unrealistic jumps
        if (cyclePeriod > 10.0)
        {
            double ratio = newPeriod / cyclePeriod;
            if (ratio > 1.5 || ratio < 0.67)  // Allow ±50% change
            {
                DBG("Period jump too large: " << cyclePeriod << " -> " << newPeriod);
                return true;  // Keep old period
            }
        }
        
        cyclePeriod = newPeriod;
        return true;
    }
    
    void shiftArraysLeft()
    {
        const int N = E.size();
        for (int i = N - 1; i > 0; i--)
        {
            E[i] = E[i - 1];
            H[i] = H[i - 1];
        }
        EH_offset--;
        int L = EH_offset;
        E[0] = computeE(L, inputBuffer, inputWritePos);
        H[0] = computeH(L, inputBuffer, inputWritePos);
    }
    
    void shiftArraysRight()
    {
        const int N = E.size();
        for (int i = 0; i < N - 1; i++)
        {
            E[i] = E[i + 1];
            H[i] = H[i + 1];
        }
        EH_offset++;
        int L = EH_offset + N - 1;
        E[N - 1] = computeE(L, inputBuffer, inputWritePos);
        H[N - 1] = computeH(L, inputBuffer, inputWritePos);
    }
    
    bool switchToDetection()
    {
        detectionMode = true;
        goodTrackingCount = 0;
        consecutiveBad = 0;
        return false;
    }
    
    // PARABOLIC INTERPOLATION: find subsample minimum for higher precision. fits a parabola through three points and finds the precise minimum location
    /*
     Given three points (x₀,y₀), (x₁,y₁), (x₂,y₂) where x₁=0, x₂=1
          * The parabola equation: y = ax² + bx + c
          * The minimum occurs at x = -b/(2a)
     
     returns interpolated subsample position of the minimum
     */
    double interpolateMinimum(int idx)
    {
        const int N = E.size();
        
        if (idx <= 0 || idx >= N - 1)
            return static_cast<double>(idx);
        
        double y0 = E[idx - 1] - 2.0 * H[idx - 1];
        double y1 = E[idx] - 2.0 * H[idx];
        double y2 = E[idx + 1] - 2.0 * H[idx + 1];
        
        double denom = 2.0 * (y0 - 2.0 * y1 + y2);
        
        if (std::abs(denom) < 1e-10)
            return static_cast<double>(idx);
        
        double offset = (y0 - y2) / denom;
        offset = std::max(-0.5, std::min(0.5, offset));
        
        return idx + offset;
    }
    
    float highPassFilter(float sample)
    {
        float filtered = sample - dcOffset;
        dcOffset = sample + 0.995f * (dcOffset - sample);
        return filtered;
    }
    
    double getCyclePeriod() const { return cyclePeriod; }
    
    double getFrequency() const
    {
        return (cyclePeriod > 0.0) ? (fs / cyclePeriod) : 0.0;
    }
    
    bool isInDetectionMode() const { return detectionMode; }
    
private:
    double fs = 44100.0;
    int maxPeriod = 1760;
    bool detectionMode = true;
    
    std::vector<float> inputBuffer;
    std::vector<float> downsampledBuffer;
    std::vector<double> E, H;
    std::vector<double> Edown, Hdown;
    
    int inputWritePos = 0;
    int downSampleWritePos = 0;
    int downSampleCounter = 0;
    float downsampleAccumulator = 0.0f;
    
    double cyclePeriod = 100.0;
    int EH_offset = 0;
    
    float dcOffset = 0.0f;
    int goodTrackingCount = 0;
    int consecutiveBad = 0;
};

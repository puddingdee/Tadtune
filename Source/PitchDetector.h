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
    void prepare(double sampleRate, int maxPeriodSamples = 3520)  // Double the max period for low frequencies
    {
        fs = sampleRate;
        maxPeriod = maxPeriodSamples; // Now supports ~12.5Hz at 44.1kHz
        
        // tracking arrays store E(L) and H(L) values for current period area
        const int N = 16;  // Larger tracking window for low frequencies
        E.resize(N, 0.0);
        H.resize(N, 0.0);
        
        // downsampled detection arrays store values for init pitch detection
        Edown.resize(200, 0.0);  // Larger for low frequency detection
        Hdown.resize(200, 0.0);
        
        // LARGER circular buffers for low frequency detection
        inputBuffer.resize(maxPeriod * 8, 0.0f);  // Now ~0.64 seconds at 44.1kHz
        downsampledBuffer.resize(2048, 0.0f);     // Larger downsampled buffer
        
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
        detectionCounter = 0;
        trackingCounter = 0;
        framesSinceLastDetection = 0;
        
        std::fill(E.begin(), E.end(), 0.0);
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(Edown.begin(), Edown.end(), 0.0);
        std::fill(Hdown.begin(), Hdown.end(), 0.0);
        std::fill(inputBuffer.begin(), inputBuffer.end(), 0.0f);
        std::fill(downsampledBuffer.begin(), downsampledBuffer.end(), 0.0f);
    }
    
    bool processSample(float sample)
    {
        // Apply high-pass filter (gentler for low frequencies)
        float filteredSample = highPassFilter(sample);
        
        // Store in circular buffer
        inputBuffer[inputWritePos] = filteredSample;
        inputWritePos = (inputWritePos + 1) % inputBuffer.size();
        
        framesSinceLastDetection++;
        
        if (detectionMode)
        {
            // Try detection periodically
            detectionCounter++;
            if (detectionCounter >= 128)  // Less frequent to save CPU
            {
                detectionCounter = 0;
                bool detected = detectPitch();
                if (detected) {
                    framesSinceLastDetection = 0;
                }
                return detected;
            }
        }
        else
        {
            // In tracking mode, update tracking less frequently
            trackingCounter++;
            if (trackingCounter >= 32)
            {
                trackingCounter = 0;
                bool stillTracking = trackPitch();
                if (!stillTracking) {
                    detectionMode = true;
                }
                return stillTracking;
            }
            
            // Force re-detection after a while
            if (framesSinceLastDetection > 22050) { // ~0.5 second
                detectionMode = true;
                framesSinceLastDetection = 0;
            }
        }
        
        return false;
    }
    
    // NORMALIZED DIFFERENCE FUNCTION
    // which measures similarity between segments of signal separated by period L
    // returns 0 for perfect match and 1 for no correlation
    double computeDifference(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0 || L > 200) return 1.0; // Reasonable limits
        
        double sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
        int bufSize = buffer.size();
        
        // Only compute for a reasonable number of samples to reduce CPU
        int computeSamples = std::min(L, 64); // Limit computation
        
        for (int k = 0; k < computeSamples; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            float x = buffer[idx1];
            float y = buffer[idx2];
            
            // accumulate stats for normalized correlation
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

    double computeE(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0 || L > 400) return 0.0;
        
        double sum = 0.0;
        int bufSize = buffer.size();
        int computeSamples = std::min(L, 64); // Limit computation
        
        for (int k = 0; k < computeSamples; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            float x1 = buffer[idx1];
            float x2 = buffer[idx2];
            
            // sum of squares gives energy of both segments
            sum += x1 * x1 + x2 * x2;
        }
        return sum / computeSamples; //normalize by length to compare across different L
    }
    
    double computeH(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0 || L > 400) return 0.0;
        
        double sum = 0.0;
        int bufSize = buffer.size();
        int computeSamples = std::min(L, 64); // Limit computation
        
        for (int k = 0; k < computeSamples; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            // cross correlation measures similarity between segments
            sum += buffer[idx1] * buffer[idx2];
        }
        return sum / computeSamples;
    }
    
    bool detectPitch()
    {
        // We need enough data in the input buffer
        if (inputWritePos < 400) return false;
        
        double minDifference = 1.0;
        int bestL = -1;
        
        // Search wider range for low frequencies using NON-downsampled data
        // 80Hz at 44.1kHz = ~551 samples, 1000Hz = ~44 samples
        for (int L = 40; L < 600; L += 2)  // Wider range, step by 2 for speed
        {
            // Skip if we don't have enough buffer history
            if (inputWritePos < L * 2) continue;
            
            double diff = computeDifference(L, inputBuffer, inputWritePos);
            
            if (diff < minDifference)
            {
                minDifference = diff;
                bestL = L;
            }
        }
        
        
        // More permissive threshold for low frequencies
        if (bestL > 0 && minDifference < 0.5)
        {
            double detectedFreq = fs / bestL;
            
            // Wider frequency range
            if (detectedFreq >= 70.0 && detectedFreq <= 1200.0)
            {
                initializeTracking(bestL);
                return true;
            }
            else
            {
                DBG("Frequency out of range: " << detectedFreq << " Hz");
            }
        }
        
        return false;
    }

    void initializeTracking(int period)
    {
        const int N = E.size();
        
        // center the tracking window around the detected period
        EH_offset = period - N / 2;
        if (EH_offset < 20) EH_offset = 20;
        
        //pre compute E and H vals for all periods in tracking window
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            if (L >= 20 && L <= 500)
            {
                E[i] = computeE(L, inputBuffer, inputWritePos);
                H[i] = computeH(L, inputBuffer, inputWritePos);
            }
        }
        
        // set init state for tracking mode
        cyclePeriod = static_cast<double>(period);
        detectionMode = false;
        goodTrackingCount = 10;  // start with some confidence
        consecutiveBad = 0;
        
        double freq = fs / period;
    }
    
    bool trackPitch()
    {
        const int N = E.size();
        
        // Update E and H arrays with current audio data
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            // Only update if L is reasonable
            if (L >= 20 && L <= 500)
            {
                E[i] = computeE(L, inputBuffer, inputWritePos);
                H[i] = computeH(L, inputBuffer, inputWritePos);
            }
        }
        
        // Find minimum of E-2H
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
            consecutiveBad++;
            if (consecutiveBad > 5) {
                return false; // Signal tracking failure
            }
            return true; // Give it a few more tries
        }
        
        // Validate the minimum
        if (E[minIdx] > 1e-8)
        {
            double ratio = minVal / E[minIdx];
            
            if (ratio < 0.4)  // Good periodicity
            {
                goodTrackingCount++;
                consecutiveBad = 0;
                
                // Update period with interpolation
                double Pmin = interpolateMinimum(minIdx);
                double newPeriod = EH_offset + Pmin;
                double newFreq = fs / newPeriod;
                
                // Validate period change
                if (cyclePeriod > 10.0)
                {
                    double periodRatio = newPeriod / cyclePeriod;
                    if (periodRatio > 1.5 || periodRatio < 0.67)
                    {
                        consecutiveBad++;
                        if (consecutiveBad > 3) {
                            return false; // Signal tracking failure
                        }
                    }
                    else
                    {
                        cyclePeriod = newPeriod;
                        consecutiveBad = 0;
                        
                        // Log significant frequency changes
                        static double lastLoggedFreq = 0.0;
                        if (std::abs(newFreq - lastLoggedFreq) > 5.0) {
                            lastLoggedFreq = newFreq;
                        }
                    }
                }
                else
                {
                    cyclePeriod = newPeriod;
                }
                
                // Adaptive window shifting
                if (minIdx <= 1 && EH_offset > 30)
                {
                    shiftArraysLeft();
                }
                else if (minIdx >= N - 2 && EH_offset < 600)
                {
                    shiftArraysRight();
                }
                
                return true;
            }
            else
            {
                consecutiveBad++;
                if (consecutiveBad > 8)
                {
                    return false; // Signal tracking failure
                }
            }
        }
        else
        {
            consecutiveBad++;
            if (consecutiveBad > 5) {
                return false;
            }
        }
        
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
    int detectionCounter = 0;
    int trackingCounter = 0;
    int framesSinceLastDetection = 0;
};

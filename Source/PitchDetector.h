#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <limits>

class PitchDetector
{
public:
    PitchDetector() = default;
    
    void prepare(double sampleRate, int maxPeriodSamples = 1760)
    {
        fs = sampleRate;
        maxPeriod = maxPeriodSamples;
        
        const int N = 12;
        E.resize(N, 0.0);
        H.resize(N, 0.0);
        
        Edown.resize(150, 0.0);
        Hdown.resize(150, 0.0);
        
        inputBuffer.resize(maxPeriod * 8, 0.0f);
        downsampledBuffer.resize(1024, 0.0f);
        
        downsampleAccumulator = 0.0f;
        
        reset();
    }
    
    void reset()
    {
        detectionMode = true;
        inputWritePos = 0;
        downSampleCounter = 0;
        downSampleWritePos = 0;
        cyclePeriod = 100.0;
        EH_offset = 0;
        downsampleAccumulator = 0.0f;
        dcOffset = 0.0f;
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
        // Apply high-pass filter
        float filteredSample = highPassFilter(sample);
        
        inputBuffer[inputWritePos] = filteredSample;
        inputWritePos = (inputWritePos + 1) % inputBuffer.size();
        
        if (detectionMode)
        {
            downsampleAccumulator += filteredSample;
            downSampleCounter++;
            
            if (downSampleCounter >= 4)
            {
                downsampledBuffer[downSampleWritePos] = downsampleAccumulator / 4.0f;
                downSampleWritePos = (downSampleWritePos + 1) % downsampledBuffer.size();
                
                downsampleAccumulator = 0.0f;
                downSampleCounter = 0;
                
                if (downSampleWritePos > 300)
                {
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
    
    /**
     * IMPROVED: Compute normalized difference function (NDF)
     * This is more robust than the E-2H approach
     */
    double computeDifference(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0) return 1.0;
        
        double sumXX = 0.0, sumYY = 0.0, sumXY = 0.0;
        int bufSize = buffer.size();
        
        for (int k = 0; k < L; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            float x = buffer[idx1];
            float y = buffer[idx2];
            
            sumXX += x * x;
            sumYY += y * y;
            sumXY += x * y;
        }
        
        // Normalized Difference Function
        double denominator = std::sqrt(sumXX * sumYY);
        if (denominator < 1e-10) return 1.0;
        
        double correlation = sumXY / denominator;
        return 1.0 - correlation;  // 0 = perfect correlation, 1 = no correlation
    }
    
    /**
     * Backward compatibility - keep original for tracking
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
            
            sum += x1 * x1 + x2 * x2;
        }
        return sum / L;
    }
    
    double computeH(int L, const std::vector<float>& buffer, int writePos)
    {
        if (L <= 0) return 0.0;
        
        double sum = 0.0;
        int bufSize = buffer.size();
        
        for (int k = 0; k < L; k++)
        {
            int idx1 = (writePos - 2 * L + k + 10 * bufSize) % bufSize;
            int idx2 = (writePos - L + k + 10 * bufSize) % bufSize;
            
            sum += buffer[idx1] * buffer[idx2];
        }
        return sum / L;
    }
    
    bool detectPitch()
    {
        // Use NDF for more robust detection
        std::vector<double> differences(110, 1.0);
        double minDifference = 1.0;
        int bestL = -1;
        
        // Compute normalized difference function
        for (int L = 20; L < 100; L++)  // Focus on reasonable range
        {
            differences[L] = computeDifference(L, downsampledBuffer, downSampleWritePos);
            if (differences[L] < minDifference)
            {
                minDifference = differences[L];
                bestL = L;
            }
        }
        
        // Check if we found a good candidate
        if (bestL > 0 && minDifference < 0.3)  // Good correlation threshold
        {
            // CRITICAL FIX: Check for octave errors by testing integer divisors
            int fundamentalL = bestL;
            double fundamentalDiff = minDifference;
            
            // Test if this could be a harmonic of a shorter period
            for (int divisor = 2; divisor <= 6; divisor++)
            {
                int testL = bestL / divisor;
                if (testL >= 15)  // Reasonable minimum period
                {
                    double testDiff = computeDifference(testL, downsampledBuffer, downSampleWritePos);
                    
                    // If the divisor period has similar or better correlation, use it
                    if (testDiff < fundamentalDiff * 1.2)  // Allow slightly worse
                    {
                        fundamentalL = testL;
                        fundamentalDiff = testDiff;
                        DBG("Octave correction: " << bestL << " -> " << testL << " (divisor " << divisor << ")");
                    }
                }
            }
            
            // Convert to full sample rate
            int fullRatePeriod = 4 * fundamentalL;
            double detectedFreq = fs / fullRatePeriod;
            
            DBG("Detection: L=" << fundamentalL << ", period=" << fullRatePeriod << ", freq=" << detectedFreq << "Hz, diff=" << fundamentalDiff);
            
            // Validate frequency range
            if (detectedFreq < 100.0 || detectedFreq > 900.0)
            {
                DBG("Frequency out of range: " << detectedFreq << "Hz");
                return false;
            }
            
            // Initialize tracking with the found period
            initializeTracking(fullRatePeriod);
            return true;
        }
        
        return false;
    }
    
    void initializeTracking(int period)
    {
        const int N = E.size();
        EH_offset = period - N / 2;
        if (EH_offset < 20) EH_offset = 20;
        
        for (int i = 0; i < N; i++)
        {
            int L = EH_offset + i;
            E[i] = computeE(L, inputBuffer, inputWritePos);
            H[i] = computeH(L, inputBuffer, inputWritePos);
        }
        
        cyclePeriod = static_cast<double>(period);
        detectionMode = false;
        goodTrackingCount = 10;  // Start with some confidence
        consecutiveBad = 0;
        
        double freq = fs / period;
        DBG("Tracking initialized: " << freq << "Hz, period=" << period);
    }
    
    bool trackPitch()
    {
        const int N = E.size();
        
        // Update E and H arrays
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
        
        // SIMPLIFIED VALIDATION
        if (E[minIdx] > 1e-8)
        {
            double ratio = minVal / E[minIdx];
            
            if (ratio < 0.3)  // Good frame
            {
                goodTrackingCount++;
                consecutiveBad = 0;
            }
            else  // Bad frame
            {
                consecutiveBad++;
                goodTrackingCount = std::max(0, goodTrackingCount - 1);
                
                if (consecutiveBad > 30)  // Much more tolerant
                {
                    DBG("Tracking lost after " << consecutiveBad << " bad frames");
                    return switchToDetection();
                }
            }
        }
        
        // Shift arrays if needed (more conservative)
        if (minIdx <= 1 && EH_offset > 30)
        {
            shiftArraysLeft();
        }
        else if (minIdx >= N - 2 && EH_offset < 600)
        {
            shiftArraysRight();
        }
        
        // Update period with validation
        double Pmin = interpolateMinimum(minIdx);
        double newPeriod = EH_offset + Pmin;
        
        // Validate period change
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

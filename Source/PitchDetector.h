#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <algorithm>

/**
 * Optimized YIN Pitch Detection Algorithm with Median Filtering
 * CPU-efficient version for real-time vocal pitch detection
 */
class PitchDetector
{
public:
    PitchDetector() = default;
    
    void prepare(double sampleRate, int bufferSize = 1024)
    {
        fs = sampleRate;
        yinBufferSize = bufferSize;
        
        // Allocate buffers
        audioBuffer.resize(yinBufferSize * 2, 0.0f);
        yinBuffer.resize(yinBufferSize / 2, 0.0f);
        
        // Voice range: 80Hz to 800Hz
        minPeriod = static_cast<int>(fs / 800.0);
        maxPeriod = static_cast<int>(fs / 80.0);
        maxPeriod = std::min(maxPeriod, yinBufferSize / 2);
        
        hopSize = 512;
        
        // Initialize pitch history for median filtering
        pitchHistory.resize(5, 0.0f);  // Last 5 pitch estimates
        historyIndex = 0;
        
        reset();
    }
    
    void reset()
    {
        writePos = 0;
        samplesInBuffer = 0;
        frameCounter = 0;
        
        currentFrequency = 0.0f;
        confidence = 0.0f;
        
        std::fill(audioBuffer.begin(), audioBuffer.end(), 0.0f);
        std::fill(yinBuffer.begin(), yinBuffer.end(), 0.0f);
        std::fill(pitchHistory.begin(), pitchHistory.end(), 0.0f);
        
        smoothedFrequency = 0.0f;
        smoothedConfidence = 0.0f;
        lastGoodFrequency = 0.0f;
        framesWithoutPitch = 0;
        historyIndex = 0;
        stableFrameCount = 0;
        
        prevInputSample = 0.0f;
        prevFilteredSample = 0.0f;
    }
    
    bool processSample(float sample)
    {
        // Apply high-pass filter
        sample = highPassFilter(sample);
        
        // Store in ring buffer
        audioBuffer[writePos] = sample;
        writePos = (writePos + 1) % audioBuffer.size();
        samplesInBuffer = std::min(samplesInBuffer + 1, (int)audioBuffer.size());
        
        frameCounter++;
        
        // Process every hopSize samples
        if (frameCounter >= hopSize)
        {
            frameCounter = 0;
            
            if (samplesInBuffer >= yinBufferSize)
            {
                return detectPitch();
            }
        }
        
        return false;
    }
    
    // Getters
    double getFrequency() const { return smoothedFrequency; }
    float getConfidence() const { return smoothedConfidence; }
    bool isInDetectionMode() const { return confidence < 0.15f; }
    bool hasValidPitch() const { return smoothedFrequency > 80.0f && smoothedFrequency < 800.0f && smoothedConfidence > 0.1f; }
    
private:
    bool detectPitch()
    {
        // Quick energy check
        if (!hasEnoughEnergy())
        {
            handleNoPitchDetected();
            return false;
        }
        
        // YIN algorithm
        calculateDifferenceFast();
        cumulativeMeanNormalizedDifference();
        
        int tauEstimate = absoluteThreshold();
        
        if (tauEstimate == -1)
        {
            handleNoPitchDetected();
            return false;
        }
        
        // Parabolic interpolation
        float betterTau = parabolicInterpolation(tauEstimate);
        float detectedFreq = fs / betterTau;
        float currentConfidence = 1.0f - yinBuffer[tauEstimate];
        
        // Validate frequency range
        if (detectedFreq < 80.0f || detectedFreq > 800.0f)
        {
            handleNoPitchDetected();
            return false;
        }
        
        // Octave error correction
        detectedFreq = correctOctaveErrors(detectedFreq);
        
        // Add to pitch history for median filtering
        pitchHistory[historyIndex] = detectedFreq;
        historyIndex = (historyIndex + 1) % pitchHistory.size();
        
        // Get median-filtered pitch (removes outliers/spikes)
        float medianFreq = getMedianPitch();
        
        // Use median if we have enough history, otherwise use detected
        float finalFreq = (stableFrameCount > 3) ? medianFreq : detectedFreq;
        
        // Additional validation against previous estimate
        if (lastGoodFrequency > 0.0f && currentConfidence < 0.6f)
        {
            float ratio = finalFreq / lastGoodFrequency;
            
            // Reject if too far from previous (unless very confident)
            if (ratio > 1.15f || ratio < 0.87f)
            {
                // Use previous frequency instead of rejecting entirely
                finalFreq = lastGoodFrequency;
                currentConfidence *= 0.7f;  // Reduce confidence
            }
        }
        
        // Update with smoothing
        updatePitchEstimate(finalFreq, currentConfidence);
        
        return true;
    }
    
    // Correct common octave errors
    float correctOctaveErrors(float detectedFreq)
    {
        if (lastGoodFrequency > 0.0f)
        {
            float ratio = detectedFreq / lastGoodFrequency;
            
            // Fix octave up
            if (ratio > 1.85f && ratio < 2.15f)
                return detectedFreq * 0.5f;
            
            // Fix octave down
            if (ratio > 0.47f && ratio < 0.54f)
                return detectedFreq * 2.0f;
            
            // Fix fifth up (3:2 ratio)
            if (ratio > 1.45f && ratio < 1.55f)
                return detectedFreq * (2.0f / 3.0f);
            
            // Fix fifth down (2:3 ratio)
            if (ratio > 0.65f && ratio < 0.70f)
                return detectedFreq * (3.0f / 2.0f);
        }
        
        return detectedFreq;
    }
    
    // Median filter to remove pitch spikes/outliers
    float getMedianPitch()
    {
        // Copy history and sort
        std::vector<float> sorted = pitchHistory;
        std::sort(sorted.begin(), sorted.end());
        
        // Remove zeros (unfilled history)
        sorted.erase(std::remove_if(sorted.begin(), sorted.end(),
                                    [](float f) { return f < 80.0f; }),
                     sorted.end());
        
        if (sorted.empty())
            return 0.0f;
        
        // Return median
        size_t mid = sorted.size() / 2;
        if (sorted.size() % 2 == 0)
            return (sorted[mid - 1] + sorted[mid]) / 2.0f;
        else
            return sorted[mid];
    }
    
    bool hasEnoughEnergy()
    {
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        float sumSquares = 0.0f;
        
        for (int i = 0; i < yinBufferSize; i += 4)
        {
            int idx = (startPos + i) % audioBuffer.size();
            float val = audioBuffer[idx];
            sumSquares += val * val;
        }
        
        float rms = std::sqrt(sumSquares / (yinBufferSize / 4));
        return rms > 0.01f;
    }
    
    void calculateDifferenceFast()
    {
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        int searchRange = std::min((int)yinBuffer.size(), maxPeriod + 10);
        
        for (int tau = 0; tau < searchRange; ++tau)
        {
            double sum = 0.0;
            int computeLength = yinBufferSize / 3;
            
            for (int i = 0; i < computeLength; ++i)
            {
                int idx1 = (startPos + i) % audioBuffer.size();
                int idx2 = (startPos + i + tau) % audioBuffer.size();
                
                float delta = audioBuffer[idx1] - audioBuffer[idx2];
                sum += delta * delta;
            }
            
            yinBuffer[tau] = sum;
        }
        
        for (int tau = searchRange; tau < yinBuffer.size(); ++tau)
        {
            yinBuffer[tau] = 1e10f;
        }
    }
    
    void cumulativeMeanNormalizedDifference()
    {
        yinBuffer[0] = 1.0f;
        double runningSum = 0.0;
        int searchRange = std::min((int)yinBuffer.size(), maxPeriod + 10);
        
        for (int tau = 1; tau < searchRange; ++tau)
        {
            runningSum += yinBuffer[tau];
            yinBuffer[tau] = (runningSum == 0.0) ? 1.0f : yinBuffer[tau] * tau / runningSum;
        }
    }
    
    int absoluteThreshold()
    {
        float threshold = 0.15f;
        int tau = minPeriod;
        
        while (tau < maxPeriod && yinBuffer[tau] >= threshold)
            tau++;
        
        if (tau >= maxPeriod)
            return -1;
        
        int minTau = tau;
        float minVal = yinBuffer[tau];
        
        while (tau < maxPeriod && yinBuffer[tau] < threshold)
        {
            if (yinBuffer[tau] < minVal)
            {
                minVal = yinBuffer[tau];
                minTau = tau;
            }
            tau++;
        }
        
        return (minVal > 0.5f) ? -1 : minTau;
    }
    
    float parabolicInterpolation(int tauEstimate)
    {
        if (tauEstimate < 1 || tauEstimate >= maxPeriod - 1)
            return static_cast<float>(tauEstimate);
        
        float s0 = yinBuffer[tauEstimate - 1];
        float s1 = yinBuffer[tauEstimate];
        float s2 = yinBuffer[tauEstimate + 1];
        
        float denom = 2.0f * (2.0f * s1 - s2 - s0);
        if (std::abs(denom) < 1e-10f)
            return static_cast<float>(tauEstimate);
        
        float adjustment = (s2 - s0) / denom;
        adjustment = std::max(-0.5f, std::min(0.5f, adjustment));
        
        return tauEstimate + adjustment;
    }
    
    void updatePitchEstimate(float newFreq, float newConfidence)
    {
        currentFrequency = newFreq;
        confidence = newConfidence;
        
        // Count stable frames
        if (lastGoodFrequency > 0.0f)
        {
            float ratio = newFreq / lastGoodFrequency;
            if (ratio > 0.98f && ratio < 1.02f)  // Within 2%
                stableFrameCount++;
            else
                stableFrameCount = 0;
        }
        
        // Confidence-based smoothing
        float smoothingFactor;
        
        if (newConfidence > 0.7f && stableFrameCount > 5)
        {
            // High confidence and stable: light smoothing
            smoothingFactor = 0.15f;
        }
        else if (newConfidence > 0.5f)
        {
            // Medium confidence: moderate smoothing
            smoothingFactor = 0.35f;
        }
        else
        {
            // Low confidence: heavy smoothing
            smoothingFactor = 0.6f;
        }
        
        if (smoothedFrequency == 0.0f)
        {
            smoothedFrequency = newFreq;
            smoothedConfidence = newConfidence;
        }
        else
        {
            smoothedFrequency = smoothingFactor * smoothedFrequency + (1.0f - smoothingFactor) * newFreq;
            smoothedConfidence = 0.6f * smoothedConfidence + 0.4f * newConfidence;
        }
        
        lastGoodFrequency = smoothedFrequency;
        framesWithoutPitch = 0;
    }
    
    void handleNoPitchDetected()
    {
        confidence = 0.0f;
        framesWithoutPitch++;
        stableFrameCount = 0;
        
        if (framesWithoutPitch > 3)
        {
            smoothedConfidence *= 0.7f;
            
            if (smoothedConfidence < 0.05f)
            {
                smoothedFrequency = 0.0f;
                smoothedConfidence = 0.0f;
                lastGoodFrequency = 0.0f;
                std::fill(pitchHistory.begin(), pitchHistory.end(), 0.0f);
            }
        }
    }
    
    float highPassFilter(float input)
    {
        const float alpha = 0.995f;
        float output = alpha * (prevFilteredSample + input - prevInputSample);
        prevInputSample = input;
        prevFilteredSample = output;
        return output;
    }
    
    // Member variables
    double fs = 44100.0;
    int yinBufferSize = 1024;
    int minPeriod = 55;
    int maxPeriod = 551;
    int hopSize = 512;
    
    std::vector<float> audioBuffer;
    std::vector<float> yinBuffer;
    
    int writePos = 0;
    int samplesInBuffer = 0;
    int frameCounter = 0;
    
    float currentFrequency = 0.0f;
    float confidence = 0.0f;
    float smoothedFrequency = 0.0f;
    float smoothedConfidence = 0.0f;
    float lastGoodFrequency = 0.0f;
    int framesWithoutPitch = 0;
    int stableFrameCount = 0;
    
    // Median filtering for pitch stability
    std::vector<float> pitchHistory;
    int historyIndex = 0;
    
    float prevInputSample = 0.0f;
    float prevFilteredSample = 0.0f;
};

#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <algorithm>

/**
 * Optimized YIN Pitch Detection Algorithm
 * CPU-efficient version for real-time vocal pitch detection
 */
class PitchDetector
{
public:
    PitchDetector() = default;
    
    void prepare(double sampleRate, int bufferSize = 1024)  // Smaller buffer
    {
        fs = sampleRate;
        yinBufferSize = bufferSize;
        
        // Allocate buffers
        audioBuffer.resize(yinBufferSize * 2, 0.0f);
        yinBuffer.resize(yinBufferSize / 2, 0.0f);
        
        // Voice range: 80Hz to 800Hz (more restricted)
        minPeriod = static_cast<int>(fs / 800.0);   // ~55 samples at 44.1kHz
        maxPeriod = static_cast<int>(fs / 80.0);    // ~551 samples at 44.1kHz
        maxPeriod = std::min(maxPeriod, yinBufferSize / 2);
        
        // Larger hop for efficiency
        hopSize = 512;  // Process less frequently (~11.6ms at 44.1kHz)
        
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
        
        smoothedFrequency = 0.0f;
        smoothedConfidence = 0.0f;
        lastGoodFrequency = 0.0f;
        framesWithoutPitch = 0;
        
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
        // Quick energy check to skip silent passages
        if (!hasEnoughEnergy())
        {
            handleNoPitchDetected();
            return false;
        }
        
        // YIN algorithm steps
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
        if (lastGoodFrequency > 0.0f)
        {
            float ratio = detectedFreq / lastGoodFrequency;
            
            // Fix common octave errors
            if (ratio > 1.85f && ratio < 2.15f)
                detectedFreq *= 0.5f;
            else if (ratio > 0.47f && ratio < 0.54f)
                detectedFreq *= 2.0f;
            
            // Reject if still too far from previous (unless high confidence)
            ratio = detectedFreq / lastGoodFrequency;
            if ((ratio > 1.2f || ratio < 0.83f) && currentConfidence < 0.6f)
            {
                handleNoPitchDetected();
                return false;
            }
        }
        
        // Update with smoothing
        updatePitchEstimate(detectedFreq, currentConfidence);
        
        return true;
    }
    
    // Fast energy check
    bool hasEnoughEnergy()
    {
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        float sumSquares = 0.0f;
        
        // Check every 4th sample for speed
        for (int i = 0; i < yinBufferSize; i += 4)
        {
            int idx = (startPos + i) % audioBuffer.size();
            float val = audioBuffer[idx];
            sumSquares += val * val;
        }
        
        float rms = std::sqrt(sumSquares / (yinBufferSize / 4));
        return rms > 0.01f;  // Noise gate
    }
    
    // Optimized difference calculation
    void calculateDifferenceFast()
    {
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        
        // Only calculate up to maxPeriod (not full buffer)
        int searchRange = std::min((int)yinBuffer.size(), maxPeriod + 10);
        
        for (int tau = 0; tau < searchRange; ++tau)
        {
            double sum = 0.0;
            
            // Reduced computation window for speed
            int computeLength = yinBufferSize / 3;  // Only use 1/3 of buffer
            
            for (int i = 0; i < computeLength; ++i)
            {
                int idx1 = (startPos + i) % audioBuffer.size();
                int idx2 = (startPos + i + tau) % audioBuffer.size();
                
                float delta = audioBuffer[idx1] - audioBuffer[idx2];
                sum += delta * delta;
            }
            
            yinBuffer[tau] = sum;
        }
        
        // Fill rest with high values (won't be selected)
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
        
        // Find first point below threshold
        while (tau < maxPeriod && yinBuffer[tau] >= threshold)
            tau++;
        
        if (tau >= maxPeriod)
            return -1;
        
        // Find minimum
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
        
        // Adaptive smoothing
        float smoothingFactor = 0.4f;
        if (newConfidence > 0.7f)
            smoothingFactor = 0.2f;  // Less smoothing when confident
        
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
        
        if (framesWithoutPitch > 3)
        {
            smoothedConfidence *= 0.7f;
            
            if (smoothedConfidence < 0.05f)
            {
                smoothedFrequency = 0.0f;
                smoothedConfidence = 0.0f;
                lastGoodFrequency = 0.0f;
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
    
    float prevInputSample = 0.0f;
    float prevFilteredSample = 0.0f;
};

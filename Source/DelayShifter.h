#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>

/**
 * Delay-Based Pitch Shifter
 * Uses variable delay with crossfading for pitch correction
 * Much simpler and lower latency than phase vocoder
 */
class DelayBasedPitchShifter
{
public:
    DelayBasedPitchShifter() = default;
    
    void prepare(double sampleRate, int maxDelaySamples = 4096)
    {
        fs = sampleRate;
        
        // Allocate delay buffer (circular buffer)
        delayBuffer.resize(maxDelaySamples, 0.0f);
        
        // Crossfade buffers for smooth transitions
        crossfadeLength = 256;  // ~5.8ms at 44.1kHz
        
        reset();
    }
    
    void reset()
    {
        writePos = 0;
        readPos1 = 0.0f;
        readPos2 = crossfadeLength / 2.0f;
        
        currentPitchRatio = 1.0f;
        targetPitchRatio = 1.0f;
        
        crossfadePhase = 0.0f;
        
        std::fill(delayBuffer.begin(), delayBuffer.end(), 0.0f);
    }
    
    void setPitchShiftRatio(float ratio)
    {
        // Clamp to reasonable range
        targetPitchRatio = std::max(0.5f, std::min(2.0f, ratio));
    }
    
    float processSample(float inputSample)
    {
        // Write input to delay buffer
        delayBuffer[writePos] = inputSample;
        
        // Smooth the pitch ratio change
        currentPitchRatio += (targetPitchRatio - currentPitchRatio) * 0.01f;
        
        // Calculate playback speed based on pitch ratio
        // ratio > 1.0 = pitch up = read faster
        // ratio < 1.0 = pitch down = read slower
        float playbackSpeed = currentPitchRatio;
        
        // Read from two positions with crossfading
        float output1 = readInterpolated(readPos1);
        float output2 = readInterpolated(readPos2);
        
        // Advance read positions
        readPos1 += playbackSpeed;
        readPos2 += playbackSpeed;
        
        // Wrap read positions
        while (readPos1 >= delayBuffer.size())
            readPos1 -= delayBuffer.size();
        while (readPos2 >= delayBuffer.size())
            readPos2 -= delayBuffer.size();
        
        // Reset read head when it gets too close to write head
        float distance1 = getDistanceFromWrite(readPos1);
        float distance2 = getDistanceFromWrite(readPos2);
        
        // If a read head gets too close, reset it far away
        float minDistance = crossfadeLength * 4;
        float resetDistance = delayBuffer.size() * 0.5f;
        
        if (distance1 < minDistance)
        {
            readPos1 = writePos - resetDistance;
            if (readPos1 < 0)
                readPos1 += delayBuffer.size();
        }
        
        if (distance2 < minDistance)
        {
            readPos2 = writePos - resetDistance;
            if (readPos2 < 0)
                readPos2 += delayBuffer.size();
        }
        
        // Crossfade between the two read heads
        crossfadePhase += 1.0f / crossfadeLength;
        if (crossfadePhase >= 1.0f)
            crossfadePhase -= 1.0f;
        
        // Use equal-power crossfade
        float gain1 = std::cos(crossfadePhase * juce::MathConstants<float>::halfPi);
        float gain2 = std::sin(crossfadePhase * juce::MathConstants<float>::halfPi);
        
        float output = output1 * gain1 + output2 * gain2;
        
        // Advance write position
        writePos = (writePos + 1) % delayBuffer.size();
        
        return output;
    }
    
private:
    // Linear interpolation for fractional delay reading
    float readInterpolated(float position)
    {
        int index = static_cast<int>(position);
        float frac = position - index;
        
        int index1 = index % delayBuffer.size();
        int index2 = (index + 1) % delayBuffer.size();
        
        return delayBuffer[index1] * (1.0f - frac) + delayBuffer[index2] * frac;
    }
    
    // Calculate distance from read position to write position
    float getDistanceFromWrite(float readPosition)
    {
        int readIdx = static_cast<int>(readPosition);
        int distance = writePos - readIdx;
        
        if (distance < 0)
            distance += delayBuffer.size();
        
        return static_cast<float>(distance);
    }
    
    // Member variables
    double fs = 44100.0;
    std::vector<float> delayBuffer;
    
    int writePos = 0;
    float readPos1 = 0.0f;    // First read head
    float readPos2 = 0.0f;    // Second read head (for crossfading)
    
    float currentPitchRatio = 1.0f;
    float targetPitchRatio = 1.0f;
    
    float crossfadePhase = 0.0f;
    int crossfadeLength = 256;
};

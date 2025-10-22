#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>

/**
 * Simple Overlap-Add Time Stretcher with Resampling
 * Based on SOLA (Synchronized Overlap-Add) principles
 */
class DelayBasedPitchShifter
{
public:
    DelayBasedPitchShifter() = default;
    
    void prepare(double sampleRate, int maxDelaySamples = 8192)
    {
        fs = sampleRate;
        
        // Input buffer
        inputBuffer.resize(maxDelaySamples, 0.0f);
        bufferSize = maxDelaySamples;
        
        // Output FIFO
        outputFifo.resize(maxDelaySamples * 2, 0.0f);
        
        // Grain parameters
        grainSize = 1024;  // ~23ms at 44.1kHz
        hopSize = 512;     // 50% overlap
        
        // Hann window
        window.resize(grainSize);
        for (int i = 0; i < grainSize; ++i)
        {
            window[i] = 0.5f * (1.0f - std::cos(2.0f * juce::MathConstants<float>::pi * i / grainSize));
        }
        
        reset();
    }
    
    void reset()
    {
        writePos = 0;
        readPos = 0.0f;
        outputReadPos = 0;
        outputWritePos = grainSize;  // Start ahead
        
        samplesUntilNextGrain = hopSize;
        
        currentPitchRatio = 1.0f;
        targetPitchRatio = 1.0f;
        
        std::fill(inputBuffer.begin(), inputBuffer.end(), 0.0f);
        std::fill(outputFifo.begin(), outputFifo.end(), 0.0f);
    }
    
    void setPitchShiftRatio(float ratio)
    {
        targetPitchRatio = std::max(0.7f, std::min(1.4f, ratio));
    }
    
    float processSample(float inputSample)
    {
        // Write input to buffer
        inputBuffer[writePos] = inputSample;
        writePos = (writePos + 1) % bufferSize;
        
        // Smooth ratio changes
        currentPitchRatio += (targetPitchRatio - currentPitchRatio) * 0.005f;
        
        // Check if it's time to process a new grain
        samplesUntilNextGrain--;
        if (samplesUntilNextGrain <= 0)
        {
            processGrain();
            
            // Calculate next hop based on pitch ratio
            // For pitch up, we want shorter input hops (read faster)
            // For pitch down, we want longer input hops (read slower)
            float inputHop = hopSize / currentPitchRatio;
            samplesUntilNextGrain = hopSize;  // Output hop is constant
            readPos += inputHop;
            
            // Wrap read position
            while (readPos >= bufferSize)
                readPos -= bufferSize;
            while (readPos < 0)
                readPos += bufferSize;
        }
        
        // Read from output FIFO
        float output = outputFifo[outputReadPos];
        outputFifo[outputReadPos] = 0.0f;  // Clear after reading
        outputReadPos = (outputReadPos + 1) % outputFifo.size();
        
        return output;
    }
    
private:
    void processGrain()
    {
        // Extract grain from input buffer
        int readStart = static_cast<int>(readPos);
        
        // Apply window and write to output with overlap-add
        for (int i = 0; i < grainSize; ++i)
        {
            int inputIdx = (readStart + i) % bufferSize;
            float sample = inputBuffer[inputIdx] * window[i];
            
            // Overlap-add to output
            int outputIdx = (outputWritePos + i) % outputFifo.size();
            outputFifo[outputIdx] += sample;
        }
        
        // Advance output write position by constant hop
        outputWritePos = (outputWritePos + hopSize) % outputFifo.size();
    }
    
    // Member variables
    double fs = 44100.0;
    int bufferSize = 8192;
    int grainSize = 1024;
    int hopSize = 512;
    
    std::vector<float> inputBuffer;
    std::vector<float> outputFifo;
    std::vector<float> window;
    
    int writePos = 0;
    float readPos = 0.0f;
    int outputReadPos = 0;
    int outputWritePos = 0;
    
    int samplesUntilNextGrain = 0;
    
    float currentPitchRatio = 1.0f;
    float targetPitchRatio = 1.0f;
};

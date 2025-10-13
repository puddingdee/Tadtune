#pragma once
#include <JuceHeader.h>
#include "PitchDetector.h"

class TadtuneAudioProcessor : public juce::AudioProcessor
{
public:
    //==========================================================================
    TadtuneAudioProcessor()
        : AudioProcessor(BusesProperties()
            .withInput("Input", juce::AudioChannelSet::stereo(), true)
            .withOutput("Output", juce::AudioChannelSet::stereo(), true))
    {
        juce::RuntimePermissions::request(
                juce::RuntimePermissions::recordAudio,
                [](bool granted)
                {
                    DBG("Microphone permission: " << (granted ? "GRANTED" : "DENIED"));
                }
            );
    }

    ~TadtuneAudioProcessor() override {}

    //==========================================================================
    // PREPARATION
    //==========================================================================
    void prepareToPlay(double sampleRate, int samplesPerBlock) override
    {
        // Initialize the pitch detector for this sample rate
        pitchDetector.prepare(sampleRate);
        
        // Initialize pitch correction with reasonable buffer sizes
        currentSampleRate = sampleRate;
        resampleRate = 1.0;
        smoothedResampleRate = 1.0;
        outputAddress = 0.0;
        
        // Larger buffer for better cycle management (but still reasonable latency)
        inputBuffer.resize(8192); // ~186ms at 44.1kHz for better low frequency handling
        std::fill(inputBuffer.begin(), inputBuffer.end(), 0.0f);
        inputWritePos = 0;
        
        // Reset detection state
        currentFrequency.store(0.0f);
        isDetecting.store(true);
        targetFrequency.store(0.0f);
        lastValidPeriod = 100.0;
        
        // Initialize crossfade buffers for click-free transitions
        crossfadeBuffer.resize(256, 0.0f);
        crossfadeActive = false;
        crossfadePosition = 0;
        
        // Initialize period change smoothing
        lastPeriod = 100.0;
        periodSmoothingCoeff = 0.95f;
    }

    void releaseResources() override
    {
        // Clean up when playback stops
    }

    //==========================================================================
    // AUDIO PROCESSING
    //==========================================================================
    void processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) override
    {
        juce::ScopedNoDenormals noDenormals;
        
        auto totalNumInputChannels  = getTotalNumInputChannels();
        auto totalNumOutputChannels = getTotalNumOutputChannels();
        
        // Clear any output channels that don't contain input data
        for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
            buffer.clear (i, 0, buffer.getNumSamples());

        // Get the left channel (we'll process mono for pitch detection)
        auto* channelData = buffer.getReadPointer(0);
        int numSamples = buffer.getNumSamples();
        
        // Process each sample through the pitch detector
        for (int i = 0; i < numSamples; ++i)
        {
            float sample = channelData[i];
            
            // Store in input buffer for pitch correction
            inputBuffer[inputWritePos] = sample;
            inputWritePos = (inputWritePos + 1) % inputBuffer.size();
            
            // Feed sample to detector
            pitchDetector.processSample(sample);
        }
        
        // Update pitch correction parameters ONCE per block
        float freq = static_cast<float>(pitchDetector.getFrequency());
        currentFrequency.store(freq);
        float rawPeriod = static_cast<float>(pitchDetector.getCyclePeriod());
        isDetecting.store(pitchDetector.isInDetectionMode());
        
        // Calculate target frequency and resample rate if we have valid pitch
        if (freq > 70.0f && freq < 1200.0f && !pitchDetector.isInDetectionMode())
        {
            float targetFreq = findNearestChromaticFrequency(freq);
            targetFrequency.store(targetFreq);
            
            float targetPeriod = currentSampleRate / targetFreq;
            
            // Smooth period changes to reduce artifacts
            float smoothedPeriod = lastPeriod * periodSmoothingCoeff + rawPeriod * (1.0f - periodSmoothingCoeff);
            lastPeriod = smoothedPeriod;
            
            if (smoothedPeriod > 30.0f && smoothedPeriod < 600.0f)
            {
                lastValidPeriod = smoothedPeriod;
                currentPeriod.store(smoothedPeriod);
                
                // Calculate new resample rate
                float newResampleRate = smoothedPeriod / targetPeriod;
                
                // TIGHTER clamping to prevent extreme values that cause distortion
                newResampleRate = std::max(0.75f, std::min(1.35f, newResampleRate));
                
                // Detect significant rate changes and trigger crossfade
                if (std::abs(newResampleRate - resampleRate) > 0.05f && !crossfadeActive)
                {
                    // Start crossfade to prevent clicks
                    startCrossfade();
                }
                
                resampleRate = newResampleRate;
                // Apply gentle smoothing to resample rate for instant but click-free correction
                smoothedResampleRate = smoothedResampleRate * 0.9f + resampleRate * 0.1f;
                
                DBG("Correction - Current: " << freq << " Target: " << targetFreq << " Rate: " << smoothedResampleRate);
            }
            else
            {
                smoothedResampleRate = 1.0f;
            }
        }
        else
        {
            // No pitch detected or in detection mode - pass through
            smoothedResampleRate = 1.0f;
            targetFrequency.store(0.0f);
        }
        
        // Apply pitch correction to the output
        applyPitchCorrection(buffer);
    }

    //==========================================================================
    // PITCH CORRECTION METHODS
    //==========================================================================
    
    /** Find nearest chromatic frequency */
    float findNearestChromaticFrequency(float inputFreq)
    {
        if (inputFreq <= 0.0f) return inputFreq;
        
        // Convert to MIDI note number
        float midiNote = 69.0f + 12.0f * std::log2(inputFreq / 440.0f);
        int roundedNote = std::round(midiNote);
        
        // Convert back to frequency
        return 440.0f * std::pow(2.0f, (roundedNote - 69.0f) / 12.0f);
    }
    
    /** Start crossfade to prevent clicks when pitch changes */
    void startCrossfade()
    {
        // Capture current output for crossfading
        int captureLength = std::min(256, (int)inputBuffer.size());
        for (int i = 0; i < captureLength; ++i)
        {
            crossfadeBuffer[i] = cubicInterpolate(outputAddress);
            outputAddress += smoothedResampleRate;
            if (outputAddress >= inputBuffer.size())
                outputAddress -= inputBuffer.size();
        }
        
        crossfadeActive = true;
        crossfadePosition = 0;
        
        // Reset output address to current input position for phase continuity
        outputAddress = static_cast<float>(inputWritePos) - lastValidPeriod * 2.0f;
        while (outputAddress < 0) outputAddress += inputBuffer.size();
        while (outputAddress >= inputBuffer.size()) outputAddress -= inputBuffer.size();
    }
    
    /** Apply pitch correction using time-domain resampling with cycle alignment */
    void applyPitchCorrection(juce::AudioBuffer<float>& buffer)
    {
        int numSamples = buffer.getNumSamples();
        int bufferSize = inputBuffer.size();
        
        // If we're in detection mode or no valid pitch, just pass through
        if (pitchDetector.isInDetectionMode() || std::abs(smoothedResampleRate - 1.0f) < 0.01f)
        {
            // Simple pass-through
            for (int channel = 0; channel < buffer.getNumChannels(); ++channel)
            {
                if (channel == 0) continue;
                buffer.copyFrom(channel, 0, buffer, 0, 0, numSamples);
            }
            return;
        }
        
        // Get current cycle period for pitch correction
        double currentPeriod = pitchDetector.getCyclePeriod();
        if (currentPeriod < 30.0 || currentPeriod > 600.0)
        {
            // Invalid period - pass through
            smoothedResampleRate = 1.0f;
            return;
        }
        
        // Process each sample with click-free pitch correction
        for (int i = 0; i < numSamples; ++i)
        {
            float outputSample;
            
            // Handle crossfade if active
            if (crossfadeActive)
            {
                float newSample = cubicInterpolate(outputAddress);
                float oldSample = crossfadeBuffer[crossfadePosition];
                
                // Linear crossfade
                float fadeAmount = static_cast<float>(crossfadePosition) / 256.0f;
                outputSample = oldSample * (1.0f - fadeAmount) + newSample * fadeAmount;
                
                crossfadePosition++;
                if (crossfadePosition >= 256)
                {
                    crossfadeActive = false;
                    crossfadePosition = 0;
                }
            }
            else
            {
                // Normal interpolation with cubic for smoother result
                outputSample = cubicInterpolate(outputAddress);
            }
            
            // Write to all output channels with CLAMPING to prevent clipping
            for (int channel = 0; channel < buffer.getNumChannels(); ++channel)
            {
                // Clamp to prevent distortion
                float clampedSample = std::max(-0.99f, std::min(0.99f, outputSample));
                buffer.getWritePointer(channel)[i] = clampedSample;
            }
            
            // Move output pointer
            outputAddress += smoothedResampleRate;
            
            // Wrap buffer address with proper cycle alignment
            while (outputAddress >= bufferSize)
            {
                outputAddress -= bufferSize;
            }
            while (outputAddress < 0.0f)
            {
                outputAddress += bufferSize;
            }
            
            // Maintain proper distance from write pointer to prevent buffer overrun/underrun
            float distanceFromWrite = inputWritePos - outputAddress;
            if (distanceFromWrite < 0) distanceFromWrite += bufferSize;
            
            // Target distance is 2 periods behind write head
            float targetDistance = lastValidPeriod * 2.0f;
            float maxDistance = bufferSize * 0.8f;
            float minDistance = lastValidPeriod * 1.5f;
            
            // Gradual correction if we drift too far
            if (distanceFromWrite > maxDistance)
            {
                // We're too far behind, speed up slightly
                outputAddress += 0.5f;
            }
            else if (distanceFromWrite < minDistance)
            {
                // We're too close, slow down slightly
                outputAddress -= 0.5f;
            }
        }
    }
    
    /** Cubic interpolation for smoother pitch shifting without aliasing */
    float cubicInterpolate(float address)
    {
        int bufferSize = inputBuffer.size();
        
        int idx1 = static_cast<int>(std::floor(address));
        float frac = address - idx1;
        
        // Get 4 points for cubic interpolation
        int idx0 = (idx1 - 1 + bufferSize) % bufferSize;
        idx1 = idx1 % bufferSize;
        int idx2 = (idx1 + 1) % bufferSize;
        int idx3 = (idx1 + 2) % bufferSize;
        
        float y0 = inputBuffer[idx0];
        float y1 = inputBuffer[idx1];
        float y2 = inputBuffer[idx2];
        float y3 = inputBuffer[idx3];
        
        // 4-point cubic interpolation (Catmull-Rom)
        float c0 = y1;
        float c1 = 0.5f * (y2 - y0);
        float c2 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
        float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
        
        return ((c3 * frac + c2) * frac + c1) * frac + c0;
    }

    //==========================================================================
    // QUERY METHODS FOR UI
    //==========================================================================
    
    /** Get the currently detected frequency in Hz */
    float getCurrentFrequency() const
    {
        return currentFrequency.load();
    }
    
    /** Get the target corrected frequency in Hz */
    float getTargetFrequency() const
    {
        return targetFrequency.load();
    }
    
    /** Get the current period in samples */
    float getCurrentPeriod() const
    {
        return currentPeriod.load();
    }
    
    /** Check if in detection mode vs tracking mode */
    bool isInDetectionMode() const
    {
        return isDetecting.load();
    }
    
    /** Get current resample rate */
    float getResampleRate() const
    {
        return smoothedResampleRate;
    }
    
    /** Convert frequency to nearest MIDI note number */
    int frequencyToMidiNote(float freq) const
    {
        if (freq <= 0.0f) return -1;
        return static_cast<int>(std::round(69.0f + 12.0f * std::log2(freq / 440.0f)));
    }
    
    /** Get note name from MIDI note number */
    juce::String midiNoteToName(int midiNote) const
    {
        if (midiNote < 0) return "---";
        
        const char* noteNames[] = { "C", "C#", "D", "D#", "E", "F",
                                   "F#", "G", "G#", "A", "A#", "B" };
        int octave = (midiNote / 12) - 1;
        int note = midiNote % 12;
        
        return juce::String(noteNames[note]) + juce::String(octave);
    }

    //==========================================================================
    // BOILERPLATE (required by JUCE)
    //==========================================================================
    
    const juce::String getName() const override
    {
        return "Tadtune";
    }

    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.0; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return {}; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock& destData) override
    {
        // Save plugin state (empty for now)
    }

    void setStateInformation(const void* data, int sizeInBytes) override
    {
        // Restore plugin state (empty for now)
    }

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

private:
    // The pitch detector instance
    PitchDetector pitchDetector;
    
    // Pitch correction variables
    double currentSampleRate = 44100.0;
    float resampleRate = 1.0f;
    float smoothedResampleRate = 1.0f;
    float outputAddress = 0.0f;
    float lastValidPeriod = 100.0f;
    float lastPeriod = 100.0f;
    float periodSmoothingCoeff = 0.95f;
    
    // Input buffer for pitch correction
    std::vector<float> inputBuffer;
    int inputWritePos = 0;
    
    // Crossfade buffers to eliminate clicks on pitch changes
    std::vector<float> crossfadeBuffer;
    bool crossfadeActive = false;
    int crossfadePosition = 0;
    
    // Thread-safe variables for communicating with UI
    std::atomic<float> currentFrequency { 0.0f };
    std::atomic<float> targetFrequency { 0.0f };
    std::atomic<float> currentPeriod { 0.0f };
    std::atomic<bool> isDetecting { true };
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessor)
};

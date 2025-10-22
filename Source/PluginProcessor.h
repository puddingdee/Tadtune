#pragma once
#include <JuceHeader.h>
#include "PitchDetector.h"
#include "PhaseVocoder.h"

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
        
        // Initialize phase vocoder with smaller FFT for lower latency (faster snap)
        phaseVocoder.prepare(sampleRate, 1024);
        
        currentSampleRate = sampleRate;
        
        // Reset detection state
        currentFrequency.store(0.0f);
        isDetecting.store(true);
        targetFrequency.store(0.0f);
        currentPitchRatio.store(1.0f);
        
        // Smoothing for pitch ratio changes
        smoothedPitchRatio = 1.0f;
        pitchRatioSmoothingCoeff = 0.01;
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
            buffer.clear(i, 0, buffer.getNumSamples());

        // Get the left channel (we'll process mono for pitch detection)
        auto* channelData = buffer.getReadPointer(0);
        int numSamples = buffer.getNumSamples();
        
        // Process each sample through the pitch detector
        for (int i = 0; i < numSamples; ++i)
        {
            float sample = channelData[i];
            pitchDetector.processSample(sample);
        }
        
        // Update pitch information ONCE per block
        float detectedFreq = static_cast<float>(pitchDetector.getFrequency());
        currentFrequency.store(detectedFreq);
        isDetecting.store(pitchDetector.isInDetectionMode());
        
        // Calculate target frequency and pitch shift ratio
        float targetPitchRatio = 1.0f;
        
        if (detectedFreq > 70.0f && detectedFreq < 1200.0f && !pitchDetector.isInDetectionMode())
        {
            float targetFreq = findNearestChromaticFrequency(detectedFreq);
            targetFrequency.store(targetFreq);
            
            // Calculate pitch shift ratio
            // Ratio > 1.0 means shift up, < 1.0 means shift down
            targetPitchRatio = targetFreq / detectedFreq;
            
            // Clamp to reasonable range
            targetPitchRatio = std::max(0.75f, std::min(1.35f, targetPitchRatio));
            
        }
        else
        {
            // No valid pitch detected - bypass
            targetFrequency.store(0.0f);
            targetPitchRatio = 1.0f;
        }
        
        // Smooth the pitch ratio to avoid artifacts
        smoothedPitchRatio = smoothedPitchRatio * pitchRatioSmoothingCoeff
                           + targetPitchRatio * (1.0f - pitchRatioSmoothingCoeff);
        
        currentPitchRatio.store(smoothedPitchRatio);
        
        // Update phase vocoder with new pitch shift ratio
        phaseVocoder.setPitchShiftRatio(smoothedPitchRatio);
        
        // Process audio through phase vocoder
        applyPhaseVocoderCorrection(buffer);
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
    
    float findNearestScaleFrequency(float inputFreq)
    {
        float midiNote = 69.0f + 12.0f * std::log2(inputFreq / 440.0f);
        int roundedNote = std::round(midiNote);
        return 440.f * std::pow(2.0f, (roundedNote - 69.0f) / 12.f);
    }
    
    /** Apply pitch correction using phase vocoder */
    void applyPhaseVocoderCorrection(juce::AudioBuffer<float>& buffer)
    {
        int numSamples = buffer.getNumSamples();
        int numChannels = buffer.getNumChannels();
        
        // Process left channel through phase vocoder
        auto* leftChannel = buffer.getWritePointer(0);
        
        for (int i = 0; i < numSamples; ++i)
        {
            float input = leftChannel[i];
            float output = phaseVocoder.processSample(input);
            
            // Soft clipping to prevent harsh distortion
            output = std::tanh(output * 0.9f);
            
            leftChannel[i] = output;
        }
        
        // Copy left to other channels
        for (int channel = 1; channel < numChannels; ++channel)
        {
            buffer.copyFrom(channel, 0, buffer, 0, 0, numSamples);
        }
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
    
    /** Check if in detection mode vs tracking mode */
    bool isInDetectionMode() const
    {
        return isDetecting.load();
    }
    
    /** Get current pitch shift ratio */
    float getPitchRatio() const
    {
        return currentPitchRatio.load();
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
    
    // Phase vocoder for pitch correction
    PhaseVocoder phaseVocoder;
    
    // Audio parameters
    double currentSampleRate = 44100.0;
    float smoothedPitchRatio = 1.0f;
    float pitchRatioSmoothingCoeff = 0.98f;
    
    
    // Thread-safe variables for communicating with UI
    std::atomic<float> currentFrequency { 0.0f };
    std::atomic<float> targetFrequency { 0.0f };
    std::atomic<float> currentPitchRatio { 1.0f };
    std::atomic<bool> isDetecting { true };
    
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessor)
};

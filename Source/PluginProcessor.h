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
        
        // Reset detection state
        currentFrequency.store(0.0f);
        isDetecting.store(true);
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
        
        // Get the left channel (we'll process mono for pitch detection)
        auto* channelData = buffer.getReadPointer(0);
        int numSamples = buffer.getNumSamples();
        
        // Process each sample through the pitch detector
        for (int i = 0; i < numSamples; ++i)
        {
            float sample = channelData[i];
            
            // Feed sample to detector
            bool detected = pitchDetector.processSample(sample);
            
            // Update our thread-safe state variables periodically
            // (don't do it every sample to reduce atomic overhead)
            if (i % 64 == 0)
            {
                currentFrequency.store(static_cast<float>(pitchDetector.getFrequency()));
                currentPeriod.store(static_cast<float>(pitchDetector.getCyclePeriod()));
                isDetecting.store(pitchDetector.isInDetectionMode());
            }
        }
        
        // For testing: pass audio through unchanged
        // (Later you'd add pitch correction here)
        auto* outputLeft = buffer.getWritePointer(0);
        auto* outputRight = buffer.getWritePointer(1);
        
        // Copy left to right if stereo
        if (buffer.getNumChannels() > 1)
        {
            buffer.copyFrom(1, 0, buffer, 0, 0, numSamples);
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
        return "PitchDetector";
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
    
    // Thread-safe variables for communicating with UI
    // (audio thread writes, UI thread reads)
    std::atomic<float> currentFrequency { 0.0f };
    std::atomic<float> currentPeriod { 0.0f };
    std::atomic<bool> isDetecting { true };
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessor)
};

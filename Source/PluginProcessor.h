#pragma once
#include <JuceHeader.h>
#include "PitchDetector.h"
#include "PhaseVocoder.h"

class TadtuneAudioProcessorEditor;

class TadtuneAudioProcessor : public juce::AudioProcessor
{
public:
    // Scale types (relative major/minor pairs)
    enum class ScaleType
    {
        Chromatic = 0,
        C_Major_A_Minor,
        Db_Major_Bb_Minor,
        D_Major_B_Minor,
        Eb_Major_C_Minor,
        E_Major_Cs_Minor,
        F_Major_D_Minor,
        Gb_Major_Eb_Minor,
        G_Major_E_Minor,
        Ab_Major_F_Minor,
        A_Major_Fs_Minor,
        Bb_Major_G_Minor,
        B_Major_Gs_Minor
    };
    
    //==========================================================================
    TadtuneAudioProcessor()
        : AudioProcessor(BusesProperties()
            .withInput("Input", juce::AudioChannelSet::stereo(), true)
            .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
          parameters(*this, nullptr, "PARAMETERS", createParameterLayout())
    {
        smoothingParam = parameters.getRawParameterValue("smoothing");
        wetGain = parameters.getRawParameterValue("wet gain");
        dryGain = parameters.getRawParameterValue("dry gain");
        scaleTypeParam = parameters.getRawParameterValue("scale");
        
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
        pitchDetector.prepare(sampleRate);
        pitchShifter.prepare(sampleRate, 2048);
        pitchShifter.setSmoothingTime(smoothingParam->load());
        
        currentSampleRate = sampleRate;
        
        currentFrequency.store(0.0f);
        isDetecting.store(true);
        targetFrequency.store(0.0f);
        currentPitchRatio.store(1.0f);
    }

    void releaseResources() override
    {
    }

    //==========================================================================
    // AUDIO PROCESSING
    //==========================================================================
    void processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) override
    {
        juce::ScopedNoDenormals noDenormals;
        
        pitchShifter.setSmoothingTime(smoothingParam->load());
        
        auto totalNumInputChannels  = getTotalNumInputChannels();
        auto totalNumOutputChannels = getTotalNumOutputChannels();
        
        for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
            buffer.clear(i, 0, buffer.getNumSamples());

        auto* channelData = buffer.getReadPointer(0);
        int numSamples = buffer.getNumSamples();
        
        for (int i = 0; i < numSamples; ++i)
        {
            pitchDetector.processSample(channelData[i]);
        }
        
        float detectedFreq = static_cast<float>(pitchDetector.getFrequency());
        float confidence = pitchDetector.getConfidence();
        
        currentFrequency.store(detectedFreq);
        isDetecting.store(pitchDetector.isInDetectionMode());
        
        float pitchRatio = 1.0f;
        
        if (detectedFreq > 80.0f && detectedFreq < 800.0f && confidence > 0.15f)
        {
            // Get current scale setting
            ScaleType scale = static_cast<ScaleType>(static_cast<int>(scaleTypeParam->load()));
            
            float targetFreq = findNearestScaleFrequency(detectedFreq, scale);
            targetFrequency.store(targetFreq);
            pitchRatio = targetFreq / detectedFreq;
            pitchRatio = std::max(0.75f, std::min(1.35f, pitchRatio));
        }
        else
        {
            targetFrequency.store(0.0f);
            pitchRatio = 1.0f;
        }
        
        currentPitchRatio.store(pitchRatio);
        pitchShifter.setPitchShiftRatio(pitchRatio);
        
        applyPitchCorrectionWithMix(buffer);
    }

    //==========================================================================
    // PARAMETER LAYOUT
    //==========================================================================
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout()
    {
        std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;
        
        params.push_back(std::make_unique<juce::AudioParameterFloat>(
            "smoothing",
            "Smoothing Time",
            juce::NormalisableRange<float>(0.0f, 200.0f, 1.0f),
            0.0f,
            juce::AudioParameterFloatAttributes()
                .withStringFromValueFunction([](float value, int) { return juce::String(value, 0) + " ms"; })
        ));
        
        params.push_back(std::make_unique<juce::AudioParameterFloat>(
            "wet gain",
            "Wet Gain",
            juce::NormalisableRange<float>(0.0f, 1.0f, 0.01f),
            1.0f,
            juce::AudioParameterFloatAttributes()
                .withStringFromValueFunction([](float value, int) { return juce::String(value * 100.0f, 0) + " %"; })
        ));
        
        params.push_back(std::make_unique<juce::AudioParameterFloat>(
            "dry gain",
            "Dry Gain",
            juce::NormalisableRange<float>(0.0f, 1.0f, 0.01f),
            0.0f,
            juce::AudioParameterFloatAttributes()
                .withStringFromValueFunction([](float value, int) { return juce::String(value * 100.0f, 0) + " %"; })
        ));
        
        params.push_back(std::make_unique<juce::AudioParameterChoice>(
            "scale",
            "Scale",
            juce::StringArray{
                "Chromatic",
                "C Major / A Minor",
                "Db Major / Bb Minor",
                "D Major / B Minor",
                "Eb Major / C Minor",
                "E Major / C# Minor",
                "F Major / D Minor",
                "Gb Major / Eb Minor",
                "G Major / E Minor",
                "Ab Major / F Minor",
                "A Major / F# Minor",
                "Bb Major / G Minor",
                "B Major / G# Minor"
            },
            0  // Default to Chromatic
        ));
        
        return { params.begin(), params.end() };
    }

    //==========================================================================
    // PITCH CORRECTION
    //==========================================================================
    
    /** Find nearest frequency in the selected scale */
    float findNearestScaleFrequency(float inputFreq, ScaleType scale)
    {
        if (inputFreq <= 0.0f) return inputFreq;
        
        // Convert to MIDI note number
        float midiNote = 69.0f + 12.0f * std::log2(inputFreq / 440.0f);
        
        if (scale == ScaleType::Chromatic)
        {
            // Chromatic - snap to nearest semitone
            int roundedNote = std::round(midiNote);
            return 440.0f * std::pow(2.0f, (roundedNote - 69.0f) / 12.0f);
        }
        else
        {
            // Get the scale degrees for this scale (uses major scale pattern)
            // All major scales have the same pattern: W-W-H-W-W-W-H (0, 2, 4, 5, 7, 9, 11)
            std::vector<int> scaleDegrees = {0, 2, 4, 5, 7, 9, 11};
            
            // Get root note offset based on scale selection
            int rootOffset = static_cast<int>(scale) - 1; // Subtract 1 because Chromatic is 0
            
            // Transpose scale degrees by root note
            for (auto& degree : scaleDegrees)
            {
                degree = (degree + rootOffset) % 12;
            }
            
            // Find nearest note in scale
            int baseNote = static_cast<int>(std::round(midiNote));
            int octave = baseNote / 12;
            int noteInOctave = baseNote % 12;
            
            int nearestNote = findNearestNoteInScale(noteInOctave, scaleDegrees);
            int targetMidiNote = octave * 12 + nearestNote;
            
            // Convert back to frequency
            return 440.0f * std::pow(2.0f, (targetMidiNote - 69.0f) / 12.0f);
        }
    }
    
    /** Find nearest note within an octave that belongs to the scale */
    int findNearestNoteInScale(int noteInOctave, const std::vector<int>& scaleDegrees)
    {
        int minDistance = 12;
        int nearestNote = noteInOctave;
        
        // Check current octave
        for (int degree : scaleDegrees)
        {
            int distance = std::abs(noteInOctave - degree);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestNote = degree;
            }
        }
        
        for (int degree : scaleDegrees)
        {
            int noteAbove = degree + 12;
            int distance = std::abs(noteInOctave - noteAbove);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestNote = noteAbove;
            }
        }
        
        for (int degree : scaleDegrees)
        {
            int noteBelow = degree - 12;
            int distance = std::abs(noteInOctave - noteBelow);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestNote = noteBelow;
            }
        }
        
        return nearestNote;
    }
    
    /** Apply pitch correction using phase vocoder */
    void applyPitchCorrectionWithMix(juce::AudioBuffer<float>& buffer)
    {
        int numSamples = buffer.getNumSamples();
        int numChannels = buffer.getNumChannels();
        
        float wet = wetGain->load();
        float dry = dryGain->load();
        
        juce::AudioBuffer<float> dryBuffer(numChannels, numSamples);
        
        for (int channel = 0; channel < numChannels; ++channel)
        {
            dryBuffer.copyFrom(channel, 0, buffer, channel, 0, numSamples);
        }
        
        auto* leftChannel = buffer.getWritePointer(0);
        
        for (int i = 0; i < numSamples; ++i)
        {
            float input = leftChannel[i];
            float output = pitchShifter.processSample(input);
            output = std::tanh(output * 0.95f);
            leftChannel[i] = output;
        }
        
        for (int channel = 1; channel < numChannels; ++channel)
        {
            buffer.copyFrom(channel, 0, buffer, 0, 0, numSamples);
        }
        
        for (int channel = 0; channel < numChannels; ++channel)
        {
            auto* wetSignal = buffer.getWritePointer(channel);
            auto* drySignal = dryBuffer.getReadPointer(channel);
            
            for (int i = 0; i < numSamples; ++i)
            {
                wetSignal[i] = (wetSignal[i] * wet * 6) + (drySignal[i] * dry);
            }
        }
    }
    
    //==========================================================================
    // QUERY METHODS FOR UI
    //==========================================================================
    
    float getCurrentFrequency() const { return currentFrequency.load(); }
    float getTargetFrequency() const { return targetFrequency.load(); }
    bool isInDetectionMode() const { return isDetecting.load(); }
    float getPitchRatio() const { return currentPitchRatio.load(); }
    juce::AudioProcessorValueTreeState& getParameters() { return parameters; }
    
    int frequencyToMidiNote(float freq) const
    {
        if (freq <= 0.0f) return -1;
        return static_cast<int>(std::round(69.0f + 12.0f * std::log2(freq / 440.0f)));
    }
    
    juce::String midiNoteToName(int midiNote) const
    {
        if (midiNote < 0) return "---";
        
        const char* noteNames[] = { "C", "C#", "D", "D#", "E", "F",
                                   "F#", "G", "G#", "A", "A#", "B" };
        int octave = (midiNote / 12) - 1;
        int note = midiNote % 12;
        
        return juce::String(noteNames[note]) + juce::String(octave);
    }
    
    // Get current scale for UI display
    ScaleType getCurrentScale() const
    {
        return static_cast<ScaleType>(static_cast<int>(scaleTypeParam->load()));
    }
    
    juce::String getScaleName() const
    {
        const char* scaleNames[] = {
            "Chromatic",
            "C Major / A Minor",
            "Db Major / Bb Minor",
            "D Major / B Minor",
            "Eb Major / C Minor",
            "E Major / C# Minor",
            "F Major / D Minor",
            "Gb Major / Eb Minor",
            "G Major / E Minor",
            "Ab Major / F Minor",
            "A Major / F# Minor",
            "Bb Major / G Minor",
            "B Major / G# Minor"
        };
        
        int index = static_cast<int>(getCurrentScale());
        return scaleNames[index];
    }

    //==========================================================================
    // BOILERPLATE
    //==========================================================================
    
    const juce::String getName() const override { return "Tadtune"; }
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
        auto state = parameters.copyState();
        std::unique_ptr<juce::XmlElement> xml(state.createXml());
        copyXmlToBinary(*xml, destData);
    }

    void setStateInformation(const void* data, int sizeInBytes) override
    {
        std::unique_ptr<juce::XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
        if (xmlState.get() != nullptr)
            if (xmlState->hasTagName(parameters.state.getType()))
                parameters.replaceState(juce::ValueTree::fromXml(*xmlState));
    }

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

private:
    PitchDetector pitchDetector;
    PhaseVocoder pitchShifter;
    
    juce::AudioProcessorValueTreeState parameters;
    std::atomic<float>* smoothingParam;
    std::atomic<float>* wetGain;
    std::atomic<float>* dryGain;
    std::atomic<float>* scaleTypeParam;
    
    double currentSampleRate = 44100.0;
    
    std::atomic<float> currentFrequency { 0.0f };
    std::atomic<float> targetFrequency { 0.0f };
    std::atomic<float> currentPitchRatio { 1.0f };
    std::atomic<bool> isDetecting { true };
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessor)
};

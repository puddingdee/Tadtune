    #pragma once
    #include <JuceHeader.h>
    #include "PitchDetector.h"
    #include "PhaseVocoder.h"

class TadtuneAudioProcessorEditor;

    class TadtuneAudioProcessor : public juce::AudioProcessor
    {
    public:
        //==========================================================================
        TadtuneAudioProcessor()
            : AudioProcessor(BusesProperties()
                .withInput("Input", juce::AudioChannelSet::stereo(), true)
                .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
              parameters(*this, nullptr, "PARAMETERS", createParameterLayout())
        {
            smoothingParam = parameters.getRawParameterValue("smoothing");
            
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
            // Initialize the pitch detector
            pitchDetector.prepare(sampleRate);
            
            // Initialize phase vocoder pitch shifter
            pitchShifter.prepare(sampleRate, 2048);
            
            // Set initial smoothing time
            pitchShifter.setSmoothingTime(smoothingParam->load());
            
            currentSampleRate = sampleRate;
            
            // Reset detection state
            currentFrequency.store(0.0f);
            isDetecting.store(true);
            targetFrequency.store(0.0f);
            currentPitchRatio.store(1.0f);
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
            
            // Update smoothing time if parameter changed
            pitchShifter.setSmoothingTime(smoothingParam->load());
            
            auto totalNumInputChannels  = getTotalNumInputChannels();
            auto totalNumOutputChannels = getTotalNumOutputChannels();
            
            for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
                buffer.clear(i, 0, buffer.getNumSamples());

            auto* channelData = buffer.getReadPointer(0);
            int numSamples = buffer.getNumSamples();
            
            // Process each sample through the pitch detector
            for (int i = 0; i < numSamples; ++i)
            {
                pitchDetector.processSample(channelData[i]);
            }
            
            // Update pitch information
            float detectedFreq = static_cast<float>(pitchDetector.getFrequency());
            float confidence = pitchDetector.getConfidence();
            
            currentFrequency.store(detectedFreq);
            isDetecting.store(pitchDetector.isInDetectionMode());
            
            // Calculate pitch shift ratio
            float pitchRatio = 1.0f;
            
            if (detectedFreq > 80.0f && detectedFreq < 800.0f && confidence > 0.15f)
            {
                float targetFreq = findNearestChromaticFrequency(detectedFreq);
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
            
            // Update pitch shifter with new ratio (will be smoothed internally)
            pitchShifter.setPitchShiftRatio(pitchRatio);
            
            // Process audio through pitch shifter
            applyPitchCorrection(buffer);
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
                50.0f,
                juce::AudioParameterFloatAttributes()
                    .withStringFromValueFunction([](float value, int) { return juce::String(value, 0) + " ms"; })
            ));
            
            return { params.begin(), params.end() };
        }

        //==========================================================================
        // PITCH CORRECTION
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
        
        /** Apply pitch correction using phase vocoder */
        void applyPitchCorrection(juce::AudioBuffer<float>& buffer)
        {
            int numSamples = buffer.getNumSamples();
            int numChannels = buffer.getNumChannels();
            
            // Process left channel through phase vocoder
            auto* leftChannel = buffer.getWritePointer(0);
            
            for (int i = 0; i < numSamples; ++i)
            {
                float input = leftChannel[i];
                float output = pitchShifter.processSample(input);
                
                // Soft clipping to prevent harsh distortion
                output = std::tanh(output * 0.95f);
                
                leftChannel[i] = output;
            }
            
            // Copy left to other channels (mono processing)
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
        
        /** Get APVTS reference */
        juce::AudioProcessorValueTreeState& getParameters() { return parameters; }
        
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
            // Save plugin state
            auto state = parameters.copyState();
            std::unique_ptr<juce::XmlElement> xml(state.createXml());
            copyXmlToBinary(*xml, destData);
        }

        void setStateInformation(const void* data, int sizeInBytes) override
        {
            // Restore plugin state
            std::unique_ptr<juce::XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));
            if (xmlState.get() != nullptr)
                if (xmlState->hasTagName(parameters.state.getType()))
                    parameters.replaceState(juce::ValueTree::fromXml(*xmlState));
        }

        juce::AudioProcessorEditor* createEditor() override;
        bool hasEditor() const override { return true; }

    private:
        // Core components
        PitchDetector pitchDetector;
        PhaseVocoder pitchShifter;
        
        // Parameters
        juce::AudioProcessorValueTreeState parameters;
        std::atomic<float>* smoothingParam;
        
        // Audio parameters
        double currentSampleRate = 44100.0;
        
        // Thread-safe variables for UI communication
        std::atomic<float> currentFrequency { 0.0f };
        std::atomic<float> targetFrequency { 0.0f };
        std::atomic<float> currentPitchRatio { 1.0f };
        std::atomic<bool> isDetecting { true };
        
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessor)
    };

    #pragma once
    #include <JuceHeader.h>
    #include "PluginProcessor.h"

    class TadtuneAudioProcessorEditor : public juce::AudioProcessorEditor,
                                        private juce::Timer
    {
    public:
        TadtuneAudioProcessorEditor(TadtuneAudioProcessor& p)
            : AudioProcessorEditor(&p),
              audioProcessor(p),
              smoothingAttachment(audioProcessor.getParameters(), "smoothing", smoothingSlider)
        {
            // Set window size
            setSize(400, 400);
            
            // Setup smoothing slider
            addAndMakeVisible(smoothingSlider);
            smoothingSlider.setSliderStyle(juce::Slider::LinearHorizontal);
            smoothingSlider.setTextBoxStyle(juce::Slider::TextBoxRight, false, 60, 20);
            
            addAndMakeVisible(smoothingLabel);
            smoothingLabel.setText("Smoothing Time:", juce::dontSendNotification);
            smoothingLabel.attachToComponent(&smoothingSlider, true);
            
            // Start timer to update UI at 30 Hz (every ~33ms)
            startTimerHz(30);
        }

        ~TadtuneAudioProcessorEditor() override
        {
            stopTimer();
        }

        void paint(juce::Graphics& g) override
        {
            // Background
            g.fillAll(juce::Colour(0xff1e1e1e));
            
            auto bounds = getLocalBounds();
            
            // Title
            g.setColour(juce::Colours::white);
            g.setFont(24.0f);
            g.drawText("Tadtune - Phase Vocoder 2048", bounds.removeFromTop(50),
                       juce::Justification::centred);
            
            // Get current values from processor
            float freq = audioProcessor.getCurrentFrequency();
            float targetFreq = audioProcessor.getTargetFrequency();
            bool detecting = audioProcessor.isInDetectionMode();
            float pitchRatio = audioProcessor.getPitchRatio();
            
            // Status indicator
            auto statusArea = bounds.removeFromTop(40);
            g.setFont(16.0f);
            
            if (detecting)
            {
                g.setColour(juce::Colours::orange);
                g.drawText("DETECTING...", statusArea, juce::Justification::centred);
            }
            else
            {
                g.setColour(juce::Colours::green);
                g.drawText("CORRECTING PITCH", statusArea, juce::Justification::centred);
            }
            
            bounds.removeFromTop(10);
            
            // Current frequency display
            auto freqArea = bounds.removeFromTop(60);
            drawFrequencyDisplay(g, freqArea, freq, "Current Frequency");
            
            // Target frequency display
            auto targetArea = bounds.removeFromTop(60);
            drawFrequencyDisplay(g, targetArea, targetFreq, "Target Frequency");
            
            bounds.removeFromTop(10);
            
            // Note name display
            auto noteArea = bounds.removeFromTop(60);
            drawNoteDisplay(g, noteArea, freq, targetFreq);
            
            bounds.removeFromTop(10);
            
            // Correction info
            auto infoArea = bounds.removeFromTop(40);
            drawCorrectionInfo(g, infoArea, pitchRatio);
            
            // Smoothing display
            auto smoothingArea = bounds.removeFromTop(40);
            drawSmoothingDisplay(g, smoothingArea);
        }

        void resized() override
        {
            auto bounds = getLocalBounds();
            
            // Title area
            bounds.removeFromTop(50);
            
            // Status area
            bounds.removeFromTop(40);
            bounds.removeFromTop(10);
            
            // Frequency displays
            bounds.removeFromTop(60); // Current freq
            bounds.removeFromTop(60); // Target freq
            bounds.removeFromTop(10);
            
            // Note display
            bounds.removeFromTop(60);
            bounds.removeFromTop(10);
            
            // Correction info
            bounds.removeFromTop(40);
            bounds.removeFromTop(10);
            
            // Smoothing control at bottom
            auto smoothingArea = bounds.removeFromTop(40).reduced(20, 5);
            smoothingSlider.setBounds(smoothingArea);
        }

    private:
        void drawFrequencyDisplay(juce::Graphics& g, juce::Rectangle<int> area, float freq, const juce::String& label)
        {
            g.setColour(juce::Colours::lightgrey);
            g.setFont(14.0f);
            g.drawText(label, area.removeFromTop(20),
                       juce::Justification::centred);
            
            if (freq > 0.0f && freq < 5000.0f)
            {
                g.setColour(juce::Colours::cyan);
                g.setFont(32.0f);
                juce::String freqText = juce::String(freq, 1) + " Hz";
                g.drawText(freqText, area, juce::Justification::centred);
            }
            else
            {
                g.setColour(juce::Colours::darkgrey);
                g.setFont(24.0f);
                g.drawText("---", area, juce::Justification::centred);
            }
        }
        
        void drawNoteDisplay(juce::Graphics& g, juce::Rectangle<int> area, float currentFreq, float targetFreq)
        {
            int currentMidiNote = audioProcessor.frequencyToMidiNote(currentFreq);
            int targetMidiNote = audioProcessor.frequencyToMidiNote(targetFreq);
            juce::String currentNoteName = audioProcessor.midiNoteToName(currentMidiNote);
            juce::String targetNoteName = audioProcessor.midiNoteToName(targetMidiNote);
            
            g.setColour(juce::Colours::lightgrey);
            g.setFont(14.0f);
            g.drawText("Note Correction:", area.removeFromTop(20),
                       juce::Justification::centred);
            
            auto noteDisplayArea = area.reduced(20, 0);
            
            if (currentMidiNote >= 0 && targetMidiNote >= 0)
            {
                g.setColour(juce::Colours::yellow);
                g.setFont(20.0f);
                g.drawText(currentNoteName, noteDisplayArea.removeFromLeft(100),
                           juce::Justification::centred);
                
                g.setColour(juce::Colours::white);
                g.setFont(24.0f);
                g.drawText("→", noteDisplayArea.withWidth(40),
                           juce::Justification::centred);
                
                g.setColour(juce::Colours::limegreen);
                g.setFont(24.0f);
                g.drawText(targetNoteName, noteDisplayArea.withLeft(noteDisplayArea.getX() + 40),
                           juce::Justification::centred);
            }
            else
            {
                g.setColour(juce::Colours::darkgrey);
                g.setFont(24.0f);
                g.drawText("--- → ---", area, juce::Justification::centred);
            }
        }
        
        void drawCorrectionInfo(juce::Graphics& g, juce::Rectangle<int> area, float pitchRatio)
        {
            g.setColour(juce::Colours::grey);
            g.setFont(12.0f);
            
            juce::String infoText;
            if (std::abs(pitchRatio - 1.0f) > 0.01f)
            {
                float cents = 1200.0f * std::log2(pitchRatio);
                infoText = "Correcting: " + juce::String(cents, 1) + " cents (ratio: "
                         + juce::String(pitchRatio, 3) + ")";
            }
            else
            {
                infoText = "No correction applied";
            }
            
            g.drawText(infoText, area, juce::Justification::centred);
        }
        
        void drawSmoothingDisplay(juce::Graphics& g, juce::Rectangle<int> area)
        {
            g.setColour(juce::Colours::lightgrey);
            g.setFont(14.0f);
            
            float smoothingMs = audioProcessor.getParameters().getRawParameterValue("smoothing")->load();
            juce::String smoothingText = "Smoothing: " + juce::String(smoothingMs, 0) + " ms";
            
            g.drawText(smoothingText, area, juce::Justification::centred);
        }
        
        void timerCallback() override
        {
            repaint();
        }

        TadtuneAudioProcessor& audioProcessor;

        juce::Slider smoothingSlider;
        juce::Label smoothingLabel;
        juce::AudioProcessorValueTreeState::SliderAttachment smoothingAttachment;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessorEditor)
    };

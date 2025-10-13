#pragma once
#include <JuceHeader.h>
#include "PluginProcessor.h"


class TadtuneAudioProcessorEditor : public juce::AudioProcessorEditor,
                                          private juce::Timer
{
public:
    TadtuneAudioProcessorEditor(TadtuneAudioProcessor& p)
        : AudioProcessorEditor(&p), audioProcessor(p)
    {
        // Set window size
        setSize(400, 350);
        
        // Start timer to update UI at 30 Hz (every ~33ms)
        startTimerHz(30);
    }

    ~TadtuneAudioProcessorEditor() override
    {
        stopTimer();
    }

    //==========================================================================
    // PAINTING
    //==========================================================================
    void paint(juce::Graphics& g) override
    {
        // Background
        g.fillAll(juce::Colour(0xff1e1e1e));
        
        auto bounds = getLocalBounds();
        
        // Title
        g.setColour(juce::Colours::white);
        g.setFont(24.0f);
        g.drawText("Tadtune - Phase Vocoder 1049", bounds.removeFromTop(50),
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
    }

    void resized() override
    {
        // Layout components here if you add sliders/buttons
    }

private:
    //==========================================================================
    // HELPER DRAWING METHODS
    //==========================================================================
    
    void drawFrequencyDisplay(juce::Graphics& g, juce::Rectangle<int> area, float freq, const juce::String& label)
    {
        g.setColour(juce::Colours::lightgrey);
        g.setFont(14.0f);
        g.drawText(label, area.removeFromTop(20),
                   juce::Justification::centred);
        
        // Draw frequency in large font
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
            // Current note
            g.setColour(juce::Colours::yellow);
            g.setFont(20.0f);
            g.drawText(currentNoteName, noteDisplayArea.removeFromLeft(100),
                       juce::Justification::centred);
            
            // Arrow
            g.setColour(juce::Colours::white);
            g.setFont(24.0f);
            g.drawText("→", noteDisplayArea.withWidth(40),
                       juce::Justification::centred);
            
            // Target note
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
            // Calculate cents from pitch ratio
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
    
    //==========================================================================
    // TIMER CALLBACK - Updates UI periodically
    //==========================================================================
    void timerCallback() override
    {
        // Trigger repaint to update display
        repaint();
    }

    // Reference to the processor
    TadtuneAudioProcessor& audioProcessor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessorEditor)
};

//==============================================================================
// Implementation of createEditor in processor
//==============================================================================
inline juce::AudioProcessorEditor* TadtuneAudioProcessor::createEditor()
{
    return new TadtuneAudioProcessorEditor(*this);
}

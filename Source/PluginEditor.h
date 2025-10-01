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
        setSize(400, 300);
        
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
        g.drawText("Pitch Detector", bounds.removeFromTop(50),
                   juce::Justification::centred);
        
        // Get current values from processor
        float freq = audioProcessor.getCurrentFrequency();
        float period = audioProcessor.getCurrentPeriod();
        bool detecting = audioProcessor.isInDetectionMode();
        
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
            g.drawText("TRACKING", statusArea, juce::Justification::centred);
        }
        
        bounds.removeFromTop(20);
        
        // Frequency display
        auto freqArea = bounds.removeFromTop(80);
        drawFrequencyDisplay(g, freqArea, freq);
        
        bounds.removeFromTop(20);
        
        // Note name display
        auto noteArea = bounds.removeFromTop(60);
        drawNoteDisplay(g, noteArea, freq);
        
        bounds.removeFromTop(10);
        
        // Period display (technical info)
        auto periodArea = bounds.removeFromTop(30);
        drawPeriodDisplay(g, periodArea, period);
    }

    void resized() override
    {
        // Layout components here if you add sliders/buttons
    }

private:
    //==========================================================================
    // HELPER DRAWING METHODS
    //==========================================================================
    
    void drawFrequencyDisplay(juce::Graphics& g, juce::Rectangle<int> area, float freq)
    {
        g.setColour(juce::Colours::lightgrey);
        g.setFont(14.0f);
        g.drawText("Frequency:", area.removeFromTop(20),
                   juce::Justification::centred);
        
        // Draw frequency in large font
        if (freq > 0.0f && freq < 5000.0f)
        {
            g.setColour(juce::Colours::cyan);
            g.setFont(48.0f);
            juce::String freqText = juce::String(freq, 1) + " Hz";
            g.drawText(freqText, area, juce::Justification::centred);
        }
        else
        {
            g.setColour(juce::Colours::darkgrey);
            g.setFont(32.0f);
            g.drawText("---", area, juce::Justification::centred);
        }
    }
    
    void drawNoteDisplay(juce::Graphics& g, juce::Rectangle<int> area, float freq)
    {
        int midiNote = audioProcessor.frequencyToMidiNote(freq);
        juce::String noteName = audioProcessor.midiNoteToName(midiNote);
        
        g.setColour(juce::Colours::lightgrey);
        g.setFont(14.0f);
        g.drawText("Note:", area.removeFromTop(20),
                   juce::Justification::centred);
        
        if (midiNote >= 0)
        {
            g.setColour(juce::Colours::yellow);
            g.setFont(36.0f);
            g.drawText(noteName, area, juce::Justification::centred);
        }
        else
        {
            g.setColour(juce::Colours::darkgrey);
            g.setFont(24.0f);
            g.drawText("---", area, juce::Justification::centred);
        }
    }
    
    void drawPeriodDisplay(juce::Graphics& g, juce::Rectangle<int> area, float period)
    {
        g.setColour(juce::Colours::grey);
        g.setFont(12.0f);
        
        if (period > 0.0f)
        {
            juce::String periodText = "Period: " + juce::String(period, 2) + " samples";
            g.drawText(periodText, area, juce::Justification::centred);
        }
        else
        {
            g.drawText("Period: ---", area, juce::Justification::centred);
        }
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

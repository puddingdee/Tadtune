#pragma once
#include <JuceHeader.h>
#include "PluginProcessor.h"

class CustomLookAndFeel : public juce::LookAndFeel_V4
{
public:
    CustomLookAndFeel()
    {
        // Load your custom font from binary data
        customTypeface = juce::Typeface::createSystemTypefaceFor(BinaryData::Catalogue_2_0_ttf,
                                                                  BinaryData::Catalogue_2_0_ttfSize);
        setDefaultSansSerifTypeface(customTypeface);
    }
    
    juce::Typeface::Ptr getTypefaceForFont(const juce::Font& font) override
    {
        return customTypeface;
    }
    
private:
    juce::Typeface::Ptr customTypeface;
};

class TadtuneAudioProcessorEditor : public juce::AudioProcessorEditor,
                                    private juce::Timer
{
public:
    TadtuneAudioProcessorEditor(TadtuneAudioProcessor& p)
        : AudioProcessorEditor(&p),
          audioProcessor(p),
          smoothingAttachment(audioProcessor.getParameters(), "smoothing", smoothingSlider),
          wetGainAttachment(audioProcessor.getParameters(), "wet gain", wetGainSlider),
          dryGainAttachment(audioProcessor.getParameters(), "dry gain", dryGainSlider)
    {
        // Set custom LookAndFeel globally for this editor
        setLookAndFeel(&customLookAndFeel);
        
        // Set window size to match your background image dimensions
        setSize(882, 671);
        
        // Setup smoothing slider as rotary knob (invisible - we'll draw our own)
        addAndMakeVisible(smoothingSlider);
        smoothingSlider.setSliderStyle(juce::Slider::RotaryHorizontalVerticalDrag);
        smoothingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
        smoothingSlider.setColour(juce::Slider::rotarySliderFillColourId, juce::Colours::transparentBlack);
        smoothingSlider.setColour(juce::Slider::rotarySliderOutlineColourId, juce::Colours::transparentBlack);
        smoothingSlider.setColour(juce::Slider::thumbColourId, juce::Colours::transparentBlack);
        
        // Setup wet/dry sliders as linear vertical (invisible - we'll draw our own)
        addAndMakeVisible(wetGainSlider);
        wetGainSlider.setSliderStyle(juce::Slider::LinearVertical);
        wetGainSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
        wetGainSlider.setColour(juce::Slider::trackColourId, juce::Colours::transparentBlack);
        wetGainSlider.setColour(juce::Slider::backgroundColourId, juce::Colours::transparentBlack);
        wetGainSlider.setColour(juce::Slider::thumbColourId, juce::Colours::transparentBlack);
        
        addAndMakeVisible(dryGainSlider);
        dryGainSlider.setSliderStyle(juce::Slider::LinearVertical);
        dryGainSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
        dryGainSlider.setColour(juce::Slider::trackColourId, juce::Colours::transparentBlack);
        dryGainSlider.setColour(juce::Slider::backgroundColourId, juce::Colours::transparentBlack);
        dryGainSlider.setColour(juce::Slider::thumbColourId, juce::Colours::transparentBlack);
        
        // Setup tuning speed label
        addAndMakeVisible(tuningSpeedLabel);
        tuningSpeedLabel.setText("tuning speed", juce::dontSendNotification);
        tuningSpeedLabel.setJustificationType(juce::Justification::centred);
        tuningSpeedLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        tuningSpeedLabel.setFont(juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font())).withHeight(14.0f));
        
        // Setup value display label
        addAndMakeVisible(valueDisplayLabel);
        valueDisplayLabel.setJustificationType(juce::Justification::centred);
        valueDisplayLabel.setColour(juce::Label::textColourId, juce::Colours::cyan);
        valueDisplayLabel.setFont(juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font())).withHeight(12.0f));
        
        // Setup wet/dry labels
        addAndMakeVisible(wetLabel);
        wetLabel.setText("wet", juce::dontSendNotification);
        wetLabel.setJustificationType(juce::Justification::centred);
        wetLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        wetLabel.setFont(juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font())).withHeight(14.0f));
        
        addAndMakeVisible(dryLabel);
        dryLabel.setText("dry", juce::dontSendNotification);
        dryLabel.setJustificationType(juce::Justification::centred);
        dryLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        dryLabel.setFont(juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font())).withHeight(14.0f));
        
        // Load background image
        backgroundImage = juce::ImageCache::getFromMemory(BinaryData::lmao_png,
                                                          BinaryData::lmao_pngSize);
        
        // Load knob image (single PNG that will rotate)
        knobImage = juce::ImageCache::getFromMemory(BinaryData::pinkpad_png,
                                                    BinaryData::pinkpad_pngSize);
        
        // Load slider images (replace with your actual binary data names)
        sliderBackgroundImage = juce::ImageCache::getFromMemory(BinaryData::stems_png,
                                                               BinaryData::stems_pngSize);
        sliderHandleImage = juce::ImageCache::getFromMemory(BinaryData::catfader_png,
                                                           BinaryData::catfader_pngSize);
        
        // Start timer to update info display at 30 Hz
        startTimerHz(30);
        
        // Add listener to update value display when slider changes
        smoothingSlider.onValueChange = [this] { updateValueDisplay(); };
        updateValueDisplay(); // Initial update
    }

    ~TadtuneAudioProcessorEditor() override
    {
        setLookAndFeel(nullptr);
        stopTimer();
    }

    void paint(juce::Graphics& g) override
    {
        // Draw background image
        if (backgroundImage.isValid())
        {
            g.drawImage(backgroundImage, getLocalBounds().toFloat());
        }
        else
        {
            // Fallback if image not loaded
            g.fillAll(juce::Colours::darkgrey);
        }
        
        // Draw custom rotating knob
        if (knobImage.isValid())
        {
            drawRotatingKnob(g);
        }
        
        // Draw wet/dry sliders
        if (sliderBackgroundImage.isValid() && sliderHandleImage.isValid())
        {
            drawSlider(g, wetGainSlider, wetSliderBounds);
            drawSlider(g, dryGainSlider, drySliderBounds);
        }
        
        // Draw info display in top right
        drawInfoDisplay(g);
        drawPluginInfo(g);
    }

    void resized() override
    {
        // Position your rotary knob
        auto knobBounds = juce::Rectangle<int>(450, 400, 200, 200);
        smoothingSlider.setBounds(knobBounds);
        
        // Position label below the knob
        tuningSpeedLabel.setBounds(knobBounds.getX(),
                                  knobBounds.getBottom() - 10,
                                  knobBounds.getWidth(),
                                  20);
        
        // Position value display below the label
        valueDisplayLabel.setBounds(knobBounds.getX(),
                                   tuningSpeedLabel.getBottom(),
                                   knobBounds.getWidth(),
                                   20);
        
        // Position wet/dry sliders in bottom right
        int sliderWidth = 120;
        int sliderHeight = 190;
        int spacing = -30;
        int rightMargin = 40;
        int bottomMargin = 70;
        
        wetSliderBounds = juce::Rectangle<int>(getWidth() - rightMargin - sliderWidth * 2 - spacing,
                                              getHeight() - sliderHeight - bottomMargin,
                                              sliderWidth,
                                              sliderHeight);
                                              
        drySliderBounds = juce::Rectangle<int>(getWidth() - rightMargin - sliderWidth,
                                              getHeight() - sliderHeight - bottomMargin,
                                              sliderWidth,
                                              sliderHeight);
        
        wetGainSlider.setBounds(wetSliderBounds);
        dryGainSlider.setBounds(drySliderBounds);
        
        // Position wet/dry labels
        wetLabel.setBounds(wetSliderBounds.getX(),
                          wetSliderBounds.getBottom() + 5,
                          wetSliderBounds.getWidth(),
                          20);
                          
        dryLabel.setBounds(drySliderBounds.getX(),
                          drySliderBounds.getBottom() + 5,
                          drySliderBounds.getWidth(),
                          20);
    }

private:
    void updateValueDisplay()
    {
        // Convert slider value to milliseconds and display
        float msValue = smoothingSlider.getValue();
        juce::String displayText = juce::String(msValue, 1) + " ms";
        valueDisplayLabel.setText(displayText, juce::dontSendNotification);
    }
    
    void drawRotatingKnob(juce::Graphics& g)
    {
        auto knobBounds = smoothingSlider.getBounds().toFloat();
        
        // Calculate rotation based on slider value
        double value = smoothingSlider.getValue();
        double minValue = smoothingSlider.getMinimum();
        double maxValue = smoothingSlider.getMaximum();
        double normalizedValue = (value - minValue) / (maxValue - minValue);
        
        // Rotate from -135 degrees to +135 degrees (270 degree range)
        float angle = -2.356f + (normalizedValue * 4.712f); // -135° to +135° in radians
        
        // Save graphics state
        g.saveState();
        
        // Move to knob position, scale down the image, and rotate around its center
        auto centerX = knobBounds.getCentreX();
        auto centerY = knobBounds.getCentreY();
        
        auto transform = juce::AffineTransform::translation(-knobImage.getWidth() * 0.5f,
                                                           -knobImage.getHeight() * 0.5f)
                        .scaled(knobBounds.getWidth() / knobImage.getWidth(),
                               knobBounds.getHeight() / knobImage.getHeight())
                        .rotated(angle)
                        .translated(centerX, centerY);
        
        // Draw rotated and scaled image
        g.drawImageTransformed(knobImage, transform);
        
        g.restoreState();
    }
    
    void drawSlider(juce::Graphics& g, juce::Slider& slider, const juce::Rectangle<int>& bounds)
    {
        auto sliderBounds = bounds.toFloat();
        
        // Draw slider background
        if (sliderBackgroundImage.isValid())
        {
            g.drawImage(sliderBackgroundImage, sliderBounds);
        }
        
        // Calculate handle position based on slider value
        double value = slider.getValue();
        double minValue = slider.getMinimum();
        double maxValue = slider.getMaximum();
        double normalizedValue = (value - minValue) / (maxValue - minValue);
        
        // Position handle (vertical slider - bottom is min, top is max)
        float handleHeight = sliderHandleImage.getHeight() * (sliderBounds.getWidth() / sliderHandleImage.getWidth());
        float handleY = sliderBounds.getY() + (sliderBounds.getHeight() - handleHeight) * (1.0f - (float)normalizedValue);
        float handleX = sliderBounds.getCentreX() - (sliderBounds.getWidth() * 0.5f);
        
        // Draw slider handle
        if (sliderHandleImage.isValid())
        {
            g.drawImage(sliderHandleImage,
                       juce::Rectangle<float>(handleX, handleY, sliderBounds.getWidth(), handleHeight));
        }
    }
    
    void drawPluginInfo(juce::Graphics& g)
    {
        g.setFont(juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font())));
        g.setColour(juce::Colours::white);
        g.drawText("tadtune", juce::Rectangle<int>(10, 10, 200, 20), juce::Justification::left);
    }
    
    void drawInfoDisplay(juce::Graphics& g)
    {
        // Use custom font for all text in the info display
        auto customFont = juce::Font(customLookAndFeel.getTypefaceForFont(juce::Font()));
        
        // Get current values from processor
        float currentFreq = audioProcessor.getCurrentFrequency();
        float targetFreq = audioProcessor.getTargetFrequency();
        float pitchRatio = audioProcessor.getPitchRatio();
        
        // Define display area in top right
        auto displayArea = juce::Rectangle<int>(515, 160, 400, 80);
        
        // Note names
        int currentMidiNote = audioProcessor.frequencyToMidiNote(currentFreq);
        int targetMidiNote = audioProcessor.frequencyToMidiNote(targetFreq);
        juce::String currentNoteName = audioProcessor.midiNoteToName(currentMidiNote);
        juce::String targetNoteName = audioProcessor.midiNoteToName(targetMidiNote);
        
        auto noteArea = displayArea.reduced(8, 5);
        
        if (currentMidiNote >= 0 && targetMidiNote >= 0)
        {
            // Current note
            g.setColour(juce::Colours::yellow);
            g.setFont(customFont.withHeight(16.0f));
            g.drawText(currentNoteName, noteArea.removeFromTop(20),
                       juce::Justification::centred);
            
            // Arrow
            g.setColour(juce::Colours::white);
            g.setFont(customFont.withHeight(14.0f));
            g.drawText("|", noteArea.removeFromTop(15),
                       juce::Justification::centred);
            
            // Target note
            g.setColour(juce::Colours::limegreen);
            g.setFont(customFont.withHeight(16.0f));
            g.drawText(targetNoteName, noteArea.removeFromTop(20),
                       juce::Justification::centred);
            
            // Cents correction
            if (std::abs(pitchRatio - 1.0f) > 0.001f)
            {
                float cents = 1200.0f * std::log2(pitchRatio);
                g.setColour(juce::Colours::cyan);
                g.setFont(customFont.withHeight(12.0f));
                juce::String centsText = juce::String(cents > 0 ? "+" : "") +
                                        juce::String(cents, 0) + " ¢";
                g.drawText(centsText, noteArea, juce::Justification::centred);
            }
        }
        else
        {
            // No valid pitch detected
            g.setColour(juce::Colours::grey);
            g.setFont(customFont.withHeight(14.0f));
            g.drawText("no correction", displayArea, juce::Justification::centred);
        }
    }
    
    void timerCallback() override
    {
        repaint();
    }

    TadtuneAudioProcessor& audioProcessor;
    
    CustomLookAndFeel customLookAndFeel;
    juce::Slider smoothingSlider;
    juce::Slider wetGainSlider;
    juce::Slider dryGainSlider;
    juce::AudioProcessorValueTreeState::SliderAttachment smoothingAttachment;
    juce::AudioProcessorValueTreeState::SliderAttachment wetGainAttachment;
    juce::AudioProcessorValueTreeState::SliderAttachment dryGainAttachment;
    juce::Label tuningSpeedLabel;
    juce::Label valueDisplayLabel;
    juce::Label wetLabel;
    juce::Label dryLabel;
    
    juce::Image backgroundImage;
    juce::Image knobImage;
    juce::Image sliderBackgroundImage;
    juce::Image sliderHandleImage;
    
    juce::Rectangle<int> wetSliderBounds;
    juce::Rectangle<int> drySliderBounds;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessorEditor)
};

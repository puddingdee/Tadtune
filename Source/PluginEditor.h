#pragma once
#include <JuceHeader.h>
#include "PluginProcessor.h"

class CustomLookAndFeel : public juce::LookAndFeel_V4
{
public:
    CustomLookAndFeel()
    {
        customTypeface = juce::Typeface::createSystemTypefaceFor(BinaryData::Catalogue_2_0_ttf,
                                                                  BinaryData::Catalogue_2_0_ttfSize);
        setDefaultSansSerifTypeface(customTypeface);
    }
    
    juce::Typeface::Ptr getTypefaceForFont(const juce::Font& font) override
    {
        return customTypeface;
    }
    
    // Override to control button font size
    juce::Font getTextButtonFont(juce::TextButton&, int buttonHeight) override
    {
        // Access buttonFontSize from the editor - for now we'll use a default
        // The font size will be controlled by the buttonFontSize constant
        return juce::Font(customTypeface).withHeight(11.0f);  // CHANGE 11.0f to adjust font size
    }
    
    // Override to draw circular buttons
    void drawButtonBackground(juce::Graphics& g, juce::Button& button, const juce::Colour& backgroundColour,
                            bool shouldDrawButtonAsHighlighted, bool shouldDrawButtonAsDown) override
    {
        auto bounds = button.getLocalBounds().toFloat();
        auto centre = bounds.getCentre();
        auto radius = juce::jmin(bounds.getWidth(), bounds.getHeight()) / 2.0f;
        
        // Apply transparency based on toggle state
        auto colour = backgroundColour;
        if (!button.getToggleState())
        {
            colour = colour.withAlpha(0.3f);  // 30% opacity when not selected
        }
        
        g.setColour(colour);
        g.fillEllipse(bounds.withSizeKeepingCentre(radius * 2, radius * 2));
        
        // Optional: draw outline
        g.setColour(colour.brighter(0.2f));
        g.drawEllipse(bounds.withSizeKeepingCentre(radius * 2, radius * 2), 1.5f);
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
        setSize(882, 671);
        
        // Setup existing controls
        addAndMakeVisible(smoothingSlider);
        smoothingSlider.setSliderStyle(juce::Slider::RotaryHorizontalVerticalDrag);
        smoothingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
        smoothingSlider.setColour(juce::Slider::rotarySliderFillColourId, juce::Colours::transparentBlack);
        smoothingSlider.setColour(juce::Slider::rotarySliderOutlineColourId, juce::Colours::transparentBlack);
        smoothingSlider.setColour(juce::Slider::thumbColourId, juce::Colours::transparentBlack);
        
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
        
        // Setup labels
        addAndMakeVisible(tuningSpeedLabel);
        tuningSpeedLabel.setText("tuning speed", juce::dontSendNotification);
        tuningSpeedLabel.setJustificationType(juce::Justification::centred);
        tuningSpeedLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        tuningSpeedLabel.setFont(juce::Font(14.0f));
        
        addAndMakeVisible(valueDisplayLabel);
        valueDisplayLabel.setJustificationType(juce::Justification::centred);
        valueDisplayLabel.setColour(juce::Label::textColourId, juce::Colours::cyan);
        valueDisplayLabel.setFont(juce::Font(12.0f));
        
        addAndMakeVisible(wetLabel);
        wetLabel.setText("wet", juce::dontSendNotification);
        wetLabel.setJustificationType(juce::Justification::centred);
        wetLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        wetLabel.setFont(juce::Font(14.0f));
        
        addAndMakeVisible(dryLabel);
        dryLabel.setText("dry", juce::dontSendNotification);
        dryLabel.setJustificationType(juce::Justification::centred);
        dryLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        dryLabel.setFont(juce::Font(14.0f));

        // Setup scale section label
        addAndMakeVisible(scaleLabel);
        scaleLabel.setText("scale select", juce::dontSendNotification);
        scaleLabel.setJustificationType(juce::Justification::centred);
        scaleLabel.setColour(juce::Label::textColourId, juce::Colours::white);
        scaleLabel.setFont(juce::Font(14.0f));

        try
        {
            setupScaleButtons();
        }
        catch (const std::exception& e)
        {
            DBG("Error setting up scale buttons: " << e.what());
        }

        // Load images
        backgroundImage = juce::ImageCache::getFromMemory(BinaryData::lmao_png, BinaryData::lmao_pngSize);
        knobImage = juce::ImageCache::getFromMemory(BinaryData::pinkpad_png, BinaryData::pinkpad_pngSize);
        sliderBackgroundImage = juce::ImageCache::getFromMemory(BinaryData::stems_png, BinaryData::stems_pngSize);
        sliderHandleImage = juce::ImageCache::getFromMemory(BinaryData::catfader_png, BinaryData::catfader_pngSize);
        
        startTimerHz(30);
        
        smoothingSlider.onValueChange = [this] { updateValueDisplay(); };
        
        updateValueDisplay();
        
        if (scaleButtons.size() > 0)
        {
            updateScaleButtonStates();
        }
        
        setLookAndFeel(&customLookAndFeel);
    }

    ~TadtuneAudioProcessorEditor() override
    {
        setLookAndFeel(nullptr);
        stopTimer();
    }

    void paint(juce::Graphics& g) override
    {
        if (backgroundImage.isValid())
        {
            g.drawImage(backgroundImage, getLocalBounds().toFloat());
        }
        else
        {
            g.fillAll(juce::Colours::darkgrey);
        }
        
        if (knobImage.isValid())
        {
            drawRotatingKnob(g);
        }
        
        if (sliderBackgroundImage.isValid() && sliderHandleImage.isValid())
        {
            drawSlider(g, wetGainSlider, wetSliderBounds);
            drawSlider(g, dryGainSlider, drySliderBounds);
        }
        
        drawInfoDisplay(g);
        drawPluginInfo(g);
    }

    void resized() override
    {
        auto knobBounds = juce::Rectangle<int>(450, 400, 200, 200);
        smoothingSlider.setBounds(knobBounds);
        tuningSpeedLabel.setBounds(knobBounds.getX(), knobBounds.getBottom() - 10, knobBounds.getWidth(), 20);
        valueDisplayLabel.setBounds(knobBounds.getX(), tuningSpeedLabel.getBottom(), knobBounds.getWidth(), 20);
        
        int sliderWidth = 120;
        int sliderHeight = 190;
        int spacing = -30;
        int rightMargin = 40;
        int bottomMargin = 70;
        
        wetSliderBounds = juce::Rectangle<int>(getWidth() - rightMargin - sliderWidth * 2 - spacing,
                                              getHeight() - sliderHeight - bottomMargin,
                                              sliderWidth, sliderHeight);
        drySliderBounds = juce::Rectangle<int>(getWidth() - rightMargin - sliderWidth,
                                              getHeight() - sliderHeight - bottomMargin,
                                              sliderWidth, sliderHeight);
        
        wetGainSlider.setBounds(wetSliderBounds);
        dryGainSlider.setBounds(drySliderBounds);
        
        wetLabel.setBounds(wetSliderBounds.getX(), wetSliderBounds.getBottom() + 5, wetSliderBounds.getWidth(), 20);
        dryLabel.setBounds(drySliderBounds.getX(), drySliderBounds.getBottom() + 5, drySliderBounds.getWidth(), 20);
        
        scaleLabel.setBounds(circleOfFifthsCenterX - 210, circleOfFifthsCenterY - circleOfFifthsRadius - 40, 150, 20);
        
        if (scaleButtons.size() > 0)
        {
            positionScaleButtons();
        }
    }

private:
    // Circle of fifths configuration
    static constexpr int circleOfFifthsCenterX = 263;
    static constexpr int circleOfFifthsCenterY = 325;
    static constexpr int circleOfFifthsRadius = 145;
    static constexpr int buttonSize = 60;
    
    void setupScaleButtons()
    {
        // Circle of fifths order: C, G, D, A, E, B, Gb/F#, Db, Ab, Eb, Bb, F, (Chr in center)
        const int circleOfFifthsOrder[] = {
            1,  // C Major / A Minor
            8,  // G Major / E Minor
            3,  // D Major / B Minor
            10, // A Major / F# Minor
            5,  // E Major / C# Minor
            12, // B Major / G# Minor
            7,  // Gb Major / Eb Minor
            2,  // Db Major / Bb Minor
            9,  // Ab Major / F Minor
            4,  // Eb Major / C Minor
            11, // Bb Major / G Minor
            6   // F Major / D Minor
        };
        
        const char* circleOfFifthsNames[] = {
            "C/a", "G/e", "D/b", "A/f#", "E/c#", "B/Cb/g#/ab",
            "F#/Gb/d#/eb", "Db/C#/bb/a#", "Ab/f", "Eb/c", "Bb/g", "F/d"
        };
        // Create the 12 circle buttons
        for (int i = 0; i < 12; ++i)
        {
            auto* button = new juce::TextButton(circleOfFifthsNames[i]);
            
            if (button == nullptr)
                continue;
            
            button->setColour(juce::TextButton::buttonColourId, juce::Colour(0x3A4702));
            button->setColour(juce::TextButton::textColourOffId, juce::Colours::white);
            button->setColour(juce::TextButton::buttonOnColourId, juce::Colours::lightcoral);
            button->setColour(juce::TextButton::textColourOnId, juce::Colours::black);
            
            // Apply custom font to button
            button->setLookAndFeel(&customLookAndFeel);
            
            button->setClickingTogglesState(true);
            button->setRadioGroupId(1001);
            
            int scaleIndex = circleOfFifthsOrder[i];
            button->onClick = [this, scaleIndex] {
                try
                {
                    auto* param = audioProcessor.getParameters().getParameter("scale");
                    if (param != nullptr)
                    {
                        float normalizedValue = static_cast<float>(scaleIndex) / 12.0f;
                        param->setValueNotifyingHost(normalizedValue);
                    }
                    updateScaleButtonStates();
                }
                catch (...) {}
            };
            
            // Position in circle
            float angle = (i / 12.0f) * juce::MathConstants<float>::twoPi - juce::MathConstants<float>::halfPi;
            int x = circleOfFifthsCenterX + static_cast<int>(circleOfFifthsRadius * std::cos(angle)) - buttonSize / 2;
            int y = circleOfFifthsCenterY + static_cast<int>(circleOfFifthsRadius * std::sin(angle)) - buttonSize / 2;
            
            // Ensure button is within bounds
            x = juce::jlimit(0, 882 - buttonSize, x);
            y = juce::jlimit(0, 671 - buttonSize, y);
            
            button->setBounds(x, y, buttonSize, buttonSize);
            
            scaleButtons.add(button);
            addAndMakeVisible(button);
        }
        
        // Create chromatic button in the center
        auto* chrButton = new juce::TextButton("chr");
        
        if (chrButton != nullptr)
        {
            chrButton->setColour(juce::TextButton::buttonColourId, juce::Colour(0x3A4702));
            chrButton->setColour(juce::TextButton::textColourOffId, juce::Colours::white);
            chrButton->setColour(juce::TextButton::buttonOnColourId, juce::Colours::lightcoral);
            chrButton->setColour(juce::TextButton::textColourOnId, juce::Colours::black);
            
            // Apply custom font to button
            chrButton->setLookAndFeel(&customLookAndFeel);
            
            chrButton->setClickingTogglesState(true);
            chrButton->setRadioGroupId(1001);
            
            chrButton->onClick = [this] {
                try
                {
                    auto* param = audioProcessor.getParameters().getParameter("scale");
                    if (param != nullptr)
                    {
                        param->setValueNotifyingHost(0.0f);
                    }
                    updateScaleButtonStates();
                }
                catch (...) {}
            };
            
            int x = circleOfFifthsCenterX - buttonSize / 2;
            int y = circleOfFifthsCenterY - buttonSize / 2;
            
            x = juce::jlimit(0, 882 - buttonSize, x);
            y = juce::jlimit(0, 671 - buttonSize, y);
            
            chrButton->setBounds(x, y, buttonSize, buttonSize);
            
            scaleButtons.add(chrButton);
            addAndMakeVisible(chrButton);
        }
    }
    
    void updateScaleButtonStates()
    {
        try
        {
            auto* param = audioProcessor.getParameters().getParameter("scale");
            if (param == nullptr || scaleButtons.size() < 13)
                return;
            
            float normalizedValue = param->getValue();
            int currentScale = juce::roundToInt(normalizedValue * 12.0f);
            
            // Map scale enum back to button index
            const int scaleToButtonIndex[] = {
                12, // Chromatic -> center button
                0,  // C
                7,  // Db
                2,  // D
                9,  // Eb
                4,  // E
                11, // F
                6,  // Gb
                1,  // G
                8,  // Ab
                3,  // A
                10, // Bb
                5   // B
            };
            
            int buttonIndex = scaleToButtonIndex[currentScale];
            
            for (int i = 0; i < scaleButtons.size(); ++i)
            {
                if (scaleButtons[i] != nullptr)
                {
                    scaleButtons[i]->setToggleState(i == buttonIndex, juce::dontSendNotification);
                }
            }
        }
        catch (...) {}
    }
    
    void positionScaleButtons()
    {
        try
        {
            // Reposition circle buttons
            for (int i = 0; i < 12 && i < scaleButtons.size(); ++i)
            {
                if (scaleButtons[i] == nullptr)
                    continue;
                    
                float angle = (i / 12.0f) * juce::MathConstants<float>::twoPi - juce::MathConstants<float>::halfPi;
                int x = circleOfFifthsCenterX + static_cast<int>(circleOfFifthsRadius * std::cos(angle)) - buttonSize / 2;
                int y = circleOfFifthsCenterY + static_cast<int>(circleOfFifthsRadius * std::sin(angle)) - buttonSize / 2;
                
                x = juce::jlimit(0, 882 - buttonSize, x);
                y = juce::jlimit(0, 671 - buttonSize, y);
                
                scaleButtons[i]->setBounds(x, y, buttonSize, buttonSize);
            }
            
            // Reposition center chromatic button
            if (scaleButtons.size() > 12 && scaleButtons[12] != nullptr)
            {
                int x = circleOfFifthsCenterX - buttonSize / 2;
                int y = circleOfFifthsCenterY - buttonSize / 2;
                
                x = juce::jlimit(0, 882 - buttonSize, x);
                y = juce::jlimit(0, 671 - buttonSize, y);
                
                scaleButtons[12]->setBounds(x, y, buttonSize, buttonSize);
            }
        }
        catch (...) {}
    }
    
    void updateValueDisplay()
    {
        float msValue = smoothingSlider.getValue();
        juce::String displayText = juce::String(msValue, 1) + " ms";
        valueDisplayLabel.setText(displayText, juce::dontSendNotification);
    }
    
    void drawRotatingKnob(juce::Graphics& g)
    {
        auto knobBounds = smoothingSlider.getBounds().toFloat();
        double value = smoothingSlider.getValue();
        double minValue = smoothingSlider.getMinimum();
        double maxValue = smoothingSlider.getMaximum();
        double normalizedValue = (value - minValue) / (maxValue - minValue);
        float angle = -2.356f + (normalizedValue * 4.712f);
        
        g.saveState();
        auto centerX = knobBounds.getCentreX();
        auto centerY = knobBounds.getCentreY();
        
        auto transform = juce::AffineTransform::translation(-knobImage.getWidth() * 0.5f, -knobImage.getHeight() * 0.5f)
                        .scaled(knobBounds.getWidth() / knobImage.getWidth(), knobBounds.getHeight() / knobImage.getHeight())
                        .rotated(angle)
                        .translated(centerX, centerY);
        
        g.drawImageTransformed(knobImage, transform);
        g.restoreState();
    }
    
    void drawSlider(juce::Graphics& g, juce::Slider& slider, const juce::Rectangle<int>& bounds)
    {
        auto sliderBounds = bounds.toFloat();
        
        if (sliderBackgroundImage.isValid())
            g.drawImage(sliderBackgroundImage, sliderBounds);
        
        double value = slider.getValue();
        double minValue = slider.getMinimum();
        double maxValue = slider.getMaximum();
        double normalizedValue = (value - minValue) / (maxValue - minValue);
        
        float handleHeight = sliderHandleImage.getHeight() * (sliderBounds.getWidth() / sliderHandleImage.getWidth());
        float handleY = sliderBounds.getY() + (sliderBounds.getHeight() - handleHeight) * (1.0f - (float)normalizedValue);
        float handleX = sliderBounds.getCentreX() - (sliderBounds.getWidth() * 0.5f);
        
        if (sliderHandleImage.isValid())
        {
            g.drawImage(sliderHandleImage, juce::Rectangle<float>(handleX, handleY, sliderBounds.getWidth(), handleHeight));
        }
    }
    
    void drawPluginInfo(juce::Graphics& g)
    {
        g.setFont(juce::Font(16.0f));
        g.setColour(juce::Colours::white);
        g.drawText("tadtune", juce::Rectangle<int>(10, 10, 200, 20), juce::Justification::left);
    }
    
    void drawInfoDisplay(juce::Graphics& g)
    {
        float currentFreq = audioProcessor.getCurrentFrequency();
        float targetFreq = audioProcessor.getTargetFrequency();
        float pitchRatio = audioProcessor.getPitchRatio();
        
        auto displayArea = juce::Rectangle<int>(515, 160, 400, 80);
        int currentMidiNote = audioProcessor.frequencyToMidiNote(currentFreq);
        int targetMidiNote = audioProcessor.frequencyToMidiNote(targetFreq);
        juce::String currentNoteName = audioProcessor.midiNoteToName(currentMidiNote);
        juce::String targetNoteName = audioProcessor.midiNoteToName(targetMidiNote);
        
        auto noteArea = displayArea.reduced(8, 5);
        
        if (currentMidiNote >= 0 && targetMidiNote >= 0)
        {
            g.setColour(juce::Colours::yellow);
            g.setFont(juce::Font(16.0f));
            g.drawText(currentNoteName, noteArea.removeFromTop(20), juce::Justification::centred);
            
            g.setColour(juce::Colours::white);
            g.setFont(juce::Font(14.0f));
            g.drawText("|", noteArea.removeFromTop(15), juce::Justification::centred);
            
            g.setColour(juce::Colours::limegreen);
            g.setFont(juce::Font(16.0f));
            g.drawText(targetNoteName, noteArea.removeFromTop(20), juce::Justification::centred);
            
            if (std::abs(pitchRatio - 1.0f) > 0.001f)
            {
                float cents = 1200.0f * std::log2(pitchRatio);
                g.setColour(juce::Colours::cyan);
                g.setFont(juce::Font(12.0f));
                juce::String centsText = juce::String(cents > 0 ? "+" : "") + juce::String(cents, 0) + " ¢";
                g.drawText(centsText, noteArea, juce::Justification::centred);
            }
        }
        else
        {
            g.setColour(juce::Colours::grey);
            g.setFont(juce::Font(14.0f));
            g.drawText("no correction", displayArea, juce::Justification::centred);
        }
    }

    void timerCallback() override
    {
        if (scaleButtons.size() > 0)
        {
            updateScaleButtonStates();
        }
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

    juce::Label scaleLabel;
    juce::OwnedArray<juce::TextButton> scaleButtons;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(TadtuneAudioProcessorEditor)
};

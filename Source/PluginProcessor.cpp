#include "PluginProcessor.h"
#include "PluginEditor.h"

// Create processor instance
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TadtuneAudioProcessor();
}


juce::AudioProcessorEditor* TadtuneAudioProcessor::createEditor()
{
    return new TadtuneAudioProcessorEditor(*this);
}

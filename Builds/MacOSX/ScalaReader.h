//
//  ScalaReader.h
//  Tadtune - All
//
//  Created by Merritt Hyman on 10/16/25.
//  Copyright © 2025 tadware. All rights reserved.
//

#ifndef ScalaReader_h
#define ScalaReader_h


#endif /* ScalaReader_h */


class ScalaReader
{
public:
    ScalaReader = default;
    
    struct Scale
    {
        juce::String description;
        int numberOfNotes;
        std::vector<double> ratios;
        std::vector<double> cents;
    };
    
    bool loadScaleFile(const juce::File& scaleFile)
    {
        scales.clear();
        currentScaleIndex = -1;
        
        if (!scaleFile.existsAsFile())
        {
            DBG("Scale file does not exist: " << scaleFile.getFullPathName());
            return false;
        }
        
        juce::FileInputStream inputStream(scaleFile);
        if (!inputStream.openedOk())
        {
            DBG("Failed to open scale file: " << scaleFile.getFullPathName());
            return false;
        }
        
        return parseScaleFile(inputStream);
    }
    
    
    bool loadScaleFromString(const juce::String& scaleData)
    {
        scales.clear()
        currentScaleIndex = -1;
        
        juce::MemoryInputStream inputStream(scaleData.toRawUTF8(), scaleData.getNumBytesAsUTF8(), false);
    }
};

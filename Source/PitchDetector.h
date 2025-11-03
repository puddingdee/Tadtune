#pragma once
#include <JuceHeader.h>
#include <vector>
#include <cmath>
#include <algorithm>

/* real time pitch detection using YIN. implementation based on "YIN, a fundamental frequency estimator for speech and music" by Alain de Cheveigne and Hideki Kawahara, implementation assistance with Claude AI
 
 1. high pass filter to remove DC offset and low frequencies
 2. buffer management (circular buffer to store incoming audio
 3. energy check skips processing on silence
 4. YIN algorithm: compute squared difference function, normalize it, find pitch period candidate
 5. post processing: improved period accuracy with parabolic interpolation, correct octave errors
 */




class PitchDetector
{
public:
    PitchDetector() = default;
    
    void prepare(double sampleRate, int bufferSize = 1024)
    {
        fs = sampleRate;
        yinBufferSize = bufferSize;
        
        // allocate buffers. audiobuffer is large for wraparound headroom. yinbuffer is small because it's bounded by 80 - 800Hz detection limits.
        audioBuffer.resize(yinBufferSize * 2, 0.0f);
        yinBuffer.resize(yinBufferSize / 2, 0.0f);
        
        // 80 - 800Hz detection range
        minPeriod = static_cast<int>(fs / 800.0);
        maxPeriod = static_cast<int>(fs / 80.0);
        maxPeriod = std::min(maxPeriod, yinBufferSize / 2);
        
        hopSize = 512;
        
        // init pitch history for median filtering
        pitchHistory.resize(5, 0.0f);  // last 5 pitch estimates
        historyIndex = 0;
        
        reset();
    }
    
    void reset()
    {
        writePos = 0;
        samplesInBuffer = 0;
        frameCounter = 0;
        
        currentFrequency = 0.0f;
        confidence = 0.0f;
        
        std::fill(audioBuffer.begin(), audioBuffer.end(), 0.0f);
        std::fill(yinBuffer.begin(), yinBuffer.end(), 0.0f);
        std::fill(pitchHistory.begin(), pitchHistory.end(), 0.0f);
        
        smoothedFrequency = 0.0f;
        smoothedConfidence = 0.0f;
        lastGoodFrequency = 0.0f;
        framesWithoutPitch = 0;
        historyIndex = 0;
        stableFrameCount = 0;
        
        prevInputSample = 0.0f;
        prevFilteredSample = 0.0f;
    }
    
    bool processSample(float sample)
    {
        // apply highpass
        sample = highPassFilter(sample);
        
        // store in ring buffer
        audioBuffer[writePos] = sample;
        writePos = (writePos + 1) % audioBuffer.size();
        samplesInBuffer = std::min(samplesInBuffer + 1, (int)audioBuffer.size());
        
        frameCounter++;
        
        // process when we reach hopsize
        if (frameCounter >= hopSize)
        {
            frameCounter = 0;
            
            if (samplesInBuffer >= yinBufferSize)
            {
                return detectPitch();
            }
        }
        
        return false;
    }
    
    // getters
    double getFrequency() const { return smoothedFrequency; }
    float getConfidence() const { return smoothedConfidence; }
    bool isInDetectionMode() const { return confidence < 0.15f; }
    bool hasValidPitch() const { return smoothedFrequency > 80.0f && smoothedFrequency < 800.0f && smoothedConfidence > 0.1f; }
    
private:
    bool detectPitch()
    {
        // energy check, dont process noise
        if (!hasEnoughEnergy())
        {
            handleNoPitchDetected();
            return false;
        }
        
        // YIN algorithm
        calculateDifferenceFast();
        cumulativeMeanNormalizedDifference();
        
        int tauEstimate = absoluteThreshold();
        
        if (tauEstimate == -1)
        {
            handleNoPitchDetected();
            return false;
        }
        
        // parabolic interpolation
        float betterTau = parabolicInterpolation(tauEstimate);
        // convert to freq
        float detectedFreq = fs / betterTau;
        float currentConfidence = 1.0f - yinBuffer[tauEstimate];
        
        // validate freq range
        if (detectedFreq < 80.0f || detectedFreq > 800.0f)
        {
            handleNoPitchDetected();
            return false;
        }
        
        // octave error correction
        detectedFreq = correctOctaveErrors(detectedFreq);
        
        // add to pitch history for median filtering
        pitchHistory[historyIndex] = detectedFreq;
        historyIndex = (historyIndex + 1) % pitchHistory.size();
        
        // remove outliers
        float medianFreq = getMedianPitch();
        
        // use median if enough data, else use regs
        float finalFreq = (stableFrameCount > 3) ? medianFreq : detectedFreq;
        
        // additional validation
        if (lastGoodFrequency > 0.0f && currentConfidence < 0.6f)
        {
            float ratio = finalFreq / lastGoodFrequency;
            
            // reject if too far from previous, unless very confident
            if (ratio > 1.15f || ratio < 0.87f)
            {
                // use prev freq instead of complete rejection
                finalFreq = lastGoodFrequency;
                currentConfidence *= 0.7f;  // reduce confidence
            }
        }
        
        // update
        updatePitchEstimate(finalFreq, currentConfidence);
        
        return true;
    }
    
    // compares new detection to the last stable pitch and nudges back to plausible value. this avoids locking onto especially present harmonics
    float correctOctaveErrors(float detectedFreq)
    {
        if (lastGoodFrequency > 0.0f)
        {
            float ratio = detectedFreq / lastGoodFrequency;
            
            // fix 8va up
            if (ratio > 1.85f && ratio < 2.15f)
                return detectedFreq * 0.5f;
            
            // fix 8va dwn
            if (ratio > 0.47f && ratio < 0.54f)
                return detectedFreq * 2.0f;
            
            // fix 5th up
            if (ratio > 1.45f && ratio < 1.55f)
                return detectedFreq * (2.0f / 3.0f);
            
            // fix 5th dwn
            if (ratio > 0.65f && ratio < 0.70f)
                return detectedFreq * (3.0f / 2.0f);
        }
        // return regs if nothing
        return detectedFreq;
    }
    
    // stabilizes output by rejecting outliers and short spikes. computes median of recent estimates in pitchHistory
    float getMedianPitch()
    {
        // sortingggg
        std::vector<float> sorted = pitchHistory;
        std::sort(sorted.begin(), sorted.end());
        
        // remove zeros and sub 80hz
        sorted.erase(std::remove_if(sorted.begin(), sorted.end(),
                                    [](float f) { return f < 80.0f; }),
                     sorted.end());
        
        if (sorted.empty())
            return 0.0f;
        
        // return median. odd count, middle elem is median. even, avg middle elems
        size_t mid = sorted.size() / 2;
        if (sorted.size() % 2 == 0)
            return (sorted[mid - 1] + sorted[mid]) / 2.0f;
        else
            return sorted[mid];
    }
    
    bool hasEnoughEnergy()
    {
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        float sumSquares = 0.0f;
        
        for (int i = 0; i < yinBufferSize; i += 4)
        {
            int idx = (startPos + i) % audioBuffer.size();
            float val = audioBuffer[idx];
            sumSquares += val * val;
        }
        
        float rms = std::sqrt(sumSquares / (yinBufferSize / 4));
        return rms > 0.01f;
    }
    
    
    // SDF for a limited set of time lags.
    void calculateDifferenceFast()
    {
        // picks the last yinBufferSize samples and wraps
        int startPos = (writePos - yinBufferSize + audioBuffer.size()) % audioBuffer.size();
        // only evaluate lags up to max expected pitch period plus some
        int searchRange = std::min((int)yinBuffer.size(), maxPeriod + 10);
        
        // for each lag tau in the search range, accumulate the squared difference btween signal and tau shifted signal
        for (int tau = 0; tau < searchRange; ++tau)
        {
            double sum = 0.0;
            int computeLength = yinBufferSize / 3; // /3 samples to save CPU
            
            for (int i = 0; i < computeLength; ++i)
            {
                int idx1 = (startPos + i) % audioBuffer.size();
                int idx2 = (startPos + i + tau) % audioBuffer.size();
                
                float delta = audioBuffer[idx1] - audioBuffer[idx2];
                sum += delta * delta;
            }
            // result. low vals are good
            yinBuffer[tau] = sum;
        }
        
        for (int tau = searchRange; tau < yinBuffer.size(); ++tau)
        {
            // cringe invalid result do not use
            yinBuffer[tau] = 1e10f;
        }
    }
    
    // converts raw SDF to cumulative mean normalized difference CMND. reduces bias towards small lags and stabilizes pitch detection, allowing for more precise threshold
    void cumulativeMeanNormalizedDifference()
    {
        // lag 0 set to bad value bc its not meaningful
        yinBuffer[0] = 1.0f;
        // define search range
        double runningSum = 0.0;
        int searchRange = std::min((int)yinBuffer.size(), maxPeriod + 10);
        
        // maintain a running sum of raw SDF. runningSum / tau is cumulative mean of SDF up to tau. dividing raw sdf by cumulative mean normalizes it.
        for (int tau = 1; tau < searchRange; ++tau)
        {
            runningSum += yinBuffer[tau];
            yinBuffer[tau] = (runningSum == 0.0) ? 1.0f : yinBuffer[tau] * tau / runningSum;
        }
    }
    
    // scans CMND curve to pick best candidate period tau for the fundamental. finds the first region where CMND drops below threshold, then selects local minimum. returns -1 if not found
    int absoluteThreshold()
    {
        float threshold = 0.15f;
        int tau = minPeriod;
        
        while (tau < maxPeriod && yinBuffer[tau] >= threshold)
            tau++;
        
        if (tau >= maxPeriod)
            return -1;
        
        int minTau = tau;
        float minVal = yinBuffer[tau];
        
        while (tau < maxPeriod && yinBuffer[tau] < threshold)
        {
            if (yinBuffer[tau] < minVal)
            {
                minVal = yinBuffer[tau];
                minTau = tau;
            }
            tau++;
        }
        
        return (minVal > 0.5f) ? -1 : minTau;
    }
    
    // refines the coarse period estimate returned by absoluteThreshold() to a subsample resolution. finds a parabola to the CMND values at tauEstimate - 1, tauEstimate, and +1 then returns location of the parabola minimum
    float parabolicInterpolation(int tauEstimate)
    {
        // need neighbors on both sides
        if (tauEstimate < 1 || tauEstimate >= maxPeriod - 1)
            return static_cast<float>(tauEstimate);
        
        // sample the three points
        float s0 = yinBuffer[tauEstimate - 1];
        float s1 = yinBuffer[tauEstimate];
        float s2 = yinBuffer[tauEstimate + 1];
        
        // find parabola's vertex offset
        float denom = 2.0f * (2.0f * s1 - s2 - s0);
        if (std::abs(denom) < 1e-10f)
            return static_cast<float>(tauEstimate);
        
        float adjustment = (s2 - s0) / denom;
        adjustment = std::max(-0.5f, std::min(0.5f, adjustment));
        
        // return refined lag, a fractional tau that's more accurate than before when converted to freq
        return tauEstimate + adjustment;
    }
    
    // update
    void updatePitchEstimate(float newFreq, float newConfidence)
    {
        currentFrequency = newFreq;
        confidence = newConfidence;
        
        // count stable frames
        if (lastGoodFrequency > 0.0f)
        {
            float ratio = newFreq / lastGoodFrequency;
            if (ratio > 0.98f && ratio < 1.02f)
                stableFrameCount++;
            else
                stableFrameCount = 0;
        }
        
        // confidence based smoothing
        float smoothingFactor;
        
        if (newConfidence > 0.7f && stableFrameCount > 5)
        {
            // high confidence and stable, light smoothing
            smoothingFactor = 0.15f;
        }
        else if (newConfidence > 0.5f)
        {
            // medium confidence, moderate smoothing
            smoothingFactor = 0.35f;
        }
        else
        {
            // low confidence, heavy smoothing
            smoothingFactor = 0.6f;
        }
        
        if (smoothedFrequency == 0.0f)
        {
            smoothedFrequency = newFreq;
            smoothedConfidence = newConfidence;
        }
        else
        {
            smoothedFrequency = smoothingFactor * smoothedFrequency + (1.0f - smoothingFactor) * newFreq;
            smoothedConfidence = 0.6f * smoothedConfidence + 0.4f * newConfidence;
        }
        
        lastGoodFrequency = smoothedFrequency;
        framesWithoutPitch = 0;
    }
    
    // avoids abrupt jumps on no pitch detected, eventually resets the state if silence or uncertainty persists
    void handleNoPitchDetected()
    {
        confidence = 0.0f;
        framesWithoutPitch++;
        stableFrameCount = 0;
        
        // gradually reduce smoothedconfidence on consecutive misses
        if (framesWithoutPitch > 3)
        {
            smoothedConfidence *= 0.7f;
            
            // if too low, clear tracked freq and history to avoid sticking somewhere we shouldnt be
            if (smoothedConfidence < 0.05f)
            {
                smoothedFrequency = 0.0f;
                smoothedConfidence = 0.0f;
                lastGoodFrequency = 0.0f;
                std::fill(pitchHistory.begin(), pitchHistory.end(), 0.0f);
            }
        }
    }
    
    // first order high pass iir
    // y[n] = α * (y[n−1] + x[n] − x[n−1])
    float highPassFilter(float input)
    {
        const float alpha = 0.995f;
        float output = alpha * (prevFilteredSample + input - prevInputSample);
        prevInputSample = input;
        prevFilteredSample = output;
        return output;
    }
    
    // member variables
    double fs = 44100.0;
    int yinBufferSize = 1024;
    int minPeriod = 55;
    int maxPeriod = 551;
    int hopSize = 512;
    
    std::vector<float> audioBuffer;
    std::vector<float> yinBuffer;
    
    int writePos = 0;
    int samplesInBuffer = 0;
    int frameCounter = 0;
    
    float currentFrequency = 0.0f;
    float confidence = 0.0f;
    float smoothedFrequency = 0.0f;
    float smoothedConfidence = 0.0f;
    float lastGoodFrequency = 0.0f;
    int framesWithoutPitch = 0;
    int stableFrameCount = 0;
    
    std::vector<float> pitchHistory;
    int historyIndex = 0;
    
    float prevInputSample = 0.0f;
    float prevFilteredSample = 0.0f;
};

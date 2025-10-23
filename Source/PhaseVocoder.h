    #pragma once
    #include <JuceHeader.h>
    #include <vector>
    #include <cmath>

    class PhaseVocoder
    {
    public:
        PhaseVocoder() = default;
        
        void prepare(double sampleRate, int fftSize)
        {
            fs = sampleRate;
            frameSize = fftSize;
            hopSize = frameSize / 4;
            
            // Initialize FFT
            fftOrder = (int)std::log2(frameSize);
            fft = std::make_unique<juce::dsp::FFT>(fftOrder);
            
            // Allocate FFT buffers
            fftData.resize(frameSize * 2, 0.0f);  // Complex data (real + imaginary)
            
            // Analysis and synthesis buffers
            inputBuffer.resize(frameSize, 0.0f);
            outputBuffer.resize(frameSize, 0.0f);
            
            // Phase tracking
            lastPhase.resize(frameSize, 0.0f);
            sumPhase.resize(frameSize, 0.0f);
            
            // Analysis and synthesis windows (Hann window)
            analysisWindow.resize(frameSize);
            synthesisWindow.resize(frameSize);
            
            for (int i = 0; i < frameSize; ++i)
            {
                float windowValue = 0.5f * (1.0f - std::cos(2.0f * juce::MathConstants<float>::pi * i / (frameSize - 1)));
                analysisWindow[i] = windowValue;
                synthesisWindow[i] = windowValue;
            }
            
            // Normalize synthesis window for perfect reconstruction with overlap-add
            float windowSum = 0.0f;
            for (int i = 0; i < frameSize; i += hopSize)
            {
                if (i < frameSize)
                    windowSum += synthesisWindow[i] * synthesisWindow[i];
            }
            
            for (int i = 0; i < frameSize; ++i)
            {
                synthesisWindow[i] /= (windowSum / hopSize);
            }
            
            // FIFO buffers for input and output
            inputFifo.resize(frameSize * 2, 0.0f);
            outputFifo.resize(frameSize * 2, 0.0f);
            
            inputFifoPos = 0;
            outputFifoPos = 0;
            
            // Pitch shifting with smoothing
            targetPitchRatio = 1.0f;
            currentPitchRatio = 1.0f;
            smoothingTimeMs = 50.0f; // Default 50ms smoothing
            calculateSmoothingCoefficient();
            
            reset();
        }
        
        void reset()
        {
            std::fill(fftData.begin(), fftData.end(), 0.0f);
            std::fill(inputBuffer.begin(), inputBuffer.end(), 0.0f);
            std::fill(outputBuffer.begin(), outputBuffer.end(), 0.0f);
            std::fill(lastPhase.begin(), lastPhase.end(), 0.0f);
            std::fill(sumPhase.begin(), sumPhase.end(), 0.0f);
            std::fill(inputFifo.begin(), inputFifo.end(), 0.0f);
            std::fill(outputFifo.begin(), outputFifo.end(), 0.0f);
            
            inputFifoPos = 0;
            outputFifoPos = 0;
            samplesSinceLastProcess = 0;
            
            targetPitchRatio = 1.0f;
            currentPitchRatio = 1.0f;
        }
        
        void setPitchShiftRatio(float ratio)
        {
            // Clamp to reasonable values
            targetPitchRatio = std::max(0.5f, std::min(2.0f, ratio));
        }
        
        void setSmoothingTime(float milliseconds)
        {
            smoothingTimeMs = std::max(0.0f, std::min(500.0f, milliseconds)); // Clamp to 0-500ms
            calculateSmoothingCoefficient();
        }
        
        float getCurrentPitchRatio() const { return currentPitchRatio; }
        float getTargetPitchRatio() const { return targetPitchRatio; }
        float getSmoothingTime() const { return smoothingTimeMs; }
        
        float processSample(float inputSample)
        {
            // Smooth pitch ratio changes
            if (std::abs(currentPitchRatio - targetPitchRatio) > 0.001f)
            {
                currentPitchRatio = currentPitchRatio + smoothingCoefficient * (targetPitchRatio - currentPitchRatio);
            }
            else
            {
                currentPitchRatio = targetPitchRatio;
            }
            
            // Write input to FIFO
            inputFifo[inputFifoPos] = inputSample;
            inputFifoPos = (inputFifoPos + 1) % inputFifo.size();
            
            samplesSinceLastProcess++;
            
            // Process a frame when we have enough samples
            if (samplesSinceLastProcess >= hopSize)
            {
                processFrame();
                samplesSinceLastProcess = 0;
            }
            
            // Read output from FIFO
            float outputSample = outputFifo[outputFifoPos];
            outputFifo[outputFifoPos] = 0.0f;  // Clear after reading
            outputFifoPos = (outputFifoPos + 1) % outputFifo.size();
            
            return outputSample;
        }
        
    private:
        void calculateSmoothingCoefficient()
        {
            // Calculate smoothing coefficient based on time constant
            // For 50ms smoothing time, we want ~63% of the way to target in 50ms
            if (smoothingTimeMs <= 0.0f || fs == 0.0f)
            {
                smoothingCoefficient = 1.0f; // No smoothing
            }
            else
            {
                float timeConstantSeconds = smoothingTimeMs / 1000.0f;
                smoothingCoefficient = 1.0f - std::exp(-1.0f / (timeConstantSeconds * fs));
                // Apply per-sample, so we need to adjust for the fact that we update every sample
                smoothingCoefficient = std::max(0.0001f, std::min(1.0f, smoothingCoefficient));
            }
        }
        
        void processFrame()
        {
            // Copy input data from FIFO to input buffer
            int readPos = (inputFifoPos - frameSize + inputFifo.size()) % inputFifo.size();
            for (int i = 0; i < frameSize; ++i)
            {
                inputBuffer[i] = inputFifo[readPos];
                readPos = (readPos + 1) % inputFifo.size();
            }
            
            // Apply analysis window
            for (int i = 0; i < frameSize; ++i)
            {
                fftData[i] = inputBuffer[i] * analysisWindow[i];
            }
            
            // Perform FFT
            fft->performRealOnlyForwardTransform(fftData.data(), true);
            
            // Process in frequency domain with phase vocoder algorithm
            float expectedPhaseInc = 2.0f * juce::MathConstants<float>::pi * hopSize / frameSize;
            
            for (int k = 0; k < frameSize / 2; ++k)
            {
                // Get magnitude and phase
                float real = fftData[k * 2];
                float imag = fftData[k * 2 + 1];
                
                float magnitude = std::sqrt(real * real + imag * imag);
                float phase = std::atan2(imag, real);
                
                // Compute phase difference
                float phaseDiff = phase - lastPhase[k];
                lastPhase[k] = phase;
                
                // Subtract expected phase advance
                phaseDiff -= k * expectedPhaseInc;
                
                // Wrap to [-pi, pi]
                while (phaseDiff > juce::MathConstants<float>::pi)
                    phaseDiff -= 2.0f * juce::MathConstants<float>::pi;
                while (phaseDiff < -juce::MathConstants<float>::pi)
                    phaseDiff += 2.0f * juce::MathConstants<float>::pi;
                
                // Compute true frequency
                float trueFreq = k * expectedPhaseInc + phaseDiff;
                
                // Scale frequency for pitch shifting using SMOOTHED ratio
                trueFreq *= currentPitchRatio;
                
                // Accumulate phase
                sumPhase[k] += trueFreq;
                
                // Compute new complex values
                float newPhase = sumPhase[k];
                fftData[k * 2] = magnitude * std::cos(newPhase);
                fftData[k * 2 + 1] = magnitude * std::sin(newPhase);
            }
            
            // Perform inverse FFT
            fft->performRealOnlyInverseTransform(fftData.data());
            
            // Apply synthesis window and overlap-add to output FIFO
            int writePos = outputFifoPos;
            for (int i = 0; i < frameSize; ++i)
            {
                float sample = fftData[i] * synthesisWindow[i] / frameSize;
                outputFifo[writePos] += sample;
                writePos = (writePos + 1) % outputFifo.size();
            }
        }
        
        double fs = 44100.0;
        int frameSize = 2048;
        int hopSize = 512;
        int fftOrder = 11;
        
        std::unique_ptr<juce::dsp::FFT> fft;
        
        std::vector<float> fftData;
        std::vector<float> inputBuffer;
        std::vector<float> outputBuffer;
        
        std::vector<float> lastPhase;
        std::vector<float> sumPhase;
        
        std::vector<float> analysisWindow;
        std::vector<float> synthesisWindow;
        
        std::vector<float> inputFifo;
        std::vector<float> outputFifo;
        
        int inputFifoPos = 0;
        int outputFifoPos = 0;
        int samplesSinceLastProcess = 0;
        
        // Pitch shifting with smoothing
        float targetPitchRatio = 1.0f;
        float currentPitchRatio = 1.0f;
        float smoothingTimeMs = 50.0f;
        float smoothingCoefficient = 0.001f;
    };

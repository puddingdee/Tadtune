    #pragma once
    #include <JuceHeader.h>
    #include <vector>
    #include <cmath>


/*
 frequency domain pitch shifting using phase vocoding
 
 built with Claude AI assistance
 */
    class PhaseVocoder
    {
    public:
        PhaseVocoder() = default;
        
        void prepare(double sampleRate, int fftSize)
        {
            fs = sampleRate;
            frameSize = fftSize;
            hopSize = frameSize / 4;
            
            // initialize juce fft
            fftOrder = (int)std::log2(frameSize);
            fft = std::make_unique<juce::dsp::FFT>(fftOrder);
            
            // allocate fft buffers
            fftData.resize(frameSize * 2, 0.0f);
           
            // analysis and resynthesis buffers
            inputBuffer.resize(frameSize, 0.0f);
            outputBuffer.resize(frameSize, 0.0f);
            
            // phase tracking
            lastPhase.resize(frameSize, 0.0f);
            sumPhase.resize(frameSize, 0.0f);
            
            // hann window init
            analysisWindow.resize(frameSize);
            synthesisWindow.resize(frameSize);
            
            // make the window
            for (int i = 0; i < frameSize; ++i)
            {
                float windowValue = 0.5f * (1.0f - std::cos(2.0f * juce::MathConstants<float>::pi * i / (frameSize - 1)));
                analysisWindow[i] = windowValue;
                synthesisWindow[i] = windowValue;
            }
            
            // scales the sysnthesis window so the overlap add makes unity gain
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
            
            // circular buffers for streaming io
            inputFifo.resize(frameSize * 2, 0.0f);
            outputFifo.resize(frameSize * 2, 0.0f);
            
            inputFifoPos = 0;
            outputFifoPos = 0;
            
            // init smoothed pitch params
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
            // clamp to reasonable values
            targetPitchRatio = std::max(0.5f, std::min(2.0f, ratio));
        }
        
        void setSmoothingTime(float milliseconds)
        {
            smoothingTimeMs = std::max(0.0f, std::min(500.0f, milliseconds)); // 0 to 500ms
            calculateSmoothingCoefficient();
        }
        
        // gets
        float getCurrentPitchRatio() const { return currentPitchRatio; }
        float getTargetPitchRatio() const { return targetPitchRatio; }
        float getSmoothingTime() const { return smoothingTimeMs; }
        
        // takes in one input sample, advances the smoothed pitch ratio toward the target, triggers pvoc process when enough samples have accumulated
        // returns one output sample by overlap adding from output FIFO
        float processSample(float inputSample)
        {
            // smooth pitch ratio
            if (std::abs(currentPitchRatio - targetPitchRatio) > 0.001f)
            {
                // exponential approach per sample toward targetpitchratio
                currentPitchRatio = currentPitchRatio + smoothingCoefficient * (targetPitchRatio - currentPitchRatio);
            }
            else
            {
                currentPitchRatio = targetPitchRatio;
            }
            
            // write input to fifo
            inputFifo[inputFifoPos] = inputSample;
            inputFifoPos = (inputFifoPos + 1) % inputFifo.size();
            
            // track how many input samples
            samplesSinceLastProcess++;
            
            // when there are enough samples, processframe
            if (samplesSinceLastProcess >= hopSize)
            {
                processFrame();
                samplesSinceLastProcess = 0;
            }
            
            // read output from fifo
            float outputSample = outputFifo[outputFifoPos];
            outputFifo[outputFifoPos] = 0.0f;  // Clear after reading
            outputFifoPos = (outputFifoPos + 1) % outputFifo.size();
            
            return outputSample;
        }
        
        
    private:
        
        // converts user facing smoothing time into per sample exponential smoothing coefficient
        void calculateSmoothingCoefficient()
        {
            // handle edge cases
            if (smoothingTimeMs <= 0.0f || fs == 0.0f)
            {
                smoothingCoefficient = 1.0f; // no smoothing
            }
            else
            {
                // convert time constant to per sample smoothing value
                float timeConstantSeconds = smoothingTimeMs / 1000.0f;
                smoothingCoefficient = 1.0f - std::exp(-1.0f / (timeConstantSeconds * fs));
                // clamp to safe range
                smoothingCoefficient = std::max(0.0001f, std::min(1.0f, smoothingCoefficient));
            }
        }
        
        // reads a frame of audio, converts to frequency domain, tracks bin phases, overlap adds resynthesized time domain frame to the output fifo
        void processFrame()
        {
            // gets most recent frameSize samples aligned with hopSize
            int readPos = (inputFifoPos - frameSize + inputFifo.size()) % inputFifo.size();
            for (int i = 0; i < frameSize; ++i)
            {
                inputBuffer[i] = inputFifo[readPos];
                readPos = (readPos + 1) % inputFifo.size();
            }
            
            // apply analysis window before fft
            for (int i = 0; i < frameSize; ++i)
            {
                fftData[i] = inputBuffer[i] * analysisWindow[i];
            }
            
            // do fft
            fft->performRealOnlyForwardTransform(fftData.data(), true);
            
            // phase vocoder algorithm in frequency domain
            // expectedPhaseInc is the linear phase advance per hop for a sinusoid at bin k
            float expectedPhaseInc = 2.0f * juce::MathConstants<float>::pi * hopSize / frameSize;
            
            // in every bin k
            for (int k = 0; k < frameSize / 2; ++k)
            {
                // get magnitude and phase
                float real = fftData[k * 2];
                float imag = fftData[k * 2 + 1];
                
                float magnitude = std::sqrt(real * real + imag * imag);
                float phase = std::atan2(imag, real);
                
                // compute phase difference
                float phaseDiff = phase - lastPhase[k];
                lastPhase[k] = phase;
                
                // subtract expected phase increment. leaves deviation caused by frequency content between bins
                phaseDiff -= k * expectedPhaseInc;
                
                // wrap to -pi, pi to ensure continuity
                while (phaseDiff > juce::MathConstants<float>::pi)
                    phaseDiff -= 2.0f * juce::MathConstants<float>::pi;
                while (phaseDiff < -juce::MathConstants<float>::pi)
                    phaseDiff += 2.0f * juce::MathConstants<float>::pi;
                
                // find true frequency with phase
                float trueFreq = k * expectedPhaseInc + phaseDiff;
                
                // do the pitch shift
                trueFreq *= currentPitchRatio;
                
                // accumulate phase
                sumPhase[k] += trueFreq;
                
                // compute new complex values
                float newPhase = sumPhase[k];
                fftData[k * 2] = magnitude * std::cos(newPhase);
                fftData[k * 2 + 1] = magnitude * std::sin(newPhase);
            }
            
            // go back to time
            fft->performRealOnlyInverseTransform(fftData.data());
            
            // synthesis window, overlap add
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
        
        float targetPitchRatio = 1.0f;
        float currentPitchRatio = 1.0f;
        float smoothingTimeMs = 50.0f;
        float smoothingCoefficient = 0.001f;
    };

#ifndef __PLUGINPROCESSOR_H_88534BAA__
#define __PLUGINPROCESSOR_H_88534BAA__
#define NUMCOMPRESSORS 3
#define FILTERSPERCOMPRESSOR 2
#define PI 3.14159265358979323846
#include "../JuceLibraryCode/JuceHeader.h"
//#include "EQFilter.h"
#include <math.h>
#include "stdio.h"

class CompressorAudioProcessor  : public AudioProcessor
{
public:
    CompressorAudioProcessor();
    ~CompressorAudioProcessor();
    
	int bufferSize;
    void prepareToPlay (double sampleRate, int samplesPerBlock);
    void releaseResources();
	void processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages);
	
	void compressor(AudioSampleBuffer &buffer, int m, int compressorIndex);// compressor functions
    void compress(AudioSampleBuffer& buffer, int index);
    void crossover(float* buffer, float* crossoverFreq, int filterType, int filterIndex);
	template <class T> const T& max ( const T& a, const T& b );

    AudioProcessorEditor* createEditor();

    bool hasEditor() const;

	AudioPlayHead::CurrentPositionInfo lastPosInfo;
 
	int round(float inn);
    const String getName() const;

    int getNumParameters();

    float getParameter (int index);
    void setParameter (int index, float newValue);

    const String getParameterName (int index);
    const String getParameterText (int index);

    const String getInputChannelName (int channelIndex) const;
    const String getOutputChannelName (int channelIndex) const;
    bool isInputChannelStereoPair (int index) const;
    bool isOutputChannelStereoPair (int index) const;
	bool silenceInProducesSilenceOut() const;
	virtual double getTailLengthSeconds() const {return 0;};
    const double EQgainFactor = 1;
	

    bool acceptsMidi() const;
    bool producesMidi() const;

    int getNumPrograms();
    int getCurrentProgram();
    void setCurrentProgram (int index);
    const String getProgramName (int index);
    void changeProgramName (int index, const String& newName);

    //==============================================================================
    void getStateInformation (MemoryBlock& destData);
    void setStateInformation (const void* data, int sizeInBytes);

	float getThreshold();
	float getRatio();
	float getGain();
	float getAttackTime();
	float getReleaseTime();
    float getCrossoverFreq();
    float getKnee();

	void setThreshold(float T, int compressorIndex);
	void setGain(float G, int compressorIndex);
	void setRatio(float R, int compressorIndex);
	void setAttackTime(float A, int compressorIndex);
	void setReleaseTime(float R, int compressorIndex);
    void setCrossoverFreq(float F, int compressorIndex);
    void setKnee(float W, int compressorIndex);
	void resetAll();

	// parameters

	bool compressorONOFF = true;
	int M;
	bool autoTime;

private:
    AudioSampleBuffer inputBuffer;
    AudioSampleBuffer buffer0;
    AudioSampleBuffer outputBuffer;
    AudioSampleBuffer monoBuffer;

//	int bufferSize;
    //these are used to persist UI's size- values are stored along with filter's other parameters, and UI component will update them when it gets resized.
	int lastUIWidth, lastUIHeight;

    HeapBlock <float> x_g, x_l, y_g, y_l, c; //input, output, control

		// parameters
	float ratio[NUMCOMPRESSORS];
    float threshold[NUMCOMPRESSORS];
    float makeUpGain[NUMCOMPRESSORS];
    float tauAttack[NUMCOMPRESSORS];
    float tauRelease[NUMCOMPRESSORS];
    float alphaAttack[NUMCOMPRESSORS];
    float alphaRelease[NUMCOMPRESSORS];
    float kneeWidth[NUMCOMPRESSORS];
    float yL_prev[NUMCOMPRESSORS];
    float crossoverFreq[NUMCOMPRESSORS-1];

	int nhost;
	int samplerate;
    

    IIRFilter** filter;
    IIRCoefficients coef;
    
    
    
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CompressorAudioProcessor);
};


#endif  // __PLUGINPROCESSOR_H_88534BAA__

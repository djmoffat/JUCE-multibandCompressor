#include "PluginProcessor.h"
#include "PluginEditor.h"
//#define NOCOMPRESSORS 3

const float tau = 200;
CompressorAudioProcessor::CompressorAudioProcessor()
	// Initializer List
	:
	inputBuffer(1,1),
    buffer0(1,1),
    outputBuffer(1,1),
    monoBuffer(1,1),
	nhost(0)
{
	lastUIWidth = 850;
    lastUIHeight = 650;
    lastPosInfo.resetToDefault();
}
CompressorAudioProcessor::~CompressorAudioProcessor()
{
    // This is very much the wrong way to do this,
    // but I cannot find any way of making 2D HeapBlock's.
    // So this will have to do for now.
    
}
//==============================================================================
void CompressorAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Compressor Initialisations
    
	M = getNumInputChannels() * (NUMCOMPRESSORS-1) * FILTERSPERCOMPRESSOR;
    int numEqFilters = M;
    filter = (IIRFilter**)malloc(numEqFilters * sizeof(IIRFilter*));
    if(filter == 0)
        numEqFilters = 0;
    else {
        for(int i = 0; i < numEqFilters; i++)
            filter[i] = new IIRFilter;
    }

    
	samplerate = (float)getSampleRate();
	bufferSize = getBlockSize();
    x_g.allocate(bufferSize, true);
    x_l.allocate(bufferSize, true);
    y_g.allocate(bufferSize, true);
    y_l.allocate(bufferSize, true);
    c.allocate(bufferSize, true);

    //yL_prev=0;
    autoTime = false;
	compressorONOFF = true;
	resetAll();
    
}
void CompressorAudioProcessor::releaseResources()
{
    // When playback stops, you can use this to free up any spare memory, etc.
}


void CompressorAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    if(compressorONOFF){
        int channels = getNumInputChannels();
        // Initialise Buffers
        // Copy input buffer and create a stereo buffer0 to process
        inputBuffer.setSize(channels,bufferSize);
        inputBuffer = buffer;
        buffer0.setSize(2,bufferSize);
        buffer.clear();
        for (int channel = 0 ; channel < channels ; channel+= 2)
        {
            // clear stereo buffer
            buffer0.clear();
            //        buffer.clear(channel, 0,bufferSize);
            //        buffer.clear(channel+1, 0,bufferSize);
            for (int i = 0; i < NUMCOMPRESSORS; i++){
                // For every compressor, make copy of current stereo channels
                buffer0.copyFrom(0,0,inputBuffer,channel,0,bufferSize);
                buffer0.copyFrom(1,0,inputBuffer,channel+1,0,bufferSize);
                float* LchannelData = buffer0.getSampleData(channel);
                float* RchannelData = buffer0.getSampleData(channel+1);
                // Pass buffer to crossover
                crossover(LchannelData, crossoverFreq, i, ((i*2*(NUMCOMPRESSORS-1)-2)+channel));
                crossover(RchannelData, crossoverFreq, i, ((i*2*(NUMCOMPRESSORS-1)-2)+channel)+1);
                //phase invert if compressor number is Even, as required by LR crossover.
                buffer0.applyGain((float)((i%2)*2)-1);
                // Apply compression to buffer, after crossover
                compress(buffer0, i);
                // add to input buffer
                buffer.addFrom(channel, 0, buffer0, channel, 0, bufferSize);
                buffer.addFrom(channel+1, 0, buffer0, channel+1, 0, bufferSize);
                //            k+=2;
            }
            // apply gain to normalise signals back to original volume
            //buffer.applyGain();
        }
    }
}
void CompressorAudioProcessor::crossover(float* channelData, float crossoverFreq[], int filterType, int filterIndex)
{
    jassert(filterType >= 0 && filterType <= NUMCOMPRESSORS);
    if(filterType == 0){ //Filter Type is Low Pass
        filter[filterIndex+2]->setCoefficients(coef.makeLRLowPass(samplerate, crossoverFreq[filterType]));
        filter[filterIndex+2]->processSamples(channelData, bufferSize);
    }
    else if (filterType == (NUMCOMPRESSORS-1)){ //Filter Type is High Pass
        filter[filterIndex]->setCoefficients(coef.makeLRHighPass(samplerate, crossoverFreq[filterType-1]));
        filter[filterIndex]->processSamples(channelData, bufferSize);
    }
    else{ //Filter Type is Band Pass
        //Apply high Pass
        filter[filterIndex]->setCoefficients(coef.makeLRHighPass(samplerate, crossoverFreq[filterType-1]));
        filter[filterIndex]->processSamples(channelData, bufferSize);
        filter[filterIndex+2]->setCoefficients(coef.makeLRLowPass(samplerate, crossoverFreq[filterType]));
        filter[filterIndex+2]->processSamples(channelData, bufferSize);
    }
}

// compress master class
void CompressorAudioProcessor::compress(AudioSampleBuffer& buffer, int compressorIndex)
{
    if ((threshold[compressorIndex]< 0 || makeUpGain[compressorIndex] != 0))
    {
        // Create mono buffer
        monoBuffer.setSize(1,bufferSize);
        monoBuffer.clear();
        // Fill monoBuffer from Buffer
        monoBuffer.addFrom(0, 0,buffer, 0,0,bufferSize, 0.5);
        monoBuffer.addFrom(0, 0,buffer, 1,0,bufferSize, 0.5);
        // Calculate compression Ratio
        compressor(monoBuffer, 0, compressorIndex);

        // apply control voltage to the audio signal
        for (int i = 0 ; i < bufferSize ; ++i)
        {
            buffer.getSampleData(0)[i] *= c[i];
            buffer.getSampleData(1)[i] *= c[i];
        }
    }
}
// compressor functions
void CompressorAudioProcessor::compressor(AudioSampleBuffer &buffer, int channel, int compressorIndex)
{
    if(tauRelease[compressorIndex] < tauAttack[compressorIndex])
        tauRelease[compressorIndex] = tauAttack[compressorIndex];
    int m = channel;
	alphaAttack[compressorIndex] = exp(-1/(0.001 * samplerate * tauAttack[compressorIndex]));
	alphaRelease[compressorIndex]= exp(-1/(0.001 * samplerate * tauRelease[compressorIndex]));
	for (int i = 0 ; i < bufferSize ; ++i)
	{
		//Level detection- estimate level using peak detector
		if (fabs(buffer.getSampleData(m)[i]) < 0.000001)
            x_g[i] =-120;
		else
            x_g[i] =20*log10(fabs(buffer.getSampleData(m)[i]));

        // Apply second order interpolation soft knee
        if ((2* fabs(x_g[i]-threshold[compressorIndex])) <= kneeWidth[compressorIndex])
            y_g[i] = x_g[i] + (1/ratio[compressorIndex] -1) * pow(x_g[i]-threshold[compressorIndex]+kneeWidth[compressorIndex]/2,2) / (2*kneeWidth[compressorIndex]);
        else if ((2*(x_g[i]-threshold[compressorIndex])) > kneeWidth[compressorIndex])
            y_g[i] = threshold[compressorIndex]+ (x_g[i] - threshold[compressorIndex]) / ratio[compressorIndex];            
        else
            y_g[i] = x_g[i];
        
        x_l[i] = x_g[i] - y_g[i];
        
        //Ballistics- smoothing of the gain
		if (x_l[0]>yL_prev[compressorIndex])
            y_l[i]=alphaAttack[compressorIndex] * yL_prev[compressorIndex]+(1 - alphaAttack[compressorIndex] ) * x_l[i] ;
		else
            y_l[i]=alphaRelease[compressorIndex]* yL_prev[compressorIndex]+(1 - alphaRelease[compressorIndex]) * x_l[i] ;
        c[i] = pow(10,(makeUpGain[compressorIndex] - y_l[i])/20);
		yL_prev[compressorIndex]=y_l[i];
	}
}
template <class T> const T& CompressorAudioProcessor::max( const T& a, const T& b )
{
  return (a < b) ? b : a;
}
void CompressorAudioProcessor::resetAll()
{
    
    crossoverFreq[0] = 1000;
    crossoverFreq[1] = 4000;

    for(int k = 0; k < NUMCOMPRESSORS; k++){
		tauAttack[k]=0.1;
        tauRelease[k] = 0.1;
		alphaAttack[k]=0;
        alphaRelease[k] = 0;
		threshold[k] = 0;
		ratio[k]= 1;
		makeUpGain[k]= 0;
		yL_prev[k]=0;
        kneeWidth[k] = 0.1;
    }
    for (int i = 0 ; i < bufferSize ; i++)
    {
        x_g[i] = 0;
        y_g[i] = 0;
        x_l[i] = 0;
        y_l[i] = 0;
          c[i] = 0;
    }

}

//////////////////////////////////////////////
float CompressorAudioProcessor::getThreshold()
{
	return *threshold;
}
float CompressorAudioProcessor::getRatio()
{
	return *ratio;
}
float CompressorAudioProcessor::getGain()
{
	return *makeUpGain;//problem?
}
float CompressorAudioProcessor::getAttackTime()
{
	return *tauAttack;
}
float CompressorAudioProcessor::getReleaseTime()
{
	return *tauRelease;
}
float CompressorAudioProcessor::getCrossoverFreq()
{
    return *crossoverFreq;
}
float CompressorAudioProcessor::getKnee()
{
    return *kneeWidth;
}
////////////////////////////////////////////////////////
void CompressorAudioProcessor::setThreshold(float T, int compressorIndex)
{
	threshold[compressorIndex]= T;
//    assert(compressorIndex != 0);

}
void CompressorAudioProcessor::setGain(float G, int compressorIndex)
{
    makeUpGain[compressorIndex]= G;
}
void CompressorAudioProcessor::setRatio(float R, int compressorIndex)
{
	ratio[compressorIndex]= R;
}
void CompressorAudioProcessor::setAttackTime(float A, int compressorIndex)
{
	tauAttack[compressorIndex] = A;
}
void CompressorAudioProcessor::setReleaseTime(float R, int compressorIndex)
{
	tauRelease[compressorIndex] = R;
}
void CompressorAudioProcessor::setCrossoverFreq(float F, int compressorIndex)
{
    crossoverFreq[compressorIndex] = F;
}
void CompressorAudioProcessor::setKnee(float W, int compressorIndex)
{
    kneeWidth[compressorIndex] = W;
}
bool CompressorAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}
AudioProcessorEditor* CompressorAudioProcessor::createEditor()
{
    return new CompressorAudioProcessorEditor (this);
}
//==============================================================================
void CompressorAudioProcessor::getStateInformation (MemoryBlock& destData)
{
//Use this to store your parameters in memory block, either as raw data, or use XML or ValueTree classes as intermediaries to make it easy to save and load complex data.
}
void CompressorAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
// Use this to restore your parameters from this memory block, whose contents will have been created by the getStateInformation() call.
}
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new CompressorAudioProcessor();
}
int CompressorAudioProcessor::round(float inn)
{
	if (inn > 0) return (int) (inn + 0.5);
	else return (int) (inn - 0.5);
}
const String CompressorAudioProcessor::getName() const
{
    return JucePlugin_Name;
}
int CompressorAudioProcessor::getNumParameters()
{
    return 0;
}
float CompressorAudioProcessor::getParameter (int index)
{
    return 0.0f;
}
void CompressorAudioProcessor::setParameter (int index, float newValue)
{
}
const String CompressorAudioProcessor::getParameterName (int index)
{
    return String::empty;
}
const String CompressorAudioProcessor::getParameterText (int index)
{
    return String::empty;
}
const String CompressorAudioProcessor::getInputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}
const String CompressorAudioProcessor::getOutputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}
bool CompressorAudioProcessor::isInputChannelStereoPair (int index) const
{
    return true;
}
bool CompressorAudioProcessor::isOutputChannelStereoPair (int index) const
{
    return true;
}
bool CompressorAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}
bool CompressorAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}
bool CompressorAudioProcessor::silenceInProducesSilenceOut() const
{
    return false;
}

int CompressorAudioProcessor::getNumPrograms()
{
    return 0;
}
int CompressorAudioProcessor::getCurrentProgram()
{
    return 0;
}
void CompressorAudioProcessor::setCurrentProgram (int index)
{
}
const String CompressorAudioProcessor::getProgramName (int index)
{
    return String::empty;
}
void CompressorAudioProcessor::changeProgramName (int index, const String& newName)
{
}



 
//  IIRCoefficient equation coppied into juce_IIRFilter

/* IIRCoefficients IIRCoefficients::makeLRLowPass (const double sampleRate,
 const double frequency) noexcept
 {
 jassert (sampleRate > 0 && frequency > 0);
 
 const double omega = M_PI * frequency;
 const double omega2 = pow(omega,2);
 const double theta = tan(omega / sampleRate);
 const double k = omega / theta;
 const double k2 = pow(k,2);
 const double delta = k2 + omega2 + 2 * k * omega;
 const double od = omega2/delta;
 const double kd = k / delta;
 
 return IIRCoefficients (od,
 od * 2.0,
 od,
 1.0,
 (-2 * kd) + (2 * od),
 (-2 * k * omega)/delta + kd + od);
 }
 
 IIRCoefficients IIRCoefficients::makeLRHighPass (const double sampleRate,
 const double frequency) noexcept
 {
 jassert (sampleRate > 0 && frequency > 0);
 
 const double omega = M_PI * frequency;
 const double omega2 = pow(omega,2);
 const double theta = tan(omega / sampleRate);
 const double k = omega / theta;
 const double k2 = pow(k,2);
 const double delta = k2 + omega2 + 2 * k * omega;
 const double od = omega2/delta;
 const double kd = k / delta;
 
 return IIRCoefficients (kd,
 kd * -2.0,
 kd,
 1.0,
 (-2 * kd) + (2 * od),
 (-2 * k * omega)/delta + kd + od);
 } 
 */
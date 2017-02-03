/*
  ==============================================================================

  This is an automatically generated GUI class created by the Introjucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Introjucer version: 3.1.0

  ------------------------------------------------------------------------------

  The Introjucer is part of the JUCE library - "Jules' Utility Class Extensions"
  Copyright 2004-13 by Raw Material Software Ltd.

  ==============================================================================
*/

#ifndef __JUCE_HEADER_388D2EBDBB646058__
#define __JUCE_HEADER_388D2EBDBB646058__

//[Headers]     -- You can add your own extra header files here --
//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Jucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class CompressorAudioProcessorEditor  : public AudioProcessorEditor,
                                        public Timer,
                                        public ButtonListener,
                                        public SliderListener
{
public:
    //==============================================================================
    CompressorAudioProcessorEditor (CompressorAudioProcessor* ownerFilter);
    ~CompressorAudioProcessorEditor();

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.
	void timerCallback();
    //[/UserMethods]

    void paint (Graphics& g);
    void resized();
    void buttonClicked (Button* buttonThatWasClicked);
    void sliderValueChanged (Slider* sliderThatWasMoved);

    // Binary resources:
    static const char* brushedMetalDark_jpg;
    static const int brushedMetalDark_jpgSize;
    static const char* c4dm_png2;
    static const int c4dm_png2Size;
    static const char* qmul_png2;
    static const int qmul_png2Size;
    static const char* knobstrip_png;
    static const int knobstrip_pngSize;
    static const char* scaleLr_png;
    static const int scaleLr_pngSize;


private:
    //[UserVariables]   -- You can add your own custom variables in this section.

    ScopedPointer<ResizableCornerComponent> resizer;
    ComponentBoundsConstrainer resizeLimits;



	AudioPlayHead::CurrentPositionInfo lastDisplayedPosition;

    CompressorAudioProcessor* getProcessor(int index) const
    {
        return static_cast <CompressorAudioProcessor*> (getAudioProcessor());
    }

    void displayPositionInfo (const AudioPlayHead::CurrentPositionInfo& pos);

    //[/UserVariables]

    //==============================================================================
    ScopedPointer<Label> label;
    ScopedPointer<TextButton> buttonONOFF;
    ScopedPointer<Slider> sliderThreshold;
    ScopedPointer<Label> label2;
    ScopedPointer<Slider> sliderRatio;
    ScopedPointer<Label> label3;
    ScopedPointer<Slider> sliderGain;
    ScopedPointer<Label> label7;
    ScopedPointer<Slider> sliderAttack;
    ScopedPointer<Label> label5;
    ScopedPointer<Slider> sliderRelease;
    ScopedPointer<Label> label6;
    ScopedPointer<Slider> sliderThreshold2;
    ScopedPointer<Slider> sliderRatio2;
    ScopedPointer<Slider> sliderGain2;
    ScopedPointer<Slider> sliderAttack2;
    ScopedPointer<Slider> sliderRelease2;
    ScopedPointer<Slider> sliderCrossFreq1;
    ScopedPointer<Label> label11;
    ScopedPointer<Slider> sliderThreshold3;
    ScopedPointer<Slider> sliderRatio3;
    ScopedPointer<Slider> sliderGain3;
    ScopedPointer<Slider> sliderAttack3;
    ScopedPointer<Slider> sliderRelease3;
    ScopedPointer<Slider> sliderCrossFreq2;
    ScopedPointer<Label> label18;
    ScopedPointer<Label> label21;
    ScopedPointer<Label> label22;
    ScopedPointer<Slider> sliderKnee;
    ScopedPointer<Label> label4;
    ScopedPointer<Slider> sliderKnee2;
    ScopedPointer<Slider> sliderKnee3;
    ScopedPointer<Label> label8;
    ScopedPointer<Label> label9;
    ScopedPointer<Label> label10;
    ScopedPointer<Label> label12;
    ScopedPointer<Label> label13;
    ScopedPointer<Label> label14;
    ScopedPointer<Label> label15;
    ScopedPointer<Label> label16;
    ScopedPointer<Label> label17;
    ScopedPointer<Label> label19;
    ScopedPointer<Label> label20;
    ScopedPointer<Label> label23;
    Image cachedImage_brushedMetalDark_jpg;


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (CompressorAudioProcessorEditor)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

#endif   // __JUCE_HEADER_388D2EBDBB646058__

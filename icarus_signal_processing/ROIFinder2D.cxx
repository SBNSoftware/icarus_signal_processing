#ifndef __icarus_signal_processing_ROIFinder2D_CXX__
#define __icarus_signal_processing_ROIFinder2D_CXX__

#include "ROIFinder2D.h"

#include "icarus_signal_processing/Filters/AdaptiveWiener.h"
#include "icarus_signal_processing/Detection/Thresholding.h"
#include "icarus_signal_processing/Filters/MiscUtils.h"
#include "icarus_signal_processing/Detection/LineDetection.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Detection/MorphologicalFunctions2D.h"
//#include "FrequencyFilters1D.h"
//#include "MorphologicalCNC.h"

icarus_signal_processing::ROIChainFilter::ROIChainFilter(size_t             FREQUENCY_THRESHOLD,
                                                         size_t             FREQUENCY_FILTER_SMOOTHNESS_ORDER,
                                                         size_t             FREQUENCY_FILTER_MODE,
                                                         char               MORPHOLOGICAL_FILTER_NAME,
                                                         const unsigned int CHANNEL_GROUPING,
                                                         const unsigned int STRUCTURING_ELEMENT_X,
                                                         const unsigned int STRUCTURING_ELEMENT_Y,
                                                         const unsigned int ROI_EXPAND_WINDOW_SIZE,
                                                         const float        MORPHOLOGICAL_THRESHOLD_FACTOR,
                                                         const size_t       THETASTEPS,
                                                         const unsigned int HOUGH_THRESHOLD,
                                                         const unsigned int NMS_WINDOW_SIZE,
                                                         const unsigned int ANGLE_WINDOW,
                                                         // float NOISE_VARIANCE = 20.0;
                                                         const unsigned int ADFILTER_SX,
                                                         const unsigned int ADFILTER_SY,
                                                         const unsigned int BINARY_DILATION_SX,
                                                         const unsigned int BINARY_DILATION_SY,
                                                         const float        GLOBAL_THRESHOLDING_FACTOR) //:
                            /*fFrequencyThreshold(FREQUENCY_THRESHOLD),
                            fFrequency_filter_smoothness_order(FREQUENCY_FILTER_SMOOTHNESS_ORDER),
                            fFrequency_filter_mode(FREQUENCY_FILTER_MODE),
                            fMorphological_filter_name(MORPHOLOGICAL_FILTER_NAME),
                            fChannel_Grouping(CHANNEL_GROUPING),
                            fStructuring_Element_X(STRUCTURING_ELEMENT_X),
                            fStructuring_Element_Y(STRUCTURING_ELEMENT_Y),
                            fROI_Expand_Window_Size(ROI_EXPAND_WINDOW_SIZE),
                            fMorphological_Threshold_Factor(MORPHOLOGICAL_THRESHOLD_FACTOR),
                            fThetaSteps(THETASTEPS),
                            fHough_Threshold(HOUGH_THRESHOLD),
                            fNMS_Window_Size(NMS_WINDOW_SIZE),
                            fAngle_Window(ANGLE_WINDOW),
                            fADFilter_SX(ADFILTER_SX),
                            fADFilter_SY(ADFILTER_SY),
                            fBinary_Dilation_SX(BINARY_DILATION_SX),
                            fBinary_Dilation_SY(BINARY_DILATION_SY),
                            fGlobal_Thresholding_Factor(GLOBAL_THRESHOLDING_FACTOR)*/

{
    return;
}


void icarus_signal_processing::ROIChainFilter::operator()(const IROIFinder2D::Array2D<float>& waveform2D,
                                                          IROIFinder2D::Array2D<float>&       fullEvent,
                                                          IROIFinder2D::Array2D<bool>&        outputROI) const
{
    // All input arrays must have the same dimensions as waveform2D
/*
    int numChannels = waveform2D.size();
    int numTicks = waveform2D.at(0).size();

    fullEvent.resize(numChannels);
    for (auto &v : fullEvent)
    {
        v.resize(numTicks);
    }

    icarus_signal_processing::MiscUtils utils;

    // 1. Remove Pedestals
    for (int i = 0; i < numChannels; ++i)
    {
        float median = utils.computeMedian(waveform2D[i]);
        for (int j = 0; j < numTicks; ++j)
        {
            fullEvent[i][j] = waveform2D[i][j] - median;
        }
    }

    // 2. Buffer for intermediate computations
    Array2D<float> buffer(numChannels);
    Array2D<bool> selectVals(numChannels);
    Array2D<bool> rois(numChannels);
    Array2D<bool> refinedSelectVals(numChannels);

    for (auto &v : buffer)
    {
        v.resize(numTicks);
    }

    for (auto &v : selectVals)
    {
        v.resize(numTicks);
    }

    for (auto &v : rois)
    {
        v.resize(numTicks);
    }

    for (auto &v : refinedSelectVals)
    {
        v.resize(numTicks);
    }

    // 3. Apply frequency high pass filters
    icarus_signal_processing::FrequencyFilters1D freqFilt;
    freqFilt.filterImage(
        fullEvent,
        FREQUENCY_THRESHOLD,
        buffer,
        FREQUENCY_FILTER_SMOOTHNESS_ORDER,
        FREQUENCY_FILTER_MODE);

    // 4. Run Coherent Noise Correction
    icarus_signal_processing::MorphologicalCNC denoiser;

    buffer.resize(numChannels);
    for (auto &v : buffer)
    {
        v.resize(numTicks);
    }
    std::cout << "1" << std::endl;
    denoiser.denoiseHough2D(
        waveLessCoherent,
        morphedWaveform2D,
        buffer,
        selectVals,
        refinedSelectVals,
        rois,
        MORPHOLOGICAL_FILTER_NAME,
        CHANNEL_GROUPING,
        STRUCTURING_ELEMENT_X,
        STRUCTURING_ELEMENT_Y,
        ROI_EXPAND_WINDOW_SIZE,
        MORPHOLOGICAL_THRESHOLD_FACTOR,
        THETASTEPS,
        HOUGH_THRESHOLD,
        NMS_WINDOW_SIZE,
        ANGLE_WINDOW);

    // 5. Run Adaptive (Incoherent) noise filtering
    std::cout << "2" << std::endl;
    icarus_signal_processing::AdaptiveWiener adFilter;
    std::cout << "3" << std::endl;
    buffer.resize(numChannels);
    for (auto &v : buffer)
    {
        v.resize(numTicks);
    }
    std::cout << "4" << std::endl;
    adFilter.MMWFStar(
        buffer,
        waveLessCoherent,
        ADFILTER_SX,
        ADFILTER_SY);
    std::cout << "5" << std::endl;
    finalErosion2D.resize(numChannels);
    for (auto &v : finalErosion2D)
    {
        v.resize(numTicks);
    }
    std::cout << "6" << std::endl;
    icarus_signal_processing::Morph2DFast morph2D;

    buffer.resize(numChannels);
    for (auto &v : buffer)
    {
        v.resize(numTicks);
    }
    std::cout << "7" << std::endl;
    if (MORPHOLOGICAL_FILTER_NAME == 'e')
    {
        morph2D.getErosion(
            buffer,
            STRUCTURING_ELEMENT_X,
            STRUCTURING_ELEMENT_Y,
            finalErosion2D);
    }
    else if (MORPHOLOGICAL_FILTER_NAME == 'd')
    {
        morph2D.getDilation(
            buffer,
            STRUCTURING_ELEMENT_X,
            STRUCTURING_ELEMENT_Y,
            finalErosion2D);
    }
    else if (MORPHOLOGICAL_FILTER_NAME == 'g')
    {
        morph2D.getGradient(
            buffer,
            STRUCTURING_ELEMENT_X,
            STRUCTURING_ELEMENT_Y,
            finalErosion2D);
    }
    else
    {
        morph2D.getGradient(
            buffer,
            STRUCTURING_ELEMENT_X,
            STRUCTURING_ELEMENT_Y,
            finalErosion2D);
    }
    std::cout << "8" << std::endl;
    // Take absolute value
    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            if (finalErosion2D[i][j] < 0)
            {
                finalErosion2D[i][j] = -finalErosion2D[i][j];
            }
        }
    }
    std::cout << "9" << std::endl;
    // 6. Run Thresholding
    icarus_signal_processing::Thresholding thresholder;
    thresholder.globalMean(
        finalErosion2D,
        rois,
        GLOBAL_THRESHOLDING_FACTOR);
    std::cout << "10" << std::endl;
    morph2D.getClosing(rois, BINARY_CLOSING_SX, BINARY_CLOSING_SY, outputROI);
    // 7. Refine ROI
    return;
*/
}

icarus_signal_processing::ROICannyFilter::ROICannyFilter(const IFFTFilterFunction*  imageFilterFunction,
                                                         const IDenoiser2D*         denoising,
                                                         const BilateralFilters*    bilateralFilter,
                                                         const EdgeDetection*       edgeDetector,
                                                         const unsigned int         ADFILTER_SX,
                                                         const unsigned int         ADFILTER_SY,
                                                         const float                sigma_x, 
                                                         const float                sigma_y, 
                                                         const float                sigma_r, 
                                                         const float                lowThreshold,
                                                         const float                highThreshold,
                                                         const unsigned int         BINARY_DILATION_SX,
                                                         const unsigned int         BINARY_DILATION_SY) :
                            fEdgeDetector(edgeDetector),
                            fADFilter_SX(ADFILTER_SX),
                            fADFilter_SY(ADFILTER_SY),
                            fSigma_x(sigma_x),
                            fSigma_y(sigma_y),
                            fSigma_r(sigma_r),
                            fLowThreshold(lowThreshold),
                            fHighThreshold(highThreshold),
                            fBinary_Dilation_SX(BINARY_DILATION_SX),
                            fBinary_Dilation_SY(BINARY_DILATION_SY)

{
    // Now create the image filter
    fImageFilter = std::make_unique<ImageFilters>(imageFilterFunction);


    return;
}


void icarus_signal_processing::ROICannyFilter::operator()(const IROIFinder2D::Array2D<float>& waveform2D,
                                                          IROIFinder2D::Array2D<float>&       fullEvent,
                                                          IROIFinder2D::Array2D<bool>&        outputROI) const
{
    // Note that this function assumes that the input image has been pedestal corrected and that any coherent noise 
    // has already been removed. 
    //
    // All input arrays must have the same dimensions as waveform2D
    size_t numChannels = waveform2D.size();
    size_t numTicks    = waveform2D[0].size();

    // 1. Make sure buffers are correct size
    fullEvent.resize(numChannels,std::vector<float>(numTicks));
    outputROI.resize(numChannels,std::vector<bool>(numTicks));

    // 2. Copy data to the output array
    for(size_t chanIdx = 0; chanIdx < numChannels; chanIdx++) std::copy(waveform2D[chanIdx].begin(),waveform2D[chanIdx].end(),fullEvent[chanIdx].begin());

    // 3. Apply the frequencyh high pass filters
    std::cout << "==> Step 3: Applying frequrency high pass filters" << std::endl;

    (*fImageFilter)(fullEvent);

    std::cout << "==> Step 4: Perform Canny Edge Detection" << std::endl;

    // Reset the rois array (the resize above may have been a no op)
    Array2D<bool> rois(numChannels,std::vector<bool>(numTicks,false));

    // 4. Apply Canny Edge Detection
    // Since we run on deconvolved waveforms, use dilation 
    fEdgeDetector->Canny(fullEvent, rois, fADFilter_SX, fADFilter_SY, fSigma_x, fSigma_y, fSigma_r, fLowThreshold, fHighThreshold, 'd'); 

    std::cout << "==> Final Step: get dilation, numChannels: " << numChannels << ", rois: " << outputROI.size() << ", output: " << outputROI.size() << std::endl;

    Dilation2D(fBinary_Dilation_SX,fBinary_Dilation_SY)(rois.begin(), numChannels, outputROI.begin());

    std::cout << "==> DONE!! returning to calling module..." << std::endl;

    return;
}

#endif

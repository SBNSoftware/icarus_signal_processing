/**
 * \file ROIFinder2D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class ROIFinder2D
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __icarus_signal_processing_ROIFinder2D_H__
#define __icarus_signal_processing_ROIFinder2D_H__

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <iterator>

#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"
#include "icarus_signal_processing/Filters/ImageFilters.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Detection/EdgeDetection.h"
#include "icarus_signal_processing/Filters/BilateralFilters.h"

namespace icarus_signal_processing
{

/**
 \class ROIFinder2D
 User defined class ROIFinder2D ... these comments are used to generate
 doxygen documentation!
**/

class IROIFinder2D
{
public:
    template <class T> using Array2D = std::vector<std::vector<T>>;
    /**
    *  @brief  Virtual Destructor
    */
    virtual ~IROIFinder2D() noexcept = default;
    /**
    *  @brief Interface functions which provided templated access
    *
    *  @param Waveform  The waveform to process
    */
    virtual void operator()(const Array2D<float>&, Array2D<float>&, Array2D<bool>&) const = 0;
};


class ROIChainFilter : virtual public IROIFinder2D
{
public:
    /// Default constructor
    ROIChainFilter(size_t             FREQUENCY_THRESHOLD,
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
                   const float        GLOBAL_THRESHOLDING_FACTOR);

    /**
     *  @brief  Destructor
     */
    ~ROIChainFilter() {}

    void operator()(
        const IROIFinder2D::Array2D<float>& waveform2D,
        IROIFinder2D::Array2D<float>&       fullEvent,
        IROIFinder2D::Array2D<bool>&        outputROI) const override;

private:
 /*   size_t             fFrequencyThreshold;
    size_t             fFrequency_filter_smoothness_order;
    size_t             fFrequency_filter_mode;
    char               fMorphological_filter_name;
    const unsigned int fChannel_Grouping;
    const unsigned int fStructuring_Element_X;
    const unsigned int fStructuring_Element_Y;
    const unsigned int fROI_Expand_Window_Size;
    const float        fMorphological_Threshold_Factor;
    const size_t       fThetaSteps;
    const unsigned int fHough_Threshold;
    const unsigned int fNMS_Window_Size;
    const unsigned int fAngle_Window;

    const unsigned int fADFilter_SX;
    const unsigned int fADFilter_SY;
    const unsigned int fBinary_Dilation_SX; 
    const unsigned int fBinary_Dilation_SY;
    const float        fGlobal_Thresholding_Factor;
*/
};

//class ROICannyFilter : virtual public IROIFinder2D
//{
//public:
//    /// Default constructor
//    ROICannyFilter(size_t             FREQUENCY_THRESHOLD,
//                   size_t             FREQUENCY_FILTER_SMOOTHNESS_ORDER,
//                   size_t             FREQUENCY_FILTER_MODE,
//                   char               MORPHOLOGICAL_FILTER_NAME,
//                   const unsigned int CHANNEL_GROUPING,
//                   const unsigned int STRUCTURING_ELEMENT_X,
//                   const unsigned int STRUCTURING_ELEMENT_Y,
//                   const unsigned int ROI_EXPAND_WINDOW_SIZE,
//                   const float        MORPHOLOGICAL_THRESHOLD_FACTOR,
//                   const size_t       THETASTEPS,
//                   const unsigned int HOUGH_THRESHOLD,
//                   const unsigned int NMS_WINDOW_SIZE,
//                   const unsigned int ANGLE_WINDOW,
//                   // float NOISE_VARIANCE = 20.0;
//                   const unsigned int ADFILTER_SX,
//                   const unsigned int ADFILTER_SY,
//                   const unsigned int BINARY_CLOSING_SX,
//                   const unsigned int BINARY_CLOSING_SY,
//                   const float        GLOBAL_THRESHOLDING_FACTOR);
//
//    /**
//     *  @brief  Destructor
//     */
//    ~ROICannyFilter() {}
//
//    void operator()(
//        const IROIFinder2D::Array2D<float> &waveform2D,
//        IROIFinder2D::Array2D<float> &fullEvent,
//        IROIFinder2D::Array2D<bool> &outputROI,
//        IROIFinder2D::Array2D<float> &waveLessCoherent,
//        IROIFinder2D::Array2D<float> &morphedWaveform2D,
//        IROIFinder2D::Array2D<float> &finalErosion2D) const override;
//        
//private:
//    size_t             fFrequencyThreshold;
//    size_t             fFrequency_filter_smoothness_order;
//    size_t             fFrequency_filter_mode;
//    char               fMorphological_filter_name;
//    const unsigned int fChannel_Grouping;
//    const unsigned int fStructuring_Element_X;
//    const unsigned int fStructuring_Element_Y;
//    const unsigned int fROI_Expand_Window_Size;
//    const float        fMorphological_Threshold_Factor;
//    const size_t       fThetaSteps;
//    const unsigned int fHough_Threshold;
//    const unsigned int fNMS_Window_Size;
//    const unsigned int fAngle_Window;
//
//    const unsigned int fADFilter_SX;
//    const unsigned int fADFilter_SY;
//    const unsigned int fBinary_Closing_SX; 
//    const unsigned int fBinary_Closing_SY;
//    const float        fGlobal_Thresholding_Factor;
//
//    std::unique_ptr<IFFTFilterFunction> fFilterFunction; ///< The filter function we use for this image filtering
//    std::unique_ptr<ImageFilters>       fImageFilter;    ///< The image filter 
//};

class ROICannyFilter : virtual public IROIFinder2D
{
public:
    /// Default constructor
    ROICannyFilter(const IFFTFilterFunction*,
//                   const IMorphologicalFunctions2D*,
                   const IDenoiser2D*,
                   const BilateralFilters*,
                   const EdgeDetection*,
                   const unsigned int ADFILTER_SX = 7,
                   const unsigned int ADFILTER_SY = 7,
                   const float        sigma_x = 5.0, 
                   const float        sigma_y = 5.0, 
                   const float        sigma_r = 30.0, 
                   const float        lowThreshold = 3.0,
                   const float        highThreshold = 15.0,
                   const unsigned int BINARY_DILATION_SX = 13,
                   const unsigned int BINARY_DILATION_SY = 13);

    /**
     *  @brief  Destructor
     */
    ~ROICannyFilter() {}

    void operator()(
        const IROIFinder2D::Array2D<float>& waveform2D,
        IROIFinder2D::Array2D<float>&       fullEvent,
        IROIFinder2D::Array2D<bool>&        outputROI) const override;
        
private:
//    const IMorphologicalFunctions2D* fMorphologyFilter;
    const IDenoiser2D*               fDenoising;
    std::unique_ptr<ImageFilters>    fImageFilter;    ///< The image filter 
 //   const BilateralFilters*          fBilateralFilter;
    const EdgeDetection*             fEdgeDetector;

    const unsigned int               fADFilter_SX;
    const unsigned int               fADFilter_SY;
    const float                      fSigma_x;
    const float                      fSigma_y;
    const float                      fSigma_r;
    const float                      fLowThreshold;
    const float                      fHighThreshold;
    const unsigned int               fBinary_Dilation_SX; 
    const unsigned int               fBinary_Dilation_SY;
    
};

} // namespace icarus_signal_processing

#endif
/** @} */ // end of doxygen group

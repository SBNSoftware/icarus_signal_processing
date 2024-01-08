/**
 * \file Denoising.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for a class Denoising
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_DENOISING_H__
#define __SIGPROC_TOOLS_DENOISING_H__

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <functional>
#include "icarus_signal_processing/Detection/MorphologicalFunctions1D.h"
#include "icarus_signal_processing/Detection/MorphologicalFunctions2D.h"
#include "ICARUSSigProcDefs.h"
#include "icarus_signal_processing/WaveformTools.h"

namespace icarus_signal_processing {

/**
   \class Denoising
   User defined class Denoising ... these comments are used to generate
   doxygen documentation!
*/
class Denoising 
{
public:
    /// Default constructor
    Denoising(bool outputStats) : fMPVec(4096) {} //fOutputStats(outputStats) {}

    /// Default destructor
    ~Denoising(){}

    void getSelectVals(ArrayFloat::const_iterator,
                       ArrayBool::iterator,
                       ArrayBool::iterator,
                       const VectorFloat&,
                       const unsigned int,
                       const unsigned int,
                       const unsigned int) const;

    void getSelectValsWPCA(ArrayFloat::const_iterator,
                           ArrayBool::iterator,
                           ArrayBool::iterator,
                           const VectorFloat&,
                           const unsigned int,
                           const unsigned int,
                           const unsigned int) const;

    void removeCoherentNoise(ArrayFloat::iterator,
                             ArrayFloat::const_iterator,
                             ArrayFloat::iterator,
                             ArrayBool::const_iterator,
                             ArrayFloat::iterator,
                             const unsigned int,
                             const unsigned int,
                             const unsigned int ) const;

    void removeCoherentNoiseV2(ArrayFloat::iterator,
                               ArrayFloat::const_iterator,
                               ArrayFloat::iterator,
                               ArrayBool::const_iterator,
                               ArrayFloat::iterator,
                               FilterFunctionVec::const_iterator,
                               const unsigned int,
                               const unsigned int ) const;

    void removeCoherentNoise(ArrayFloat::iterator,
                             ArrayFloat::const_iterator,
                             ArrayFloat::iterator,
                             ArrayBool::iterator,
                             ArrayFloat::iterator,
                             const unsigned int,
                             const unsigned int ) const;
    
    float getMedian(        std::vector<float>::iterator, std::vector<float>::iterator) const;
    float getMostProbable(  std::vector<float>::iterator, std::vector<float>::iterator) const;
    float getMode(          std::vector<float>::iterator, std::vector<float>::iterator) const;
    float getIteratedMedian(std::vector<float>::iterator, std::vector<float>::iterator, int&, int&) const;
    float getIteratedMean(  std::vector<float>::iterator, std::vector<float>::iterator, float) const;
    float getTruncatedMean( std::vector<float>::iterator, std::vector<float>::iterator, float) const;

    bool  getPredictedCorrections(const PointCloud<float>&, const WavePoint<float>&, const Eigen::Matrix2f&, float, float, VectorFloat&) const;

private:
    // The code for the most probable calculation will need a std vector
    // We don't wnat to allocated/deallocate each call so have master copy here
    mutable std::vector<int> fMPVec;
   // bool                     fOutputStats;
  
};

class IDenoiser1D
{
public:
   /**
   *  @brief  Virtual Destructor
   */
   virtual ~IDenoiser1D() noexcept = default;
 
   /**
   *  @brief Interface to the 1D denoising
   *
   *  @param Waveform  The waveform to process
   */
    virtual void operator()(ArrayFloat::iterator,
                            ArrayFloat::const_iterator,
                            ArrayFloat::iterator,
                            ArrayFloat::iterator,
                            ArrayBool::iterator,
                            ArrayBool::iterator,
                            ArrayFloat::iterator,
                            FilterFunctionVec::const_iterator,
                            const VectorFloat&,
                            const unsigned int,
                            const unsigned int,
                            const unsigned int,
                            const unsigned int ) const = 0;
};

class IDenoiser2D
{
public:
   /**
   *  @brief  Virtual Destructor
   */
   virtual ~IDenoiser2D() noexcept = default;
 
   /**
   *  @brief Interface to the 1D denoising
   *
   *  @param Waveform  The waveform to process
   */
    virtual void operator()(ArrayFloat::iterator,
                            ArrayFloat::const_iterator,
                            ArrayFloat::iterator,
                            ArrayFloat::iterator,
                            ArrayBool::iterator,
                            ArrayBool::iterator,
                            ArrayFloat::iterator,
                            const unsigned int ) const = 0;
};

class Denoiser1D : virtual public IDenoiser1D, public Denoising {
public:
    /// Default constructor
    Denoiser1D(bool outputStats=false) : Denoising(outputStats), fOutputStats(outputStats) {}

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    FilterFunctionVec::const_iterator,
                    const VectorFloat&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int ) const override;

private:
    bool fOutputStats;

};

class Denoiser1D_Ave : virtual public IDenoiser1D, public Denoising {
public:
    /// Default constructor
    Denoiser1D_Ave(bool outputStats=false) : Denoising(outputStats), fOutputStats(outputStats) {}

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    FilterFunctionVec::const_iterator,
                    const VectorFloat&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int ) const override;

private:
    bool fOutputStats;

};

class Denoiser1D_PCA : virtual public IDenoiser1D, public Denoising {
public:
    /// Default constructor
    Denoiser1D_PCA(bool outputStats=false) : Denoising(outputStats), fOutputStats(outputStats) {}

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    FilterFunctionVec::const_iterator,
                    const VectorFloat&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int ) const override;

private:
    bool fOutputStats;

};

class Denoiser1D_NoCoherent : virtual public IDenoiser1D, public Denoising {
public:
    /// Default constructor
    Denoiser1D_NoCoherent(bool outputStats=false) : Denoising(outputStats), fOutputStats(outputStats) {}

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    FilterFunctionVec::const_iterator,
                    const VectorFloat&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int ) const override;

private:

};

class Denoiser1D_Protect : virtual public IDenoiser1D, public Denoising {
public:
    /// Default constructor
    Denoiser1D_Protect(bool outputStats=false) : Denoising(outputStats), fOutputStats(outputStats) {}

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    FilterFunctionVec::const_iterator,
                    const VectorFloat&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int ) const override;

private:
    bool fOutputStats;

};

class Denoiser2D : virtual public IDenoiser2D, public Denoising {
public:
    /// Default constructor
    Denoiser2D(const IMorphologicalFunctions2D*,         // Filter function to apply for finding protected regions
               const VectorFloat&,                       // Threshold to apply
               unsigned int,                             // Coherent noise grouping (# of channels)
               unsigned int,                             // Window for morphological filter
               bool outputStats=false);                  // If on will activate some timing statistics

    void operator()(ArrayFloat::iterator,                // Output coherent noise corrected waveforms
                    ArrayFloat::const_iterator,          // Input pedestal subtracted waveforms
                    ArrayFloat::iterator,                // Output morphological waveforms
                    ArrayFloat::iterator,                // Output rms of the coherent noise waveforms
                    ArrayBool::iterator,                 // Selected values
                    ArrayBool::iterator,                 // ROI's used in coherent noise subtractio  n
                    ArrayFloat::iterator,                // Medians used in coherent noise correction
                    const unsigned int)                  // Number of channels in the input image
                    const override;

private:
    const IMorphologicalFunctions2D* fFilterFunction;
    // const VectorFloat&               fThresholdVec;
    unsigned int                     fCoherentNoiseGrouping;
    // unsigned int                     fCoherentNoiseGroupingOffset;
    // unsigned int                     fMorphologicalWindow;
    bool                             fOutputStats;
};

class Denoiser2D_Hough : virtual public IDenoiser2D, public Denoising {
public:
    /// Default constructor
    Denoiser2D_Hough(const IMorphologicalFunctions2D*,    // Filter function to apply for finding protected regions
                     const VectorFloat&,                  // Thresholds to apply
                     unsigned int,                        // Coherent noise grouping (# of channels)
                     unsigned int,                        // Channel grouping offset from 0 channel #.
                     unsigned int,                        // Window for morphological filter
                     bool outputStats=false);             // If on will activate some timing statistics

    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    const unsigned int ) const override;

private:
    const IMorphologicalFunctions2D* fFilterFunction;
    const VectorFloat&               fThresholdVec;
    unsigned int                     fCoherentNoiseGrouping;
   //  unsigned int                     fCoherentNoiseGroupingOffset;
    unsigned int                     fMorphologicalWindow;
    bool fOutputStats;

};


class Denoiser2D_RestrictedHough : virtual public IDenoiser2D, public Denoising {
public:

    Denoiser2D_RestrictedHough(const IMorphologicalFunctions2D*,    // Filter Functions
                               const VectorFloat&,                  // Thresholds
                               unsigned int,                        // Coherent Noise Grouping (# of Channels)
                               unsigned int,                        // Channel grouping offset from 0 to channel #
                               unsigned int,                        // Morphological filter window
                               float,                               // Maximum angular deviation from isochronous line.
                               unsigned int thetaSteps=100,         // number of angular spacings (100 is enough usually)
                               unsigned int houghThreshold=500,     // Threshold for determining lines in hough space
                               bool outputStats=false);             // If on, activate timing statistics


    void operator()(ArrayFloat::iterator,
                    ArrayFloat::const_iterator,
                    ArrayFloat::iterator,
                    ArrayFloat::iterator,
                    ArrayBool::iterator,
                    ArrayBool::iterator,
                    ArrayFloat::iterator,
                    const unsigned int ) const override;


private:
    const IMorphologicalFunctions2D* fFilterFunction;
    const VectorFloat&               fThresholdVec;
    unsigned int                     fCoherentNoiseGrouping;
    // unsigned int                     fCoherentNoiseGroupingOffset;
    unsigned int                     fMorphologicalWindow;
    float                            fMaxAngleDev;
    unsigned int                     fThetaSteps;
    int                              fHoughThreshold;
    bool fOutputStats;

};


}

#endif
/** @} */ // end of doxygen group 


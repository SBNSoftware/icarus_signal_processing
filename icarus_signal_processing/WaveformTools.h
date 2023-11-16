///////////////////////////////////////////////////////////////////////
///
/// \file   WaveformTools.h
///
/// \brief  This is the interface class for a tool to aid in the analysis
///         of waveforms.
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWaveformTool_H
#define IWaveformTool_H

#include "ICARUSSigProcDefs.h"

#include <memory>
#include <tuple>

#include <numeric>
#include <cmath>
#include <complex>

// Eigen includes
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/Jacobi"

#include <functional>
#include <fstream>
#include <iomanip>

#include <iostream>

namespace icarus_signal_processing
{

template <class T> using Waveform   = std::vector<T>;
template <class T> using WavePoint  = Eigen::Vector<T,2>;
template <class T> using PointCloud = std::vector<WavePoint<T>>;

template <class T> class WaveformTools
{
public:
    WaveformTools(int = 10);
    ~WaveformTools() {};
   
    using PeakTuple    = std::tuple<size_t,size_t,size_t>;   // first bin, peak bin, last bin
    using PeakTupleVec = std::vector<PeakTuple>;             // The collection of candidate peaks

    void triangleSmooth(     const Waveform<T>&, Waveform<T>&, size_t = 0)                                                                   const;
    void medianSmooth(       const Waveform<T>&, Waveform<T>&, size_t = 3)                                                                   const;
    void averageSmooth(      const Waveform<T>&, Waveform<T>&, size_t = 3)                                                                   const;
    void truncAveSmooth(     const Waveform<T>&, Waveform<T>&, size_t = 3)                                                                   const;
    void principalComponents(Waveform<T>&, Eigen::Vector<T,2>&, Eigen::Matrix<T,2,2>&, Eigen::Vector<T,2>&, T = T(6), T = T(3), bool = true) const;
    void principalComponents(PointCloud<T>&, 
                             Eigen::Vector<T,2>&, 
                             Eigen::Matrix<T,2,2>&, 
                             Eigen::Vector<T,2>&, 
                             T = T(6), T = T(3), 
                             bool = true) const;
    void firstDerivative(    const Waveform<T>&,  Waveform<T>&)                                                                              const;
    void findPeaks(typename Waveform<T>::iterator, typename Waveform<T>::iterator, PeakTupleVec&, T, size_t)                                 const;

    void getErosionDilationAverageDifference(const Waveform<T>&,
                                             int,
                                             Waveform<T>&,
                                             Waveform<T>&,
                                             Waveform<T>&,
                                             Waveform<T>&) const;

    void getOpeningAndClosing(const Waveform<T>&,  const Waveform<T>&,  int, Waveform<T>&,  Waveform<T>&)  const;

    void getMedian(typename Waveform<T>::iterator, typename Waveform<T>::iterator, T&)                     const;

    void getMedian(const PointCloud<T>&, WavePoint<T>&)                                                    const;
    void getMedian(const Waveform<T>&, T&)                                                                 const;
    void getTruncatedMean(const Waveform<T>&, T&, int&, int&)                                              const;
    void getTruncatedRMS(const Waveform<T>&, T, T&, T&, int&)                                              const;
    void getTruncatedMeanRMS(const Waveform<T>&, T, T&, T&, T&, int&, int&)                                const;
    void getPedestalCorrectedWaveform(const Waveform<T>&, Waveform<T>&, T, T&, T&, T&, int&, int&)         const;
   
private:
    bool computePCA(PointCloud<T>&, Eigen::Vector<T,2>&, Eigen::Matrix<T,2,2>&, Eigen::Vector<T,2>&, T&)     const;
 
    int fMinRange;              //< Minimum range to use when calculating the mean
};

    
//----------------------------------------------------------------------
// Constructor.
template <typename T>  inline WaveformTools<T>::WaveformTools(int minRange) :
    fMinRange(minRange)
{
    return;
}

template <typename T>  inline void WaveformTools<T>::triangleSmooth(const std::vector<T>& inputVec, std::vector<T>& smoothVec, size_t lowestBin) const
{
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());

    if (inputVec.size() - lowestBin > 4)
    {
        typename std::vector<T>::iterator       curOutItr = smoothVec.begin();
        typename std::vector<T>::const_iterator curInItr  = inputVec.begin() + lowestBin;

        // Handle the edge cases outside the loop
        *curOutItr++ = (3. * *curInItr + 2. * *(curInItr + 1) + *(curInItr + 2)) / 6.;
        curInItr++;
        *curOutItr++ = (2. * *(curInItr - 1) + 3. * *curInItr + 2. * *(curInItr + 1) + *(curInItr + 2)) / 8.;

        while(curInItr++ != inputVec.end() - 3)
        {
            // Take the weighted average of five consecutive points centered on current point
            T newVal = (*(curInItr - 2) + 2. * *(curInItr - 1) + 3. * *curInItr + 2. * *(curInItr + 1) + *(curInItr + 2)) / 9.;

            *curOutItr++ = newVal;
        }

        // Handle edge cases again
        *curOutItr++ = (*(curInItr - 2) + 2. * *(curInItr - 1) + 3. * *curInItr + 2. * *(curInItr + 1)) / 8.;
        curInItr++;
        *curOutItr = (*(curInItr - 2) + 2. * *(curInItr - 1) + 3. * *curInItr) / 6.;
    }

    return;
}

template <typename T>  inline void WaveformTools<T>::medianSmooth(const Waveform<T>& inputVec, Waveform<T>& smoothVec, size_t nBins) const
{
    // For our purposes, nBins must be odd
    if (nBins % 2 == 0) nBins++;
    
    // Make sure the input vector is right sized
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    // Basic set up
    typename Waveform<T>::const_iterator medianBinItr = inputVec.begin();
    typename Waveform<T>::const_iterator stopItr      = inputVec.end();

    size_t smoothBin(0);

    Waveform<T> medianVec(nBins);
    
    while(std::distance(medianBinItr,stopItr) > 0)
    {
        size_t lowBins = std::min(nBins/2,size_t(std::distance(inputVec.begin(),medianBinItr)));
        size_t hiBins  = std::min(nBins/2,size_t(std::distance(medianBinItr,inputVec.end())));

        typename Waveform<T>::const_iterator startItr = medianBinItr - lowBins;
        typename Waveform<T>::const_iterator endItr   = medianBinItr + hiBins;

        std::copy(startItr, endItr, medianVec.begin());

        T medianVal = T(0);

        getMedian(medianVec.begin(),medianVec.begin()+std::distance(startItr,endItr),medianVal);
        
        smoothVec[smoothBin++] = medianVal;
        
        medianBinItr++;
    }
    
    return;
}

template <typename T>  inline void WaveformTools<T>::averageSmooth(const Waveform<T>& inputVec, Waveform<T>& smoothVec, size_t nBins) const
{
    // Make sure the input vector is right sized
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    // Basic set up
    typename Waveform<T>::const_iterator currentBinItr = inputVec.begin();
    typename Waveform<T>::const_iterator stopItr       = inputVec.end();

    size_t smoothBin(0);
    
    while(std::distance(currentBinItr,stopItr) > 0)
    {
        size_t lowBins = std::min(nBins/2,size_t(std::distance(inputVec.begin(),currentBinItr)));
        size_t hiBins  = std::min(nBins/2,size_t(std::distance(currentBinItr,inputVec.end())));

        typename Waveform<T>::const_iterator startItr = currentBinItr - lowBins;
        typename Waveform<T>::const_iterator endItr   = currentBinItr + hiBins;
        
        smoothVec[smoothBin++] = std::accumulate(startItr,endItr,T(0)) / T(std::distance(startItr,endItr));
        
        currentBinItr++;
    }
    
    return;
}

template <typename T>  inline void WaveformTools<T>::truncAveSmooth(const Waveform<T>& inputVec, Waveform<T>& smoothVec, size_t nBins) const
{
    // Make sure the input vector is right sized
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    // Basic set up
    typename Waveform<T>::const_iterator currentBinItr = inputVec.begin();
    typename Waveform<T>::const_iterator stopItr       = inputVec.end();

    size_t smoothBin(0);
    
    while(std::distance(currentBinItr,stopItr) > 0)
    {
        size_t lowBins = std::min(nBins/2,size_t(std::distance(inputVec.begin(),currentBinItr)));
        size_t hiBins  = std::min(nBins/2,size_t(std::distance(currentBinItr,inputVec.end())));

        Waveform<T> locWaveform(currentBinItr - lowBins, currentBinItr + hiBins);
        
        T    truncMean;
        int  nTrunc;
        int  range;

        getTruncatedMean(locWaveform, truncMean, nTrunc, range);

        smoothVec[smoothBin++] = truncMean;
        
        currentBinItr++;
    }
    
    return;
}

template <typename T> inline void WaveformTools<T>::getMedian(typename Waveform<T>::iterator startItr, 
                                                              typename Waveform<T>::iterator stopItr, 
                                                              T&                             median)  const
{
    median = T(0);

    size_t nVals = std::distance(startItr,stopItr);

    if (nVals > 1)
    {
        size_t middleElem = nVals / 2;

        // Processing slightly different depending on whether even or odd number of values
        if (nVals % 2 == 0)
        {
            std::nth_element(startItr, startItr + middleElem,     stopItr);
            std::nth_element(startItr, startItr + middleElem - 1, stopItr);

            median = (*(startItr + middleElem) + *(startItr + middleElem - 1)) / 2.;
        }
        else
        {
            std::nth_element(startItr, startItr + middleElem, stopItr);

            median = *(startItr + middleElem);
        }
    }

    return;
}


template <typename T> inline void WaveformTools<T>::getMedian(const PointCloud<T>& waveformIn, WavePoint<T>& median) const
{
    // We need to make to "Waveform<T>" objects and get median. 
    Waveform<T> firstVec(waveformIn.size());
    Waveform<T> secondVec(waveformIn.size());

    for(size_t idx = 0; idx < waveformIn.size(); idx++)
    {
        firstVec[idx]  = waveformIn[idx][0];
        secondVec[idx] = waveformIn[idx][1];
    }

    getMedian(firstVec,median.first);
    getMedian(secondVec,median.second);

    return;
}

template <typename T> inline void WaveformTools<T>::getMedian(const Waveform<T>& waveformIn, T& median) const
{
    median = T(0);

    // Require at least 3 for median
    if (waveformIn.size() > 2) 
    {
        Waveform<T> waveform = waveformIn;

        size_t nVals = waveform.size();

        // Differentiate even/odd case for size of vector
        if (nVals % 2 == 0) 
        {
            size_t middleElem = nVals / 2;

            std::nth_element(waveform.begin(), waveform.begin() + middleElem - 1, waveform.end());
            std::nth_element(waveform.begin(), waveform.begin() + middleElem,     waveform.end());

            median = (waveform[middleElem] + waveform[middleElem -1]) / 2;
        } 
        else 
        {
            std::nth_element(waveform.begin(), waveform.begin() + nVals/2, waveform.end());
            median = waveform[nVals/2];
        }
    }
    // Otherwise, if non empty, just average
    else if (!waveformIn.empty())
    {
        median = std::accumulate(waveformIn.begin(),waveformIn.end(),T(0)) / T(waveformIn.size());
    }

    return;
}

template <typename T> inline void WaveformTools<T>::getTruncatedMean(const std::vector<T>& waveformIn, T& mean, int& nTrunc, int& range) const
{
    // Try smoothing first to really find the baseline
    Waveform<T> waveform;

//    medianSmooth(waveformIn, waveform, 7);
    waveform = waveformIn;
    
    // We need to get a reliable estimate of the mean and can't assume the input waveform will be ~zero mean...
    // Basic idea is to find the most probable value in the ROI presented to us
    // From that we can develop an average of the true baseline of the ROI.
    std::pair<typename std::vector<T>::const_iterator,typename std::vector<T>::const_iterator> minMaxValItr = std::minmax_element(waveform.begin(),waveform.end());
    
    T   minVal    = *minMaxValItr.first;
    int minValInt = std::floor(minVal);
    T   maxVal    = *minMaxValItr.second;
    int maxValInt = std::ceil(maxVal);

    range = maxValInt - minValInt + 1;

    if (range > 1)
    {
        std::vector<T> frequencyVec(range, 0);
        int            mpCount(0);
        int            mpVal(0);

        //for(const auto& val : waveform)
        for(size_t waveIdx=0; waveIdx<waveform.size(); waveIdx++)
        {
            T val = waveform[waveIdx];

            if (isnan(val))
            {
                std::cout << "getTruncatedMean finds nan for waveform val! waveform size: " << waveform.size() << ", idx:" << waveIdx << " --> Skipping" << std::endl;
                continue;
            }

            int intVal = std::max(std::min(int(std::round(val - minVal)),range-1),0);

            frequencyVec[intVal]++;

            if (frequencyVec[intVal] > mpCount)
            {
                mpCount = frequencyVec[intVal];
                mpVal   = intVal;
            }
        }

        // take a weighted average of two neighbor bins
        int meanSum  = 0;
        int binRange = std::min(fMinRange, range/2+1);
        int startVal = std::max(0,mpVal-binRange);
        int stopVal  = std::min(range-1,mpVal+binRange);

        nTrunc = 0;

        for(int idx = startVal; idx <= stopVal; idx++)
        {
            meanSum += idx * frequencyVec[idx];
            nTrunc  += frequencyVec[idx];
        }

        mean = T(meanSum) / T(nTrunc) + T(minVal);
    }
    else
    {
        mean   = minVal;
        nTrunc = waveform.size();
    }
    
    return;
}

template <typename T>  inline void WaveformTools<T>::getTruncatedRMS(const std::vector<T>& waveformIn, T nSig, T& rmsFull, T& rmsTrunc, int& nTrunc) const
{
    // This function assumes the input waveform has already been zero-suppressed
    // do rms calculation - the old fashioned way and over all adc values
    // recalculate the rms for truncation
    // Need a local copy of the input waveform (sadly... a copy operation)
    Waveform<T> waveform = waveformIn;

    // Now get the 
    rmsFull = std::inner_product(waveform.begin(), waveform.end(), waveform.begin(), 0.);
    rmsFull = std::sqrt(std::max(T(0.),rmsFull / T(waveform.size())));
    
    // Alternative to sorting is to remove elements outside of a cut range...
    T rmsCut = nSig * rmsFull;
    
    typename std::vector<T>::iterator newEndItr = std::remove_if(waveform.begin(),waveform.end(),[rmsCut](const T& val){return std::abs(val) > rmsCut;});

    rmsTrunc = std::inner_product(waveform.begin(), newEndItr, waveform.begin(), 0.);
    nTrunc   = std::distance(waveform.begin(),newEndItr);
    rmsTrunc = std::sqrt(std::max(T(0.),rmsTrunc / T(nTrunc)));

    return;
}

template <typename T>  inline void WaveformTools<T>::getTruncatedMeanRMS(const std::vector<T>& waveformIn, T nSig, T& mean, T& rmsFull, T& rmsTrunc, int& nTrunc, int& range) const
{
    getTruncatedMean(waveformIn, mean, nTrunc, range);

    std::vector<T> zeroSuppressed(waveformIn.size());

    std::transform(waveformIn.begin(),waveformIn.end(),zeroSuppressed.begin(),std::bind(std::minus<T>(),std::placeholders::_1,mean));

    getTruncatedRMS(zeroSuppressed, nSig, rmsFull, rmsTrunc, nTrunc);
    
    return;
}


template <typename T> inline void WaveformTools<T>::principalComponents(std::vector<T>&       waveform, 
                                                                        Eigen::Vector<T,2>&   meanPos,  
                                                                        Eigen::Matrix<T,2,2>& eigenVectors,
                                                                        Eigen::Vector<T,2>&   eigenValues,
                                                                        T                     adcLimit,
                                                                        T                     nEigenVals,
                                                                        bool                  correctBaseline) const
{
    // Assume we have an input waveform with continuous indices from 0 to size of waveform
    // Step one is to convert to a point cloud
    PointCloud<T> pointCloud(waveform.size());

    for(size_t idx = 0; idx < waveform.size(); idx++) pointCloud[idx] = WavePoint<T>(T(idx),waveform[idx]);

    // Now call function to do PCA
    principalComponents(pointCloud, meanPos, eigenVectors, eigenValues, adcLimit, nEigenVals, correctBaseline);

    // If correcting the baseline then copy back to output
    if (correctBaseline) 
        for(size_t idx = 0; idx < waveform.size(); idx++) waveform[idx] = pointCloud[idx][1];

    return;
}

template <typename T> inline void WaveformTools<T>::principalComponents(PointCloud<T>&        pointCloud, 
                                                                        Eigen::Vector<T,2>&   meanPos,  
                                                                        Eigen::Matrix<T,2,2>& eigenVectors,
                                                                        Eigen::Vector<T,2>&   eigenValues,
                                                                        T                     adcLimit,
                                                                        T                     nEigenVals,
                                                                        bool                  correctBaseline) const
{
    // For first computation insure all ADC values are considered
    T maxADCVal = 5000.;

    if (computePCA(pointCloud, meanPos, eigenVectors, eigenValues, maxADCVal)) 
    {
        // Should we refine to make sure we get a straight baseline?
        if (std::sqrt(eigenValues(0)) > adcLimit)
        {
            maxADCVal = nEigenVals * std::sqrt(eigenValues(0));

            // Note that if the calculation "fails" then you will have the original eigenvectors and values
            computePCA(pointCloud, meanPos, eigenVectors, eigenValues, maxADCVal);
        }
        
        // Correct the baseline?
        if (correctBaseline)
        {
            T yVelocity = eigenVectors.row(1)[1];

            for(size_t idx = 0; idx < pointCloud.size(); idx++)
            {
                T adcCorrection = meanPos[1] + (pointCloud[idx][0] - meanPos[0]) * yVelocity;

                pointCloud[idx][1] -= adcCorrection;
            }
        }
    }

    return;
}

template <typename T>  inline bool WaveformTools<T>::computePCA(PointCloud<T>&        pointCloud,
                                                                Eigen::Vector<T,2>&   meanPos,  
                                                                Eigen::Matrix<T,2,2>& eigenVectors,
                                                                Eigen::Vector<T,2>&   eigenValues,
                                                                T&                    maxADCVal) const
{
    bool goodPCA(false);

    T    yVelocity = eigenVectors.row(1)[1];

    // Start by computing the mean position
    // If the maxADCVal is large then we are not going to reject outliers so this is a simple accumulate
    if (maxADCVal > 1000.)
    {
        // Function to accumulate
        auto pointCloudAdd = [](WavePoint<T> sum, WavePoint<T> val) {sum[0] += val[0]; sum[1] += val[1]; return sum;};

        WavePoint<T> meanVal = std::accumulate(pointCloud.begin(),pointCloud.end(),WavePoint<T>(0,0),pointCloudAdd);

        meanVal[0] /= T(pointCloud.size());
        meanVal[1] /= T(pointCloud.size());

        meanPos = Eigen::Vector<T,2>(meanVal[0], meanVal[1]);
    }
    // Otherwise we do not want to include rejected points in the computation
    else
    {
        WavePoint<T> meanSum(0,0);
        int          nADCVals(0);

        for(size_t waveIdx =0; waveIdx < pointCloud.size(); waveIdx++)
        {
            T x = pointCloud[waveIdx][0] - meanPos(0);
            T y = pointCloud[waveIdx][1] - meanPos(1);

            if (abs(y - x*yVelocity) > maxADCVal) continue;

            meanSum[0] += pointCloud[waveIdx][0];
            meanSum[1] += pointCloud[waveIdx][1];
            nADCVals++;
        }

        meanPos = Eigen::Vector<T,2>(meanSum[0]/T(nADCVals),meanSum[1]/T(nADCVals));
    }

    // Define elements of our covariance matrix
    T xi2(0.);
    T xiyi(0.);
    T yi2(0.0);

    // Reset counter
    int nADCVals(0);

    for (size_t waveIdx = 0; waveIdx < pointCloud.size(); waveIdx++)
    {
        T x = pointCloud[waveIdx][0] - meanPos(0);
        T y = pointCloud[waveIdx][1] - meanPos(1);

        if (abs(y - x*yVelocity) > maxADCVal) continue;

        nADCVals++;

        xi2  += x * x;
        xiyi += x * y;
        yi2  += y * y;
    }

    // Using Eigen package
    Eigen::Matrix<T,2,2> sig;

    sig << xi2, xiyi, xiyi, yi2;

    sig *= 1./(T(nADCVals-1));

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,2,2>> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success) 
    {
        // Now copy output
        eigenVectors = eigenMat.eigenvectors().transpose();
        eigenValues  = eigenMat.eigenvalues();

        goodPCA = true;
    }

    return goodPCA;
}

template <typename T>  inline void WaveformTools<T>::getPedestalCorrectedWaveform(const Waveform<T>& inputWaveform, 
                                                                                  Waveform<T>&       outputWaveform, 
                                                                                  T                  nSig,
                                                                                  T&                 mean,
                                                                                  T&                 rmsFull, 
                                                                                  T&                 rmsTrunc, 
                                                                                  int&               nTrunc,
                                                                                  int&               range) const
{
    // First determine the mean/pedestal based on an average around the most probable value
//    getTruncatedMean(inputWaveform, mean, nTrunc, range);

    Waveform<T> waveform;

    // Try this with the median instead...
    medianSmooth(inputWaveform, waveform, 7);

    std::pair<typename Waveform<T>::const_iterator,typename Waveform<T>::const_iterator> minMaxValItr = std::minmax_element(waveform.begin(),waveform.end());

    range  = std::floor(*minMaxValItr.first) - std::ceil(*minMaxValItr.second) + 1;
    nTrunc = waveform.size();

//    getMedian(waveform, mean);
    getTruncatedMean(waveform,mean,nTrunc,range);
 
    // Do the pedestal correction
    std::transform(inputWaveform.begin(),inputWaveform.end(),outputWaveform.begin(),std::bind(std::minus<T>(),std::placeholders::_1,mean));

    // Now get the RMS values
    getTruncatedRMS(outputWaveform, nSig, rmsFull, rmsTrunc, nTrunc);

    return;
}
    
template <typename T>  inline void WaveformTools<T>::firstDerivative(const std::vector<T>& inputVec, std::vector<T>& derivVec) const
{
    derivVec.resize(inputVec.size(), 0.);
    
    for(size_t idx = 1; idx < derivVec.size() - 1; idx++)
        derivVec.at(idx) = 0.5 * (inputVec.at(idx + 1) - inputVec.at(idx - 1));
    
    return;
}

template <typename T>  inline void WaveformTools<T>::findPeaks(typename std::vector<T>::iterator startItr,
                                                               typename std::vector<T>::iterator stopItr,
                                                               PeakTupleVec&                     peakTupleVec,
                                                               T                                 threshold,
                                                               size_t                            firstTick) const
{
    // Need a minimum distance or else nothing to do
    if (std::distance(startItr,stopItr) > 4)
    {
        // This is a divide and conquer algorithm, start by finding the maximum element.
        typename std::vector<T>::iterator firstItr = std::max_element(startItr,stopItr,[](float left, float right){return std::fabs(left) < std::fabs(right);});

        // Are we over threshold?
        if (std::fabs(*firstItr) > threshold)
        {
            // What am I thinking?
            // First task is to find the "other" lobe max point
            // Set one to the "first", the other to the "second"
            // Search backward from first to find start point, forward from second to find end point
            // Set mid point between first and second as "peak"?
            typename std::vector<T>::iterator secondItr = firstItr;
        
            // Assume if max bin is positive then second lobe is later
            if (*firstItr > 0)
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != stopItr)
                {
                    if (*++tempItr < -threshold)
                    {
                        if (*tempItr < *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            }
            // Otherwise it goes the other way
            else
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != startItr)
                {
                    if (*--tempItr > threshold)
                    {
                        if (*tempItr > *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            
                std::swap(firstItr,secondItr);
            }
        
            // It might that no real pulse was found
            if (firstItr != secondItr)
            {
                // Get the "peak" position
                size_t peakBin = std::distance(startItr,firstItr) + std::distance(firstItr,secondItr) / 2;
        
                // Advance (forward or backward) the first and second iterators to get back to zero crossing
                while(firstItr  != startItr) if (*--firstItr  < 0.) break;
                while(secondItr != stopItr)  if (*++secondItr > 0.) break;
        
                size_t firstBin = std::distance(startItr,firstItr);
                size_t lastBin  = std::distance(startItr,secondItr);
        
                // Find leading peaks
                findPeaks(startItr, firstItr, peakTupleVec, threshold, firstTick);
        
                // Save this peak
                peakTupleVec.push_back(PeakTuple(firstBin+firstTick,peakBin+firstTick,lastBin+firstTick));
        
                // Find downstream peaks
                findPeaks(secondItr, stopItr, peakTupleVec, threshold, firstTick + std::distance(startItr,secondItr));
            }
        }
    }

    return;
}

template <typename T>  inline void WaveformTools<T>::getErosionDilationAverageDifference(const Waveform<T>& inputWaveform,
                                                                                         int                structuringElement,
                                                                                         Waveform<T>&       erosionVec,
                                                                                         Waveform<T>&       dilationVec,
                                                                                         Waveform<T>&       averageVec,
                                                                                         Waveform<T>&       differenceVec) const
{
    // Set the window size
    int halfWindowSize(structuringElement/2);
    
    // Initialize min and max elements
    std::pair<typename Waveform<T>::const_iterator,typename Waveform<T>::const_iterator> minMaxItr =
            std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);
    
    typename Waveform<T>::const_iterator minElementItr = minMaxItr.first;
    typename Waveform<T>::const_iterator maxElementItr = minMaxItr.second;
    
    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    averageVec.resize(inputWaveform.size());
    differenceVec.resize(inputWaveform.size());
    
    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator minItr = erosionVec.begin();
    typename Waveform<T>::iterator maxItr = dilationVec.begin();
    typename Waveform<T>::iterator aveItr = averageVec.begin();
    typename Waveform<T>::iterator difItr = differenceVec.begin();
    
    for (typename Waveform<T>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
        *aveItr++ = 0.5 * (*maxElementItr + *minElementItr);
        *difItr++ = *maxElementItr - *minElementItr;
    }
        
    return;
}

template <typename T>  inline void WaveformTools<T>::getOpeningAndClosing(const Waveform<T>& erosionVec,
                                                                          const Waveform<T>& dilationVec,
                                                                          int                structuringElement,
                                                                          Waveform<T>&       openingVec,
                                                                          Waveform<T>&       closingVec)  const
{
    // Set the window size
    int halfWindowSize(structuringElement/2);
    
    // Start with the opening, here we get the max element in the input erosion vector
    typename Waveform<T>::const_iterator maxElementItr = std::max_element(erosionVec.begin(),erosionVec.begin()+halfWindowSize);
    
    // Initialize the opening vector
    openingVec.resize(erosionVec.size());
    
    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator maxItr = openingVec.begin();
    
    for (typename Waveform<T>::const_iterator inputItr = erosionVec.begin(); inputItr != erosionVec.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,erosionVec.end()) > halfWindowSize)
        {
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *maxItr++ = *maxElementItr;
    }
    
    // Now go with the closling, here we get the min element in the input dilation vector
    typename Waveform<T>::const_iterator minElementItr = std::min_element(dilationVec.begin(),dilationVec.begin()+halfWindowSize);
    
    // Initialize the opening and closing vectors
    closingVec.resize(dilationVec.size());
    
    // Now loop through remaining elements and complete the vectors
    typename Waveform<T>::iterator minItr = closingVec.begin();
    
    for (typename Waveform<T>::const_iterator inputItr = dilationVec.begin(); inputItr != dilationVec.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,dilationVec.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
    }

    return;
}


}

#endif

#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"
#include "icarus_signal_processing/Detection/MorphologicalFunctions2D.h"
#include "icarus_signal_processing/Detection/LineDetection.h"
#include "WaveformTools.h"

#include <chrono>


void icarus_signal_processing::Denoising::getSelectVals(ArrayFloat::const_iterator  morphedWaveformsItr,
                                                        ArrayBool::iterator         selectValsItr,
                                                        ArrayBool::iterator         roiItr,
                                                        const VectorFloat&          thresholdVec,
                                                        const unsigned int          numChannels,
                                                        const unsigned int          grouping,
                                                        const unsigned int          window) const
{
    auto nTicks = morphedWaveformsItr[0].size();

    // Set a protection width
    int halfWidth = std::max(int(window/8),4);

//    std::cout << "***** getSelectVals ***** numChannels: " << numChannels << ", grouping: " << grouping << std::endl;

    for (size_t i=0; i<numChannels; ++i) 
    {
        std::vector<float> localVec = morphedWaveformsItr[i];

        float median = getMedian(localVec.begin(), localVec.end());

        std::vector<float> baseVec;
        baseVec.resize(localVec.size());

        for (size_t j=0; j<baseVec.size(); ++j) baseVec[j] = morphedWaveformsItr[i][j] - median;

//        float rms       = std::sqrt(std::inner_product(baseVec.begin(), baseVec.end(), baseVec.begin(), 0.) / float(baseVec.size()));
//        float threshold = thresholdVec[i / grouping] * rms;
//
//        float maxValue  = *std::max_element(baseVec.begin(),baseVec.end());
//        float minValue  = *std::min_element(baseVec.begin(),baseVec.end());
//
//        std::cout << "   - i: " << i << ", median: " << median << ", rms: " << rms << ", threshold: " << threshold << ", max/min/range: " << maxValue << "/" << minValue << "/" << maxValue-minValue << std::endl;

//        threshold = std::min(threshold,float(12.));
        float threshold = thresholdVec[i / grouping];

        // make sure the roi vector is initialized
        std::fill(roiItr[i].begin(),roiItr[i].end(),false);

        int lb(-2 * halfWidth);

        for (size_t j=0; j<nTicks; ++j) 
        {
            // Are we over threshold?
            if (std::fabs(baseVec[j]) > threshold) 
            {
                // Check Bounds
                selectValsItr[i][j] = true;

                // Check ROI limits
                if (lb < -halfWidth) lb = j - halfWidth;
            } 
            else 
            {
                selectValsItr[i][j] = false;

                // Check if we had gone over threshold and need to set our ROI
                if (lb > -(halfWidth+1))
                {
                    int    ub         = j + halfWidth;
                    size_t lowerBound = std::max(lb, 0);
                    size_t upperBound = std::min(ub, (int) nTicks);

                    for (size_t k=lowerBound; k<upperBound; ++k) roiItr[i][k] = true;
                }

                lb = -2 * halfWidth;
            }
        }
    }

    return;
}

void icarus_signal_processing::Denoiser1D::operator()(ArrayFloat::iterator              waveLessCoherentItr,
                                                      ArrayFloat::const_iterator        filteredWaveformsItr,
                                                      ArrayFloat::iterator              morphedWaveformsItr,
                                                      ArrayFloat::iterator              intrinsicRMSItr,
                                                      ArrayBool::iterator               selectValsItr,
                                                      ArrayBool::iterator               roiItr,
                                                      ArrayFloat::iterator              correctedMediansItr,
                                                      FilterFunctionVec::const_iterator filterFunctionsItr,
                                                      const VectorFloat&                thresholdVec,
                                                      const unsigned int                numChannels,
                                                      const unsigned int                grouping,
                                                      const unsigned int                groupingOffset,
                                                      const unsigned int                window) const
{
//    auto nTicks  = filteredWaveformsItr->size();
//    auto nGroups = numChannels / grouping;
    std::chrono::high_resolution_clock::time_point funcStartTime;
    std::chrono::high_resolution_clock::time_point morphStart;
    std::chrono::high_resolution_clock::time_point morphStop;
    std::chrono::high_resolution_clock::time_point selStart;
    std::chrono::high_resolution_clock::time_point selStop;
    std::chrono::high_resolution_clock::time_point noiseStart;

    if (fOutputStats)
    {
        funcStartTime = std::chrono::high_resolution_clock::now();
        morphStart    = funcStartTime;
    }

    for(size_t funcIdx = 0; funcIdx < numChannels; funcIdx++)
    {
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctionsItr + funcIdx)).get();

        if (!func) 
        {
            std::cout << "Found null function for funcIdx " << funcIdx << " of " << numChannels << std::endl;
            continue;
        }

        (*func)(*(filteredWaveformsItr + funcIdx), *(morphedWaveformsItr + funcIdx));
    }

    if (fOutputStats)
    {
        morphStop  = std::chrono::high_resolution_clock::now();
        selStart   = morphStop;
    }

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdVec, numChannels, grouping, window);

    if (fOutputStats)
    {
        selStop    = std::chrono::high_resolution_clock::now();
        noiseStart = selStop;
    }

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, grouping, groupingOffset);

    if (fOutputStats)
    {
        std::chrono::high_resolution_clock::time_point noiseStop    = std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
    
        std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
        std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
        std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
        std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);

        std::cout << "*** Denoising ***  - # channels: " << numChannels << ", ticks: " << filteredWaveformsItr->size() << ", groups: " << numChannels / grouping << std::endl;
        std::cout << "                   - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
    }
  
    return;
}

void icarus_signal_processing::Denoiser1D_Ave::operator()(ArrayFloat::iterator              waveLessCoherentItr,
                                                          ArrayFloat::const_iterator        filteredWaveformsItr,
                                                          ArrayFloat::iterator              morphedWaveformsItr,
                                                          ArrayFloat::iterator              intrinsicRMSItr,
                                                          ArrayBool::iterator               selectValsItr,
                                                          ArrayBool::iterator               roiItr,
                                                          ArrayFloat::iterator              correctedMediansItr,
                                                          FilterFunctionVec::const_iterator filterFunctionsItr,
                                                          const VectorFloat&                thresholdVec,
                                                          const unsigned int                numChannels,
                                                          const unsigned int                grouping,
                                                          const unsigned int                offset,
                                                          const unsigned int                window) const
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    // Plan here will be form morphed waveforms from the average smoothed waveform of the group
    // To start need a holder for the average waveform...
    VectorFloat aveWaveform(nTicks);
    size_t      channelIdx(0);
    size_t      morphIdx(0);

    // get an instance of the waveform tools
    icarus_signal_processing::WaveformTools<float> waveformTools;
    
    // Initiate the outer loop over channels which will actually be over groups
    while(channelIdx < numChannels)
    {
        std::fill(aveWaveform.begin(),aveWaveform.end(),0.);

        for(size_t innerIdx = 0; innerIdx < grouping; innerIdx++)
        {
            std::transform(aveWaveform.begin(),aveWaveform.end(),(*(filteredWaveformsItr + channelIdx+innerIdx)).begin(), aveWaveform.begin(), std::plus<float>());
        }

        float normFactor(1./float(grouping));

        std::transform(aveWaveform.begin(),aveWaveform.end(),aveWaveform.begin(),std::bind(std::multiplies<float>(),std::placeholders::_1,normFactor));

        // Smooth to try to take out high frequency fluctuations
        waveformTools.medianSmooth(aveWaveform,aveWaveform,5);

        // Now run morphological filter
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctionsItr + channelIdx)).get();

        (*func)(aveWaveform, *(morphedWaveformsItr + morphIdx++));

        // Make sure the channelIdx keeps pace
        channelIdx += grouping;
    }

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdVec, nGroups, grouping, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    VectorFloat v(grouping);
    ArrayFloat  smoothWave(numChannels);

    for (size_t i=0; i<nTicks; ++i) 
    {
        for (size_t j=0; j<nGroups; ++j) 
        {
            size_t group_start = j * grouping;
            size_t group_end = (j+1) * grouping;

            // Compute median.
            size_t idxV(0);

            for (size_t c=group_start; c<group_end; ++c) 
            {
                if (i==0) waveformTools.triangleSmooth(filteredWaveformsItr[c],smoothWave[c]);
                if (!selectValsItr[j][i] && std::abs(filteredWaveformsItr[c][i]) < 20) v[idxV++] = smoothWave[c][i]; //filteredWaveformsItr[c][i];
                //if (!selectValsItr[j][i]) v[idxV++] = smoothWave[c][i]; //filteredWaveformsItr[c][i];
            }

            if (idxV > 0)
            {
                float median   = getMedian(v.begin(),v.begin()+idxV);
//                float mostProb = getMostProbable(v,idxV);
//                float average  = std::accumulate(v.begin(),v.end(),0.) / float(idxV);

                std::transform(v.begin(),v.begin()+idxV,v.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
                float rms = std::sqrt(std::inner_product(v.begin(),v.begin()+idxV,v.begin(),0.) / float(v.size()));

                std::sort(v.begin(),v.begin()+idxV,[](const auto& left,const auto& right){return std::abs(left) < std::abs(right);});

                while(idxV > 0)
                {
                    if (std::abs(v[idxV-1]) < rms) break;
                    idxV--;
                }

                float tempMedian = getMedian(v.begin(),v.begin()+idxV);

//                std::cout << ">tick: " << i << ", group: " << j << ", median: " << median << ", mostProb: " << mostProb << ", average: " << average << ", rms: " << rms << ", tempMedian: " << tempMedian << std::endl;

                median += tempMedian;

                correctedMediansItr[j][i] = median;
                for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
            }
        }
    }

    float rms(0.);
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
            intrinsicRMSItr[i][j] = rms;
        }
    }

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);

    if (fOutputStats)
    {
        std::cout << "*** Denoising ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
        std::cout << "                   - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
    }
  
    return;
}

void icarus_signal_processing::Denoiser1D_Protect::operator()(ArrayFloat::iterator              waveLessCoherentItr,
                                                              ArrayFloat::const_iterator        filteredWaveformsItr,
                                                              ArrayFloat::iterator              morphedWaveformsItr,
                                                              ArrayFloat::iterator              intrinsicRMSItr,
                                                              ArrayBool::iterator               selectValsItr,
                                                              ArrayBool::iterator               roiItr,
                                                              ArrayFloat::iterator              correctedMediansItr,
                                                              FilterFunctionVec::const_iterator filterFunctionsItr,
                                                              const VectorFloat&                thresholdVec,
                                                              const unsigned int                numChannels,
                                                              const unsigned int                grouping,
                                                              const unsigned int                groupingOffset,
                                                              const unsigned int                window) const
{
//    auto nTicks  = filteredWaveformsItr->size();
//    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    for(size_t funcIdx = 0; funcIdx < numChannels; funcIdx++)
    {
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctionsItr + funcIdx)).get();

        if (!func) 
        {
            std::cout << "Found null function for funcIdx " << funcIdx << " of " << numChannels << std::endl;
            continue;
        }

        (*func)(*(filteredWaveformsItr + funcIdx), *(morphedWaveformsItr + funcIdx));
    }

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdVec, numChannels, grouping, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoiseV2(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, filterFunctionsItr, numChannels, grouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);

    if (fOutputStats)
    {
        std::cout << "*** Denoising V2 ***  - # channels: " << numChannels << ", ticks: " << filteredWaveformsItr->size() << ", groups: " << numChannels / grouping << std::endl;
        std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
    }
  
    return;
}

float icarus_signal_processing::Denoising::getMedian(std::vector<float>::iterator vecStart, std::vector<float>::iterator vecEnd) const
{
    float median(0.);

    size_t nVals = std::distance(vecStart,vecEnd);

    if (nVals > 0) 
    {
        if (nVals < 3)
        {
            median = std::accumulate(vecStart,vecEnd,0.) / float(nVals);
        }
        else if (nVals % 2 == 0) 
        {
            const auto m1 = vecStart + nVals / 2 - 1;
            const auto m2 = vecStart + nVals / 2;
            std::nth_element(vecStart, m1, vecEnd);
            const auto e1 = *m1;
            std::nth_element(vecStart, m2, vecEnd);
            const auto e2 = *m2;
            median = (e1 + e2) / 2.0;
        } 
        else 
        {
            const auto m = vecStart + nVals / 2;
            std::nth_element(vecStart, m, vecEnd);
            median = *m;
        }
    }

    return median;
}

float icarus_signal_processing::Denoising::getMostProbable(std::vector<float>::iterator vecStart, std::vector<float>::iterator vecEnd) const
{
    float mostProbable(0.);

    size_t nVals = std::distance(vecStart,vecEnd);

    // Do a simple average if only a few bins
    if (nVals < 5)
    {
        mostProbable = std::accumulate(vecStart,vecStart+nVals,0.) / float(nVals);
    }
    // Otherwise try to form value around MP
    else
    {
        auto minVItr = std::min_element(vecStart,vecStart+nVals);
        auto maxVItr = std::max_element(vecStart,vecStart+nVals);

        float minValue = *minVItr;
        float maxValue = *maxVItr;

        size_t numBins = size_t(maxValue - minValue + 1);

        std::fill(fMPVec.begin(),fMPVec.begin() + numBins,0);
  
        for(std::vector<float>::const_iterator vItr = vecStart; vItr != vecStart + nVals; vItr++)
        {
            int binIdx = int(std::round((*vItr - minValue)));
  
            fMPVec[binIdx]++;
        }
  
        std::vector<int>::iterator maxItr = std::max_element(fMPVec.begin(),fMPVec.begin()+numBins);
  
        int count = *maxItr;
        float   mpVal = float(count) * float(std::distance(fMPVec.begin(),maxItr));

        auto cntItr = maxItr;

        while(cntItr != fMPVec.begin())
        {
            if (*(cntItr - 1) < 0.5 * *maxItr) break;

            count += *(cntItr - 1);
            mpVal += float(*(cntItr - 1)) * float(std::distance(fMPVec.begin(),cntItr - 1));

            cntItr--;
        }

        cntItr = maxItr;

        while(cntItr + 1 != fMPVec.end())
        {
            if (*(cntItr + 1) < 0.5 * *maxItr) break;

            count += *(cntItr + 1);
            mpVal += float(*(cntItr + 1)) * float(std::distance(fMPVec.begin(),cntItr + 1));
            cntItr++;
        }
  
        mostProbable = mpVal / float(count) + minValue;
    }

    return mostProbable;
}

float icarus_signal_processing::Denoising::getMode(std::vector<float>::iterator vecStart, std::vector<float>::iterator vecEnd) const
{
    float mode(0.);

    size_t nVals = std::distance(vecStart,vecEnd);

    // Do a simple average if only a few bins
    if (nVals < 5)
    {
        mode = std::accumulate(vecStart,vecStart+nVals,0.) / float(nVals);
    }
    // Otherwise, try to calculate the mode
    else
    {
        // Note our input data is in float format and one can't really calculate a mode from that
        // Our program is to convert the input vals to integer values where we multiply by two so 
        // we can effectively have numbers on the scale of 0.5 ADC counts. 
        std::map<int,int> occuranceMap;

//        for(const auto& val : vals) occuranceMap[int(std::round(2 * val))]++;
        for(std::vector<float>::const_iterator itr = vecStart; itr != vecEnd; itr++) occuranceMap[int(*itr)]++;

        // Find the value with the most occurances
        int maxOccuranceVal(0);
        int maxOccuranceCnt(-1);

        for(const auto& mapItr : occuranceMap)
        {
            if (mapItr.second > maxOccuranceCnt)
            {
                maxOccuranceVal = mapItr.first;
                maxOccuranceCnt = mapItr.second;
            }
        }

        if (maxOccuranceCnt > 0)
        {
            float modeCnt = maxOccuranceCnt;

            mode = maxOccuranceCnt * float(maxOccuranceVal);

            // Averaging with two nearest neigbors
            if (occuranceMap.find(maxOccuranceVal-1) != occuranceMap.end())
            {
                mode    += occuranceMap[maxOccuranceVal-1] * (maxOccuranceVal - 1);
                modeCnt += occuranceMap[maxOccuranceVal-1];
            }

            if (occuranceMap.find(maxOccuranceVal+1) != occuranceMap.end())
            {
                mode    += occuranceMap[maxOccuranceVal+1] * (maxOccuranceVal + 1);
                modeCnt += occuranceMap[maxOccuranceVal+1];
            }

            mode /= modeCnt;
        }
        else
        {
            std::cout << "****> Not able to compute the mode! But this cannot happen? " << maxOccuranceCnt << std::endl;
        }
    }

    return mode;
}

float icarus_signal_processing::Denoising::getIteratedMedian(std::vector<float>::iterator vecBegin, std::vector<float>::iterator vecEnd, int& range, int& coreRange) const
{
    float median   = 0.;
    size_t numBins = std::distance(vecBegin,vecEnd);

    range     = 5000;
    coreRange = 5000;

    // Realistically, we should have enough bins to get a reliable result, arbitrarily choose that to be 8 bins
    if (numBins > 8)
    {
         // Now sort into ascending order, values at the beginning may be negative, values at the end positive
         std::sort(vecBegin,vecEnd);

         range     = std::round(*(vecEnd-1) - *vecBegin);
         coreRange = std::round(*(vecEnd-4)- *(vecBegin+3));
         median    = getMedian(vecBegin,vecEnd);

         // Make sure we don't have an infinite loop
//         int maxLoops = numBins / 4;

         // Try a refinement by tossing out the largest deviation (which will be at one of the two ends)
//         while(maxLoops-- && range > 2 * coreRange)
         while(std::distance(vecBegin,vecEnd) > 8 && range > 2 * coreRange)
         {
             // Which end is largest deviation? 
             // By definition the median will be more positive than first element, less than the last element
             if (median - *vecBegin > *(vecEnd-1) - median) vecBegin++;
             else                                           vecEnd--;

            range     = std::round(*(vecEnd-1) - *vecBegin);
            coreRange = std::round(*(vecEnd-4)- *(vecBegin+3));
            median    = getMedian(vecBegin,vecEnd);
         }
    }

    return median;
}


icarus_signal_processing::Denoiser2D::Denoiser2D(const IMorphologicalFunctions2D* filterFunction,   // Filter function to apply for finding protected regions
                                                 const VectorFloat&               thresholdVec,     // Threshold to apply
                                                 unsigned int                     coherentGrouping, // Coherent noise grouping (# of channels)
                                                 unsigned int                     morphWindow,      // Window for morphological filter
                                                 bool                             outputStats)      // If on will activate some timing statistics
            : Denoising(outputStats), 
              fFilterFunction(filterFunction),
              // fThresholdVec(thresholdVec),
              fCoherentNoiseGrouping(coherentGrouping),
              // fMorphologicalWindow(morphWindow),
              fOutputStats(outputStats)
{}

void icarus_signal_processing::Denoiser2D::operator()(ArrayFloat::iterator       waveLessCoherentItr,
                                                      ArrayFloat::const_iterator filteredWaveformsItr,
                                                      ArrayFloat::iterator       morphedWaveformsItr,
                                                      ArrayFloat::iterator       intrinsicRMSItr,
                                                      ArrayBool::iterator        selectValsItr,
                                                      ArrayBool::iterator        roiItr,
                                                      ArrayFloat::iterator       correctedMediansItr,
                                                      const unsigned int         numChannels) const
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / fCoherentNoiseGrouping;

 //   std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    // std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*fFilterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

   // std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    // std::chrono::high_resolution_clock::time_point selStart = morphStop;

//    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, fThresholdVec, numChannels, fCoherentNoiseGrouping, fMorphologicalWindow);

    //std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    // std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, fCoherentNoiseGrouping);

    // std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    // std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
   // std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
   // std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
   // std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
   // std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
    if (fOutputStats) std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;

    std::vector<float> v(fCoherentNoiseGrouping);

    for (size_t i=0; i<nTicks; ++i) 
    {
        for (size_t j=0; j<nGroups; ++j) 
        {
            size_t group_start = j * fCoherentNoiseGrouping;
            size_t group_end = (j+1) * fCoherentNoiseGrouping;
            // Compute median.
            size_t idxV(0);

            for (size_t c=group_start; c<group_end; ++c) 
            {
                if (!selectValsItr[c][i]) v[idxV++] = filteredWaveformsItr[c][i];
            }

            float median = getMedian(v.begin(),v.begin()+idxV);

            correctedMediansItr[j][i] = median;
            for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
        }
    }

    float rms(0.);
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            for (size_t k=i*fCoherentNoiseGrouping; k<(i+1)*fCoherentNoiseGrouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
            intrinsicRMSItr[i][j] = rms;
        }
    }
 
//    std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;

    return;
}

icarus_signal_processing::Denoiser2D_Hough::Denoiser2D_Hough(const IMorphologicalFunctions2D* filterFunction,   // Filter function to apply for finding protected regions
                                                             const VectorFloat&               thresholdVec,     // Threshold to apply
                                                             unsigned int                     coherentGrouping, // Coherent noise grouping (# of channels)
                                                             unsigned int                     groupingOffset,   // The collection and middle induction planes are shifted by 23 channels in the beginning.
                                                             unsigned int                     morphWindow,      // Window for morphological filter
                                                             bool                             outputStats)      // If on will activate some timing statistics
            : Denoising(outputStats), 
              fFilterFunction(filterFunction),
              fThresholdVec(thresholdVec),
              fCoherentNoiseGrouping(coherentGrouping),
              // fCoherentNoiseGroupingOffset(groupingOffset),
              fMorphologicalWindow(morphWindow),
              fOutputStats(outputStats)
{}

void icarus_signal_processing::Denoiser2D_Hough::operator()(ArrayFloat::iterator             waveLessCoherentItr,
                                                            ArrayFloat::const_iterator       filteredWaveformsItr,
                                                            ArrayFloat::iterator             morphedWaveformsItr,
                                                            ArrayFloat::iterator             intrinsicRMSItr,
                                                            ArrayBool::iterator              selectValsItr,
                                                            ArrayBool::iterator              roiItr,
                                                            ArrayFloat::iterator             correctedMediansItr,
                                                            const unsigned int               numChannels) const
{
    auto nTicks  = filteredWaveformsItr->size();
//    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*fFilterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    ArrayBool roughSelectVals(numChannels,VectorBool(nTicks));

    getSelectVals(morphedWaveformsItr, roughSelectVals.begin(), roiItr, fThresholdVec, numChannels, fCoherentNoiseGrouping, fMorphologicalWindow);

    ArrayBool localSelectVals = roughSelectVals;

    LineDetection lineModule;

    lineModule.refineSelectVals(roughSelectVals, localSelectVals);

    for(size_t channelIdx = 0; channelIdx < numChannels; channelIdx++)
        *(selectValsItr + channelIdx) = localSelectVals[channelIdx];

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, fCoherentNoiseGrouping, fCoherentNoiseGrouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);

    if (fOutputStats)
    {
        std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << numChannels / fCoherentNoiseGrouping << std::endl;
        std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
    }

    return;
}

icarus_signal_processing::Denoiser2D_RestrictedHough::Denoiser2D_RestrictedHough(const IMorphologicalFunctions2D* filterFunction,   // Filter function to apply for finding protected regions
                                                                                 const VectorFloat&               thresholdVec,     // Threshold to apply
                                                                                 unsigned int                     coherentGrouping, // Coherent noise grouping (# of channels)
                                                                                 unsigned int                     groupingOffset,   // The collection and middle induction planes are shifted by 23 channels in the beginning.
                                                                                 unsigned int                     morphWindow,      // Window for morphological filter
                                                                                 float                            maxAngleDev,      // Maximum angular deviation from isochronous line
                                                                                 unsigned int                     thetaSteps,       // Spacing of angular dimension of accumulator array
                                                                                 unsigned int                     houghThreshold,   // 
                                                                                 bool                             outputStats)      // If on will activate some timing statistics
            : Denoising(outputStats), 
              fFilterFunction(filterFunction),
              fThresholdVec(thresholdVec),
              fCoherentNoiseGrouping(coherentGrouping),
              // fCoherentNoiseGroupingOffset(groupingOffset),
              fMorphologicalWindow(morphWindow),
              fMaxAngleDev(maxAngleDev),
              fThetaSteps(thetaSteps),
              fHoughThreshold(houghThreshold),
              fOutputStats(outputStats)
{}

void icarus_signal_processing::Denoiser2D_RestrictedHough::operator()(ArrayFloat::iterator             waveLessCoherentItr,
                                                                      ArrayFloat::const_iterator       filteredWaveformsItr,
                                                                      ArrayFloat::iterator             morphedWaveformsItr,
                                                                      ArrayFloat::iterator             intrinsicRMSItr,
                                                                      ArrayBool::iterator              selectValsItr,
                                                                      ArrayBool::iterator              roiItr,
                                                                      ArrayFloat::iterator             correctedMediansItr,
                                                                      const unsigned int               numChannels) const
{
    auto nTicks  = filteredWaveformsItr->size();
//    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*fFilterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    ArrayBool roughSelectVals(numChannels,VectorBool(nTicks));

    // Rough protection regions from morphological thresholding
    getSelectVals(morphedWaveformsItr, roughSelectVals.begin(), roiItr, fThresholdVec, numChannels, fCoherentNoiseGrouping, fMorphologicalWindow);

    // Initialize accumulator array for restricted hough transform space.
    ArrayBool localSelectVals(numChannels,VectorBool(nTicks));
    ArrayInt accumulator2D;

    // Vectors to store (i,j) coordinates of the new refined selectVals (signal protection region)
    VectorInt interceptIndex;
    VectorInt angleIndex;

    LineDetection lineModule;

    float maxAngle = (fMaxAngleDev / 180.0) * M_PI;
    float dtheta = 2.0 * maxAngle / (2.0 * (float) fThetaSteps);

    // Run a Cartesian Hough Transform instead of the usual polar form, since we only have to
    // search over the restricted space of near-vertical track signals
    int padding = lineModule.CartesianHoughTransform(roughSelectVals, accumulator2D, fMaxAngleDev, fThetaSteps);

    // Non-maximal suppression gives the local maxima of the accumulator array
    // Use 7 x 7 structuring element for NMS
    lineModule.simpleFastNMS(accumulator2D, interceptIndex, angleIndex, fHoughThreshold, 7, 7);

    // Draw Lines to new refined protection region array
    for (size_t i = 0; i < interceptIndex.size(); ++i) {
        float angle = ((float) angleIndex[i] - (float) fThetaSteps) * dtheta;
        lineModule.drawLine2(localSelectVals, interceptIndex[i], angle, padding);
    }

    icarus_signal_processing::Dilation2D dilationFilter(7,7);
    dilationFilter(localSelectVals.begin(), numChannels, roughSelectVals.begin()); // Store dilation in roughSelectVals

    for(size_t channelIdx = 0; channelIdx < numChannels; channelIdx++)
        *(selectValsItr + channelIdx) = roughSelectVals[channelIdx];

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, fCoherentNoiseGrouping, fCoherentNoiseGrouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);

    if (fOutputStats)
    {
        std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << numChannels / fCoherentNoiseGrouping << std::endl;
        std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
    }

    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise(ArrayFloat::iterator       waveLessCoherentItr,
                                                              ArrayFloat::const_iterator filteredWaveformsItr,
                                                              ArrayFloat::iterator       intrinsicRMSItr,
                                                              ArrayBool::const_iterator  selectValsItr,
                                                              ArrayFloat::iterator       correctedMediansItr,
                                                              const unsigned int         numChannels,
                                                              const unsigned int         grouping,
                                                              const unsigned int         groupingOffset) const
{
    size_t nTicks  = filteredWaveformsItr->size();
    size_t nGroups = (numChannels - groupingOffset) / grouping;   // Assuming numChannels > groupingOffset

    // get an instance of the waveform tools
    icarus_signal_processing::WaveformTools<float> waveformTools;

    VectorFloat vl(grouping/2);

    VectorFloat vu(grouping/2);

    for (size_t i=0; i<nTicks; ++i) 
    {
        for (size_t j=0; j<nGroups; ++j) 
        {
            size_t group_start = j * grouping;
            size_t group_end = (j+1) * grouping;

            size_t group_mid = group_start + grouping/2;

            // Compute median.
            size_t idxL(0);
            size_t idxU(0);

            for(size_t c = group_start; c < group_end; c++)
            {
                // Allow for signal protection?
                // When subtracting in groups of 64 this does not seem to be working as one might expect
                // So removing for now (12/6/2021)
//                if (!selectValsItr[c][i])
//                {
                    if (c < group_mid) vl[idxL++] = filteredWaveformsItr[c][i];
                    else               vu[idxU++] = filteredWaveformsItr[c][i];
//                }
            }

            float median(0.);

            // If we have values which are not "protected" then compute the median
            if (idxL > 8 || idxU > 8)
            {
                float medianL(0.);
                int   rangeL(1000);
                int   coreRangeL(1000);

                if (idxL > 8)
                {
                    // Should we try smoothing? Triangle smoothing requires at least 5 bins
                    if (idxL > 4) 
                    {
                        // Fill out the end of the vector
                        if (idxL < vl.size()) std::fill(vl.begin()+idxL,vl.end(),*(vl.begin()+idxL));

                        VectorFloat tempVec = vl;

                        waveformTools.triangleSmooth(tempVec,vl);
                    }

                    medianL = getIteratedMedian(vl.begin(), vl.begin()+idxL, rangeL, coreRangeL);
                }

                float medianU(0.);
                int   rangeU(1000);
                int   coreRangeU(1000);

                if (idxU > 8)
                {
                   // Should we try smoothing? Triangle smoothing requires at least 5 bins
                    if (idxU > 4) 
                    {
                        // Fill out the end of the vector
                        if (idxU < vu.size()) std::fill(vu.begin()+idxU,vu.end(),*(vl.begin()+idxU));

                        VectorFloat tempVec = vu;

                         waveformTools.triangleSmooth(tempVec,vu);
                    }

                    medianU = getIteratedMedian(vu.begin(), vu.begin()+idxU, rangeU, coreRangeU);
                }

                // Lots of special cases here... for example, if we have "ghost" channels in a grouop then the 
                // range will be zero... so we want to avoid that if possible.
                if (coreRangeU > 1. && coreRangeL > 1.)
                {
                    // Generally, form the "median" as the weighted average between the two groups
                    if (abs(coreRangeL - coreRangeU) < 5) 
                    {
                        float weightL = 1./float(coreRangeL);
                        float weightU = 1./float(coreRangeU);

                        median = (weightL * medianL + weightU * medianU) / (weightL + weightU);
                    }
                    // Otherwise, pick the one from the smallest range
                    else if (coreRangeL < coreRangeU) median = medianL;
                    else                              median = medianU;
                }
                else if (rangeU > 1.) median = medianU;
                else                  median = medianL;
            }
                
            // Add correction
            for (auto k=group_start; k<group_end; ++k) 
            {
                correctedMediansItr[k][i] = median;
                waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
            }
        }
    }

    // Compensate for offset in channel groupings
//    if (groupingOffset > 0) {
//        for (size_t i=0; i<nTicks; ++i) {
//
//            size_t idxV(0);
//
//            for (size_t c=0; c<groupingOffset; ++c)
//            {
//                if (!selectValsItr[c][i]) v[idxV++] = filteredWaveformsItr[c][i];
//            }
//
//            float median(0.);
//
//            // If we have values which are not "protected" then compute the median
//            if (idxV > 0)
//            {
//                std::fill(v.begin()+idxV,v.end(),v.back());
//
//                if (idxV > 3) waveformTools.triangleSmooth(v,v);
//
//                median   = getMedian(v,idxV);
//
//                // Try to improve by throwing out the values at the extremes
//                std::transform(v.begin(),v.begin()+idxV,v.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
//                float rms = std::sqrt(std::inner_product(v.begin(),v.begin()+idxV,v.begin(),0.) / float(idxV));
//
//                std::sort(v.begin(),v.begin()+idxV,[](const auto& left,const auto& right){return std::abs(left) < std::abs(right);});
//
//                while(idxV > 0)
//                {
//                    if (std::abs(v[idxV-1]) < 2.0 * rms) break;
//                    idxV--;
//                }
//
//                // Try to get the improved value for the median. Note we have to add to the previously calculated quantity since it
//                // was subtracted from the vector of values already. 
////                if (idxV > 5) median += std::accumulate(v.begin(),v.begin()+idxV,0.) / float(idxV);
//                if (idxV > 5) median += getMedian(v,idxV);
//            }
//
//            for (unsigned int k=0; k<groupingOffset; ++k) 
//            {
//                correctedMediansItr[k][i] = median;
//                waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
//            }
//
//        }
//    }

    // Now compute the rms for the corrected waveforms
    float rms(0.);
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
            std::vector<float> v(grouping,0.);
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(idxV));
            intrinsicRMSItr[i][j] = rms;
        }
    }
  
    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoiseV2(ArrayFloat::iterator              waveLessCoherentItr,
                                                                ArrayFloat::const_iterator        filteredWaveformsItr,
                                                                ArrayFloat::iterator              intrinsicRMSItr,
                                                                ArrayBool::const_iterator         selectValsItr,
                                                                ArrayFloat::iterator              correctedMediansItr,
                                                                FilterFunctionVec::const_iterator filterFunctionsItr,
                                                                const unsigned int                numChannels,
                                                                const unsigned int                grouping) const
{
    size_t nTicks  = filteredWaveformsItr->size();
    size_t nGroups = (numChannels) / grouping;   

    // get an instance of the waveform tools
    icarus_signal_processing::WaveformTools<float> waveformTools;

    VectorFloat vl(grouping/2);

    VectorFloat vu(grouping/2);

    for (size_t groupIdx=0; groupIdx<nGroups; ++groupIdx) 
    {
        size_t group_start = groupIdx * grouping;
        size_t group_end = (groupIdx+1) * grouping;

        size_t group_mid = group_start + grouping/2;

        for (size_t tickIdx=0; tickIdx<nTicks; ++tickIdx) 
        {
            // Compute median.
            size_t idxL(0);
            size_t idxU(0);

            for(size_t c = group_start; c < group_end; c++)
            {
                if (c < group_mid) vl[idxL++] = filteredWaveformsItr[c][tickIdx];
                else               vu[idxU++] = filteredWaveformsItr[c][tickIdx];
            }

            float median(0.);

            // If we have values which are not "protected" then compute the median
            if (idxL > 8 || idxU > 8)
            {
                float medianL(0.);
                int   rangeL(1000);
                int   coreRangeL(1000);

                if (idxL > 8)
                {
                    // Should we try smoothing? Triangle smoothing requires at least 5 bins
                    if (idxL > 4) 
                    {
                        // Fill out the end of the vector
                        if (idxL < vl.size()) std::fill(vl.begin()+idxL,vl.end(),*(vl.begin()+idxL));

                        VectorFloat tempVec = vl;

                        waveformTools.triangleSmooth(tempVec,vl);
                    }

                    medianL = getIteratedMedian(vl.begin(), vl.begin()+idxL, rangeL, coreRangeL);
                }

                float medianU(0.);
                int   rangeU(1000);
                int   coreRangeU(1000);

                if (idxU > 8)
                {
                   // Should we try smoothing? Triangle smoothing requires at least 5 bins
                    if (idxU > 4) 
                    {
                        // Fill out the end of the vector
                        if (idxU < vu.size()) std::fill(vu.begin()+idxU,vu.end(),*(vl.begin()+idxU));

                        VectorFloat tempVec = vu;

                        waveformTools.triangleSmooth(tempVec,vu);
                    }

                    medianU = getIteratedMedian(vu.begin(), vu.begin()+idxU, rangeU, coreRangeU);
                }

                // Lots of special cases here... for example, if we have "ghost" channels in a grouop then the 
                // range will be zero... so we want to avoid that if possible.
                if (rangeU > 2. && rangeL > 2.)
                {
                    // Generally, form the "median" as the weighted average between the two groups
                    if (abs(coreRangeL - coreRangeU) < 5) 
                    {
                        float weightL = 1./float(coreRangeL);
                        float weightU = 1./float(coreRangeU);

                        median = (weightL * medianL + weightU * medianU) / (weightL + weightU);
                    }
                    // Otherwise, pick the one from the smallest range
                    else if (coreRangeL < coreRangeU) median = medianL;
                    else                              median = medianU;
                }
                // In this case use the larger
                else if (rangeU > rangeL) median = medianU;
                else                      median = medianL;
            }
                
            // Add correction
            for (auto k=group_start; k<group_end; ++k) correctedMediansItr[k][tickIdx] = median;
        }

        // ***********
        // Now can search correction waveform for potential signal. Only need to look at one of them
        const icarus_signal_processing::IMorphologicalFunctions1D* func = (*(filterFunctionsItr + group_start)).get();

        if (!func) 
        {
            std::cout << "Found null function for funcIdx " << group_start << " of " << numChannels << std::endl;
            continue;
        }

        VectorFloat morphedWaveform(nTicks,0.);

//        (*func)(correctedMediansItr[group_start], morphedWaveform);
        std::copy(correctedMediansItr[group_start].begin(),correctedMediansItr[group_start].end(),morphedWaveform.begin());

//        float median = getMedian(morphedWaveform.begin(), morphedWaveform.end());
//
//        // Subtract the median and find the rms
//        std::transform(morphedWaveform.begin(),morphedWaveform.end(),morphedWaveform.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
//        float rms = std::sqrt(std::inner_product(morphedWaveform.begin(),morphedWaveform.end(),morphedWaveform.begin(),0.) / float(nTicks));

//        // Go through and search for potential signal to protect
//        float  hitThreshold(1000.);
//        size_t tickIdx(0);
//
//        while(tickIdx < nTicks)
//        {
//            // Might there be some signal in the correction?
//            if (std::abs(morphedWaveform[tickIdx]) > hitThreshold) //3. * rms)
//            {
//                // Set the stop index since we'll be modifying it
//                size_t startIdx(tickIdx);
//                size_t stopIdx(tickIdx);
//
//                while(stopIdx < nTicks && std::abs(morphedWaveform[stopIdx]) > hitThreshold) stopIdx++;
//
//                if (stopIdx - startIdx > 3)
//                {
//                    // The goal is to move start/stop to the potential outside edges of the found region
//                    // But this will depend on whether the peak found is positive or negative... 
//                    if (morphedWaveform[startIdx] > 0.)
//                    {
//                        // back the start up to find zero or plateau
//                        while(startIdx > 0 && morphedWaveform[startIdx] > 0 && morphedWaveform[startIdx] > morphedWaveform[startIdx-1]) startIdx--;
//
//                        // Similarly, move stop ahead
//                        while(stopIdx < nTicks-1 && morphedWaveform[stopIdx] > 0 && morphedWaveform[stopIdx] > morphedWaveform[stopIdx+1]) stopIdx++;
//                    }
//                    else
//                    {
//                        // back the start up to find zero or plateau
//                        while(startIdx > 0 && morphedWaveform[startIdx] < 0 && morphedWaveform[startIdx] < morphedWaveform[startIdx-1]) startIdx--;
//
//                        // Similarly, move stop ahead
//                        while(stopIdx < nTicks-1 && morphedWaveform[stopIdx] < 0 && morphedWaveform[stopIdx] < morphedWaveform[stopIdx+1]) stopIdx++;
//                    }
//
//                    // Now we require that the width is "pulse like" (meaning more than a simple spike)
//                    if (stopIdx - startIdx > 10)
//                    {
//                        float startBinVal  = morphedWaveform[startIdx];
//                        float morphedSlope = (morphedWaveform[stopIdx] - startBinVal) / float(stopIdx - startIdx);
//
//                        for(size_t corTickIdx=startIdx; corTickIdx<stopIdx; corTickIdx++)
//                        {
//                            for(size_t wireIdx=group_start; wireIdx<group_end; wireIdx++) correctedMediansItr[wireIdx][corTickIdx] = startBinVal + float(corTickIdx - startIdx) * morphedSlope;
//
////                            std::cout << "Found possible signal with value: " << morphedWaveform[((startIdx+stopIdx)/2)] << ", rms: " << rms << ", start/stop: " << startIdx << "/" << stopIdx << ", group: " << group_start << "/" << group_end << std::endl;
//                            std::cout << "Found possible signal with value: " << morphedWaveform[corTickIdx] << ", tickIdx: " << corTickIdx << ", group: " << group_start << "/" << group_end << std::endl;
//                        }
//                    }
//                }
//
//                tickIdx = stopIdx;
//            }
//
//            tickIdx++;
//        }

        // ***********
        // Ok, now can apply correction
        for(size_t k=group_start; k<group_end; k++)
            std::transform(filteredWaveformsItr[k].begin(), filteredWaveformsItr[k].end(), correctedMediansItr[k].begin(), waveLessCoherentItr[k].begin(), std::minus<float>());

        // And now compute the intrinsic rms
        for(size_t tickIdx=0; tickIdx<nTicks; tickIdx++)
        {
            size_t idxV(0);
            std::vector<float> v(grouping,0.);
            for (size_t k=group_start; k<group_end; ++k) v[idxV++] = waveLessCoherentItr[k][tickIdx];
            float rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(idxV));
            intrinsicRMSItr[groupIdx][tickIdx] = rms;
        }
    }
  
    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise(ArrayFloat::iterator       waveLessCoherentItr,
                                                              ArrayFloat::const_iterator filteredWaveformsItr,
                                                              ArrayFloat::iterator       intrinsicRMSItr,
                                                              ArrayBool::iterator        selectValsItr,
                                                              ArrayFloat::iterator       correctedMediansItr,
                                                              const unsigned int         numChannels,
                                                              const unsigned int         grouping) const
{
    auto   nTicks  = filteredWaveformsItr->size();
    size_t nGroups = numChannels / grouping;

    // get an instance of the waveform tools
    icarus_signal_processing::WaveformTools<float> waveformTools;

    VectorFloat v(grouping);

    for (size_t i=0; i<nTicks; ++i) 
    {
        for (size_t j=0; j<nGroups; ++j) 
        {
            size_t group_start = j * grouping;
            size_t group_end = (j+1) * grouping;

            // Compute median.
            size_t idxV(0);

            for (size_t c=group_start; c<group_end; ++c) 
            {
                if (!selectValsItr[c][i]) v[idxV++] = filteredWaveformsItr[c][i];
            }

            float median(0.);

            // If we have values which are not "protected" then compute the median
            if (idxV > 0)
            {
                std::fill(v.begin()+idxV,v.end(),v.back());

                if (idxV > 3) waveformTools.triangleSmooth(v,v);

                median   = getMedian(v.begin(),v.begin()+idxV);

                // Try to improve by throwing out the values at the extremes
                std::transform(v.begin(),v.begin()+idxV,v.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
                float rms = std::sqrt(std::inner_product(v.begin(),v.begin()+idxV,v.begin(),0.) / float(idxV));

                std::sort(v.begin(),v.begin()+idxV,[](const auto& left,const auto& right){return std::abs(left) < std::abs(right);});

                while(idxV > 0)
                {
                    if (std::abs(v[idxV-1]) < 1.75 * rms) break;
                    idxV--;
                }

                // Try to get the improved value for the median. Note we have to add to the previously calculated quantity since it
                // was subtracted from the vector of values already. 
                if (idxV > 5) median += getMedian(v.begin(),v.begin()+idxV);
            }

            for (size_t k = group_start; k < group_end; k++) correctedMediansItr[k][i] = median;
        }
    }

    // Now compute the rms for the corrected waveforms
    for (size_t i=0; i<nGroups; ++i) 
    {
        for (size_t j=0; j<nTicks; ++j) 
        {
            size_t idxV(0);
                
            // Add correction
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) 
            {
                float median = correctedMediansItr[k][j];

                waveLessCoherentItr[k][j] = filteredWaveformsItr[k][j] - median;
                v[idxV++]                 = waveLessCoherentItr[k][j];
            }

            intrinsicRMSItr[i][j] = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
        }
    }
  
    return;
}

#endif

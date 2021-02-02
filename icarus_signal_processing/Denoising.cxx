#ifndef __SIGPROC_TOOLS_DENOISING_CXX__
#define __SIGPROC_TOOLS_DENOISING_CXX__

#include "Denoising.h"
#include "icarus_signal_processing/Detection/LineDetection.h"
#include "WaveformTools.h"

#include <chrono>


void icarus_signal_processing::Denoising::getSelectVals(ArrayFloat::const_iterator  morphedWaveformsItr,
                                                        ArrayBool::iterator         selectValsItr,
                                                        ArrayBool::iterator         roiItr,
                                                        VectorFloat::const_iterator thresholdItr,
                                                        const unsigned int          numChannels,
                                                        const unsigned int          window)
{
    auto nTicks = morphedWaveformsItr[0].size();

    // Set a protection width
    int halfWidth = std::max(int(window/8),1);

    for (size_t i=0; i<numChannels; ++i) 
    {
        std::vector<float> localVec = morphedWaveformsItr[i];

        float median = getMedian(localVec, localVec.size());

        std::vector<float> baseVec;
        baseVec.resize(localVec.size());

        for (size_t j=0; j<baseVec.size(); ++j) baseVec[j] = morphedWaveformsItr[i][j] - median;

        float rms       = std::sqrt(std::inner_product(baseVec.begin(), baseVec.end(), baseVec.begin(), 0.) / float(baseVec.size()));
        float threshold = (*(thresholdItr + i)) * rms;

//        threshold = std::min(threshold,float(12.));

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

void icarus_signal_processing::Denoising::removeCoherentNoise1D(ArrayFloat::iterator              waveLessCoherentItr,
                                                                ArrayFloat::const_iterator        filteredWaveformsItr,
                                                                ArrayFloat::iterator              morphedWaveformsItr,
                                                                ArrayFloat::iterator              intrinsicRMSItr,
                                                                ArrayBool::iterator               selectValsItr,
                                                                ArrayBool::iterator               roiItr,
                                                                ArrayFloat::iterator              correctedMediansItr,
                                                                FilterFunctionVec::const_iterator filterFunctionsItr,
                                                                VectorFloat::const_iterator       thresholdItr,
                                                                const unsigned int                numChannels,
                                                                const unsigned int                grouping,
                                                                const unsigned int                window)
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

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdItr, numChannels, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, grouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
//    std::cout << "*** Denoising ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
//    std::cout << "                   - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
  
    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise1D_Ave(ArrayFloat::iterator              waveLessCoherentItr,
                                                                    ArrayFloat::const_iterator        filteredWaveformsItr,
                                                                    ArrayFloat::iterator              morphedWaveformsItr,
                                                                    ArrayFloat::iterator              intrinsicRMSItr,
                                                                    ArrayBool::iterator               selectValsItr,
                                                                    ArrayBool::iterator               roiItr,
                                                                    ArrayFloat::iterator              correctedMediansItr,
                                                                    FilterFunctionVec::const_iterator filterFunctionsItr,
                                                                    VectorFloat::const_iterator       thresholdItr,
                                                                    const unsigned int                numChannels,
                                                                    const unsigned int                grouping,
                                                                    const unsigned int                window)
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

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdItr, nGroups, window);

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
                float median   = getMedian(v,idxV);
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

                float tempMedian = getMedian(v,idxV);

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
  
//    std::cout << "*** Denoising ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
//    std::cout << "                   - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;
  
    return;
}

float icarus_signal_processing::Denoising::getMedian(std::vector<float>& vals, const unsigned int nVals) const
{
    float median(0.);

    if (nVals > 0) 
    {
        if (nVals < 3)
        {
            median = std::accumulate(vals.begin(),vals.begin()+nVals,0.) / float(nVals);
        }
        else if (nVals % 2 == 0) 
        {
            const auto m1 = vals.begin() + nVals / 2 - 1;
            const auto m2 = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m1, vals.begin() + nVals);
            const auto e1 = *m1;
            std::nth_element(vals.begin(), m2, vals.begin() + nVals);
            const auto e2 = *m2;
            median = (e1 + e2) / 2.0;
        } 
        else 
        {
            const auto m = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m, vals.begin() + nVals);
            median = *m;
        }
    }

    return median;
}

float icarus_signal_processing::Denoising::getMostProbable(std::vector<float>& vals, const unsigned int nVals)
{
    float mostProbable(0.);

    // Do a simple average if only a few bins
    if (nVals < 5)
    {
        mostProbable = std::accumulate(vals.begin(),vals.begin()+nVals,0.) / float(nVals);
    }
    // Otherwise try to form value around MP
    else
    {
        auto minVItr = std::min_element(vals.begin(),vals.begin()+nVals);
        auto maxVItr = std::max_element(vals.begin(),vals.begin()+nVals);

        float minValue = *minVItr;
        float maxValue = *maxVItr;

        size_t numBins = size_t(maxValue - minValue + 1);

        std::fill(fMPVec.begin(),fMPVec.begin() + numBins,0);
  
        for(typename std::vector<float>::iterator vItr = vals.begin(); vItr != vals.begin() + nVals; vItr++)
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

void icarus_signal_processing::Denoising::removeCoherentNoise2D(ArrayFloat::iterator                                       waveLessCoherentItr,
                                                                ArrayFloat::const_iterator                                 filteredWaveformsItr,
                                                                ArrayFloat::iterator                                       morphedWaveformsItr,
                                                                ArrayFloat::iterator                                       intrinsicRMSItr,
                                                                ArrayBool::iterator                                        selectValsItr,
                                                                ArrayBool::iterator                                        roiItr,
                                                                ArrayFloat::iterator                                       correctedMediansItr,
                                                                const icarus_signal_processing::IMorphologicalFunctions2D* filterFunction,
                                                                VectorFloat::const_iterator                                thresholdItr,
                                                                const unsigned int                                         numChannels,
                                                                const unsigned int                                         grouping,
                                                                const unsigned int                                         window)
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*filterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    getSelectVals(morphedWaveformsItr, selectValsItr, roiItr, thresholdItr, numChannels, window);

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, grouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
//    std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;

    std::vector<float> v(grouping);

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

            float median = getMedian(v,idxV);

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
            for (size_t k=i*grouping; k<(i+1)*grouping; ++k) v[idxV++] = waveLessCoherentItr[k][j];
            rms = std::sqrt(std::inner_product(v.begin(), v.begin()+idxV, v.begin(), 0.) / float(v.size()));
            intrinsicRMSItr[i][j] = rms;
        }
    }
 
//    std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;

    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoiseHough(ArrayFloat::iterator                                       waveLessCoherentItr,
                                                                   ArrayFloat::const_iterator                                 filteredWaveformsItr,
                                                                   ArrayFloat::iterator                                       morphedWaveformsItr,
                                                                   ArrayFloat::iterator                                       intrinsicRMSItr,
                                                                   ArrayBool::iterator                                        selectValsItr,
                                                                   ArrayBool::iterator                                        roiItr,
                                                                   ArrayFloat::iterator                                       correctedMediansItr,
                                                                   const icarus_signal_processing::IMorphologicalFunctions2D* filterFunction,
                                                                   VectorFloat::const_iterator                                thresholdItr,
                                                                   const unsigned int                                         numChannels,
                                                                   const unsigned int                                         grouping,
                                                                   const unsigned int                                         window)
{
    auto nTicks  = filteredWaveformsItr->size();
//    auto nGroups = numChannels / grouping;

    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    std::chrono::high_resolution_clock::time_point morphStart = funcStartTime;

    (*filterFunction)(filteredWaveformsItr, numChannels, morphedWaveformsItr);

    std::chrono::high_resolution_clock::time_point morphStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point selStart = morphStop;

    ArrayBool roughSelectVals(numChannels,VectorBool(nTicks));

    getSelectVals(morphedWaveformsItr, roughSelectVals.begin(), roiItr, thresholdItr, numChannels, window);

    ArrayBool localSelectVals = roughSelectVals;

    LineDetection lineModule;

    lineModule.refineSelectVals(roughSelectVals, localSelectVals);

    for(size_t channelIdx = 0; channelIdx < numChannels; channelIdx++)
        *(selectValsItr + channelIdx) = localSelectVals[channelIdx];

    std::chrono::high_resolution_clock::time_point selStop  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point noiseStart = selStop;

    removeCoherentNoise(waveLessCoherentItr, filteredWaveformsItr, intrinsicRMSItr, selectValsItr, correctedMediansItr, numChannels, grouping);

    std::chrono::high_resolution_clock::time_point noiseStop = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> funcTime   = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> morphTime  = std::chrono::duration_cast<std::chrono::duration<double>>(morphStop - morphStart);
    std::chrono::duration<double> selTime    = std::chrono::duration_cast<std::chrono::duration<double>>(selStop - selStart);
    std::chrono::duration<double> noiseTime  = std::chrono::duration_cast<std::chrono::duration<double>>(noiseStop - noiseStart);
  
//    std::cout << "*** Denoising 2D ***  - # channels: " << numChannels << ", ticks: " << nTicks << ", groups: " << nGroups << std::endl;
//    std::cout << "                      - morph: " << morphTime.count() << ", sel: " << selTime.count() << ", noise: " << noiseTime.count() << ", total: " << funcTime.count() << std::endl;

    return;
}

void icarus_signal_processing::Denoising::removeCoherentNoise(ArrayFloat::iterator       waveLessCoherentItr,
                                                              ArrayFloat::const_iterator filteredWaveformsItr,
                                                              ArrayFloat::iterator       intrinsicRMSItr,
                                                              ArrayBool::iterator        selectValsItr,
                                                              ArrayFloat::iterator       correctedMediansItr,
                                                              const unsigned int         numChannels,
                                                              const unsigned int         grouping)
{
    auto nTicks  = filteredWaveformsItr->size();
    auto nGroups = numChannels / grouping;

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

                median   = getMedian(v,idxV);

                // Try to improve by throwing out the values at the extremes
                std::transform(v.begin(),v.begin()+idxV,v.begin(),std::bind(std::minus<float>(),std::placeholders::_1,median));
                float rms = std::sqrt(std::inner_product(v.begin(),v.begin()+idxV,v.begin(),0.) / float(v.size()));

                std::sort(v.begin(),v.begin()+idxV,[](const auto& left,const auto& right){return std::abs(left) < std::abs(right);});

                while(idxV > 0)
                {
                    if (std::abs(v[idxV-1]) < 2.0 * rms) break;
                    idxV--;
                }

                // Try to get the improved value for the median. Note we have to add to the previously calculated quantity since it
                // was subtracted from the vector of values already. 
//                if (idxV > 5) median += std::accumulate(v.begin(),v.begin()+idxV,0.) / float(idxV);
                if (idxV > 5) median += getMedian(v,idxV);
            }
                
            // Add correction
            correctedMediansItr[j][i] = median;
            for (auto k=group_start; k<group_end; ++k) waveLessCoherentItr[k][i] = filteredWaveformsItr[k][i] - median;
        }
    }

    // Now compute the rms for the corrected waveforms
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
  
    return;
}

#endif

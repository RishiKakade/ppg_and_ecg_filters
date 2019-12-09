import numpy as np
import scipy.signal as sg
import scipy as sp
import pandas as pd
import os

# note difference bwtween filtfilt and lfilter
# filtfilt runs filter both forward and backwars resulting in 0 phase response
# is more suited for prerecorded data
# lfilter more resembles the behaivour of electronic filters, applying a filter forward in time
# and causing phase shift in signal

# these filteres are useful for reducing EMG noise and baseline wander, notch for power line interference

def cheby_lowpass(wp, ws, fs, gpass, gstop):
    wp = wp/fs
    ws = ws/fs
    order, wn = sg.cheb2ord(wp, ws, gpass, gstop)
    b, a = sg.cheby2(order, gstop, wn)
    return b, a

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
    
def cheby_lowpass_filter(data, cutoff, fs, gpass, gstop):
    b, a = sg.cheby_lowpass(cutoff[0], cutoff[1], fs, gpass, gstop)
    y = sg.lfilter(b, a, data)
    return y

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = sg.butter_lowpass(cutoff, fs, order=order)
    y = sg.lfilter(b, a, data)
    return y

def highpass_filter(data, cutoff, order=2):
    nyq = 0.5 * 250
    cutoff = cutoff / nyq
    b, a = sg.butter(order, cutoff, btype='highpass', analog=False)
    return sg.filtfilt(b, a, data.tolist())

def bandpass_filter(data, cutoff_low, cutoff_high, order=2):
    nyq = 0.5 * 250 # the samplin rate for my sensor was 250 hz
    normal_cutoff_low = cutoff_low / nyq
    normal_cutoff_high = cutoff_high / nyq
    b, a = sg.butter(order, [normal_cutoff_low, normal_cutoff_high], btype='bandpass', analog=False)
    return sg.filtfilt(b, a, data.tolist())

#QRS peak Detection

def r_peaks(time,ecg):
    r_max=max(ecg)
    r_max_index=[]
    r_peak_array=[]
    x_max=[]
    y_max=[] 

    for i in range(len(ecg)):
        if ecg[i]>r_max*0.5: #0.5 is depended on situation.
            r_max_index.append(i)

    tempIndex=[r_max_index[0]]
    for i in range(1,len(r_max_index)+1):

        if i<len(r_max_index):
            if r_max_index[i]-r_max_index[i-1]<50:
                tempIndex.append(r_max_index[i])
                continue
            elif len(tempIndex)==0:
                tempIndex.append(r_max_index[i])

        tempECG=[]
        for ind in tempIndex:
            tempECG.append(ecg[ind])

        tempRPeakIndex=tempECG.index(max(tempECG))
        r_peak_array.append(tempIndex[tempRPeakIndex])
        tempIndex=[]

    for j in r_peak_array:
        x_max.append(np.array(time)[j])
        y_max.append(np.array(ecg)[j])

    return x_max,y_max,r_peak_array

def findQS(time,ecg,r_peak_array):
    fs=250
    l=len(ecg)
    qx_min=[]
    qy_min=[]
    sx_min=[]
    sy_min=[]

    for ind in r_peak_array:

        w1=fs/10
        w2=fs/5
        if ind<w1:
            w1=ind
        if ind < w2:
            w2=ind
        if ind+w1>l:
            w1=l-ind-1
        if ind+w2>l:
            w2=l-ind-1
        w1TempQ=[]
        w1TempS=[]
        w2TempQ=[]
        w2TempS=[]
        for i in range(int(w1)):

            w1TempQ.append([ecg[ind-i],ind-i])
            w1TempS.append([ecg[ind+i],ind+i])

        for j in range(int(w2)):

            w2TempQ.append([ecg[ind-i],ind-i])
            w2TempS.append([ecg[ind+i],ind+i])

        w1Qmin=min(w1TempQ)
        w1Smin=min(w1TempS)
        w2Qmin=min(w2TempQ)
        w2Smin=min(w2TempS)

        wQmin=min([w1Qmin,w2Qmin])
        qx_min.append(time[wQmin[1]])
        qy_min.append(ecg[wQmin[1]])
        QIndex.append(wQmin[1])

        wSmin=min([w1Smin,w2Smin])
        sx_min.append(time[wSmin[1]])
        sy_min.append(ecg[wSmin[1]])
        SIndex.append(wSmin[1])

    return qx_min,qy_min,sx_min,sy_min


# removes datapoints beyond saturation limits
def correct_saturation(data):
    max_val = np.max(data)
    signal_diff = np.diff(data)
    signal_diff_end = len(signal_diff)-1
    max_slope = 50000 #specify sat threshold
    
    for k in range(len(signal_diff)):
        if signal_diff[k] > max_slope:
            data[k+1:signal_diff_end] = data[k+1:signal_diff_end] - max_val
        if signal_diff[k] < -1*max_slope:
            data[k+1:signal_diff_end] = data[k+1:signal_diff_end] + max_val

# apply moving median filter, simplifies ppg filtering
def remove_outliers(data):
    data = np.array(medfilt(data.tolist(), 3))
    return data

# Calculate the heart rate and the beat variance 
def get_BPM_and_peak_variance(self, r_peaks):
    # Focus on standard HR frequency range
    segment = bandpass_filter(0.8, 2.5)

    peak_variance = np.finfo(np.float32).max
    if len(r_peaks) > 1:
        peak_variance = np.var(np.diff(r_peaks))

    time_difference = len(data)/250 #time range, given sampling rate 250hz
    time_difference_in_minutes = float(time_difference.seconds + float(time_difference.microseconds)/10**6)/60.0
    BPM = float(len(r_peaks)) / time_difference_in_minutes

    return BPM, peak_variance
from Prefiltering import*
from FiducialPoints import*
from Biomarkers2 import*
from Summary import*
from Statistics import*
from Biomarkers import *

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time
import wfdb
import pickle as pkl

import matplotlib.mlab


def create_annot(hr, fs):
    s = DotMap()
    s.v=hr
    s.fs=fs

    # print('creating prefiltering')
    s.filt_sig, s.filt_d1, s.filt_d2, s.filt_d3 = Prefiltering(s)
    fiducials = getFiducialsPoints(s)

    #######
    # print('creating biomarkers')
    ppg_biomarkers = Biomarkers2(s, fiducials)
    ######
    
    # # print('creating summary')
    # ppg_summary = Summary(s.v, fiducials['pk'], fiducials['os'], s.fs)
    # #######
    # # print('creating statistics')
    # ppg_statistics = Statistics(fiducials['pk'], fiducials['os'], ppg_biomarkers)

    # # create a dictionary of all the results and export to pkl file

    # results = {'hr': hr, 'fs': fs, 'fiducials': fiducials, 'biomarkers': ppg_biomarkers, 'summary': ppg_summary, 'statistics': ppg_statistics}
    
    #* results version which contains only the biomarkers
    results = {'hr': hr, 'fs': fs, 'fiducials': fiducials, 'biomarkers': ppg_biomarkers}
    return results


def create_annot_aux(num, folder_path, condition):
    print('creating AF annots')
    for n in range(1, num+1):
        print(n)
        try:
            str_num = '00' + str(n) if n < 10 else '0' + str(n)
            file_name = folder_path + f'mimic_perform_{condition}_{str_num}'
            signals, fields = wfdb.rdsamp(file_name)
            
            hr = signals[:, 0]
            mask = np.isnan(hr)
            hr[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), hr[~mask])
            hr = (np.rint(signals[:, 0] * 1000)).astype(int)
            fs  = fields['fs']
            results = create_annot(hr, fs)
            fiducials = results['fiducials']
            # print(len(fiducials['pk']))
            biomarkers = results['biomarkers']
            # print(len(biomarkers['BM_OSignal']))


            osignal_data = biomarkers['BM_OSignal'].ffill().bfill()
            osignal_data['sample_idx'] = (osignal_data['Tpi'] * fs).cumsum().astype(int)
            osignal_data['os'] = osignal_data['sample_idx'] + fiducials['os'][0] # add _set to the os column

            result = {'fs': fs, 'condition': {condition}, 'sample_num': str_num, 'fiducials': fiducials, 'biomarkers': biomarkers, 'osignal_data': osignal_data, 'hr': hr}
            
            # print(len(osignal_data), len(fiducials['os'])-1)
            if len(osignal_data) != len(fiducials['os'])-1:
                raise Exception('osignal_data length is not equal to fiducials os length')
        except:
            print(f'failed to process file {condition} {str_num}')
            continue
        
        with open(f'/home/meiri.yoav/biomed_proj/data/annotated/annot_{condition}_{str_num}.pkl', 'wb') as f:
            #! note that the osignal_data is not saved in the pkl file
            pkl.dump(result, f)

def create_all_annots():
    AF_num = 19
    non_AF_num = 16
    create_annot_aux(AF_num, '/home/meiri.yoav/biomed_proj/data/mimic_perform_af_wfdb/', 'af')
    create_annot_aux(non_AF_num, '/home/meiri.yoav/biomed_proj/data/mimic_perform_non_af_wfdb/', 'non_af')

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':
    create_all_annots()
    
    # # sig_path = 'D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/GODA_pyPPG/sample_data/PPG_sample_00.mat'
    # # sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl")])
    # sig_path = '/data/home/meiri.yoav/biomed_proj/data/mimic_perform_af_data.mat'

    # sig_format=sig_path[len(sig_path)-sig_path[::-1].index('.'):]
    # if sig_format=='mat':
    #     input_sig = scipy.io.loadmat(sig_path)
    #     # print(len(input_sig.get("data")[0]))
    #     print(input_sig.get("data")[0][1])
    #     hr = np.float64(input_sig.get("data")[0][0][0][0][0][0].squeeze())[0:]
    #     fs = input_sig.get("data")[0][0][0][0][0][1].squeeze()
    # elif sig_format=='csv':
    #     input_sig = np.loadtxt(sig_path, delimiter=',').astype(int)
    #     hr = input_sig
    #     fs = 125
    # elif sig_format == 'edf':
    #     input_sig = mne.io.read_raw_edf(sig_path)
    #     hr=-input_sig[22][0][0]
    #     fs = 125
    
    # go over patients 001 to 019 with AF:
    
    # failed to process file AF 004
    # failed to process file AF 005
    # failed to process file AF 006
    # failed to process file AF 007
    # failed to process file non_AF 004
    # failed to process file non_AF 006
    # failed to process file non_AF 008
    # failed to process file non_AF 012
    
    #-------------------------------------------------------------------------
    # folder_path = '/data/home/meiri.yoav/biomed_proj/data/mimic_perform_af_wfdb/' 
    # str_num = '005'
    # file_name = folder_path + f'mimic_perform_af_{str_num}'
    # signals, fields = wfdb.rdsamp(file_name)

    # hr = signals[:, 0]
    # mask = np.isnan(hr)
    # hr[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), hr[~mask])
    # hr = (np.rint(signals[:, 0] * 1000))

    # fs  = fields['fs']
    # results = create_annot(hr, fs)
    # fiducials = results['fiducials']
    # print(len(fiducials['pk']))
    # biomarkers = results['biomarkers']
    # print(len(biomarkers['BM_OSignal']))

    
    # osignal_data = biomarkers['BM_OSignal'].ffill().bfill()
    # osignal_data['sample_idx'] = (osignal_data['Tpi'] * fs).cumsum().astype(int)
    # osignal_data['os'] = osignal_data['sample_idx'] + fiducials['os'][0]
    # print(len(fiducials['os']-1), len(osignal_data['os']))
    # # for i in range(len(fiducials['os'])-1):
    # #     ith_in_data = osignal_data['os'].tolist()[i]
    # #     ith_in_original = fiducials['os'][i+1]
    # #     print(f'({ith_in_data}, {ith_in_original})')
        
     #-------------------------------------------------------------------------
    
    # plt.figure(figsize=(85,10))
    # plt.plot(hr,'k',linewidth=0.7)
    # plt.plot(fiducials['pk'].dropna().values, hr[fiducials['pk'].dropna().values.astype(int)], 'ro')
    # plt.plot(fiducials['os'].dropna().values, hr[fiducials['os'].dropna().values.astype(int)], 'bs')
    # plt.plot(fiducials['dn'].dropna().values, hr[fiducials['dn'].dropna().values.astype(int)], 'm*')

    # plt.xlim(right=fiducials['os'].dropna().values.astype(int)[37], left=0)  # adjust the right leaving left unchanged
    # plt.title(f'PPG signal with fiducials - {str_num}')
    # plt.legend(['signal', 'peak','onset','dic.notch'])
    # plt.xlabel('Sample (Fs='+str(fs)+' Hz)')
    # plt.ylabel('Amplitude')
    # plt.savefig("/data/home/meiri.yoav/biomed_proj/figure.png")

    # print('Program finished')

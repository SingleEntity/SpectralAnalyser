#make functional first
#check for number of columns, change behaviour depending

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import math

def read_spectral_data(spectral_data_path):
    gas_list = []
    df_final = []
    files = os.listdir(spectral_data_path)
    for file_index in range(0,len(files),1):

        #Reading each file containing OT for each IR gas
        gas = str(files[file_index])[0:len(files[file_index])-4]
        gas_list.append(gas)
        if file_index == 0: 
            df_final = pd.read_csv(spectral_data_path + str(files[file_index]))
            df_final.columns = ['wavenumber', gas]
        if file_index >0:
            df = pd.read_csv(spectral_data_path + str(files[file_index]))
            df.columns = ['wavenumber', gas]
            df[gas] = (1.0 - np.exp(-df[gas]))   #For fractional absorption
            df_final[gas] = df[gas]
    
    return (df_final, gas_list)

def read_instrument_bands(instrument_bands_path, spectral_res_nm):
    
    #Reading bands
    bands = pd.read_csv(instrument_bands_path)
    
    #Checking to see if there is a band width (2 values) or center point only
    if len(bands.columns) == 2:
                
        bands.columns = ['id', 'wl (nm)']
    
        #Converting to wavenumbers
        bands['wl (cm-1)'] = 1e7 / bands['wl (nm)']
    
        #Computing band windows
        bands['window start (cm-1)'] = bands['wl (cm-1)'] + spectral_res_nm
        bands['window end (cm-1)'] = bands['wl (cm-1)'] - spectral_res_nm
    elif len(bands.columns) == 3:
        bands.columns = ['id', 'window start (nm)', 'window end (nm)']
        
        #Converting to wavenumbers
        bands['wl (cm-1)'] = 1e7 / ((bands['window start (nm)'] + bands['window end (nm)']) / 2.0)
        bands['window start (cm-1)'] = 1e7 / bands['window start (nm)']
        bands['window end (cm-1)'] = 1e7 / bands['window end (nm)']      
    else:
        sys.exit("Error message")
  
    return(bands)

def compute_band_absorption(gas_list, spectral_data):
    
    #For each gas and each band we compute the spectral absorption (and mean absorption)
    for gas_index in range(len(gas_list)):
        absorbance_list = []
        nested_list = []
        for band_index in range(len(bands)):
            filter_start = spectral_data['wavenumber'] < (bands['window start (cm-1)'][band_index])
            filter_end = spectral_data['wavenumber'] > (bands['window end (cm-1)'][band_index])
            absorbance = np.mean(spectral_data[gas_list[gas_index]][filter_start & filter_end])
            if absorbance > 0.000001:
                absorbance_list.append(absorbance)
                nested_list.append(spectral_data[gas_list[gas_index]][filter_start & filter_end].tolist())
            else:
                absorbance_list.append(0.0)
                nested_list.append(spectral_data[gas_list[gas_index]][filter_start & filter_end].tolist())
        band_mean_absorption[gas_list[gas_index]] = absorbance_list
        band_spectral_absorption[gas_list[gas_index]] = nested_list
        
    #Writing absorption data frame to file
    band_mean_absorption.to_csv(output_path + 'Band-Species Absorbance Matrix.csv')

    return (band_mean_absorption, band_spectral_absorption)


if __name__ == '__main__':
    
    #Collecting command line arguments
    parser = argparse.ArgumentParser(description='Compute spectral response within instrument bands')

    parser.add_argument('--spectral_data_path', default='./IR and NIR/',
                    help='Folder containing spectral data from spectraplot.com')
    parser.add_argument('--instr_bands_file', default='WV3.csv',
                    help='File name for the instrument band csv')
    parser.add_argument('--spectral_res', default="2.5",
                    help='Spectral resolution in nm')
    parser.add_argument('--gas', default="CH4",
                    help='Gas species of interest e.g. CH4')
        
    args = parser.parse_args()
    
    spectral_data_path = args.spectral_data_path
    instrument_bands_path = args.instr_bands_file
    spectral_res_nm = float(args.spectral_res)
    gas_of_interest = args.gas

    output_path = "./output/"
 
    #Creating output folder if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    spectral_data_and_gas_list = read_spectral_data(spectral_data_path)
    spectral_data = spectral_data_and_gas_list[0]
    gas_list = spectral_data_and_gas_list[1]
     
    bands = read_instrument_bands(instrument_bands_path, spectral_res_nm)
    
    band_mean_absorption = pd.DataFrame(columns = ['band id'])
    band_mean_absorption['band id'] = bands['id']
    
    band_spectral_absorption = pd.DataFrame(columns = ['band id'])
    band_spectral_absorption['band id'] = bands['id']
    
    
    #Generating band to absorption data sets
    spectral_absorption = compute_band_absorption(gas_list, spectral_data)
    band_mean_absorption = spectral_absorption[0]
    band_spectral_absorption = spectral_absorption[1]    
    
    
    #Plotting each band ROI in frequency space and showing relavent gasious
    #species absorption spectra only
    for band_index in range(len(bands)):
        handles_list = []
        for gas_index in range(len(gas_list)):
            filter_start = spectral_data['wavenumber'] < (bands['window start (cm-1)'][band_index])
            filter_end = spectral_data['wavenumber'] > (bands['window end (cm-1)'][band_index])
            absorption = np.array(band_spectral_absorption[gas_list[gas_index]][band_index]) * 100
            if len(absorption) > 0:
                if gas_list[gas_index] == gas_of_interest:               
                    y_height = spectral_data[gas_of_interest][filter_start & filter_end].max() * 100
                    if not(math.isnan(y_height)): 
                        axes = plt.gca()
                        axes.set_ylim([0.0,y_height])
                if np.mean(absorption) > 0.00001:
                    plt.plot(spectral_data['wavenumber'][filter_start & filter_end], absorption)
                    handles_list.append(gas_list[gas_index])
        if len(handles_list) > 0:
            plt.legend(handles_list)
            plt.ylabel('I/Io %')
            plt.xlabel('wavenumber (cm-1)')
            plt.title('Band ' + str(bands['id'][band_index]))
            plt.axvline(x=bands['window start (cm-1)'][band_index],color='g',ls='--',linewidth=0.5)
            plt.axvline(x=bands['window end (cm-1)'][band_index],color='g',ls='--',linewidth=0.5)
            plt.savefig(output_path + 'Absorbance-Band-' + str(bands['id'][band_index]) + '.png')
            plt.clf()
    


    #Plotting just the gas of interest for each band
    if gas_of_interest != "": 
        for band_index in range(len(bands)):
            handles_list = []
            for gas_index in range(len(gas_list)):
                if gas_list[gas_index] == gas_of_interest:
                    filter_start = spectral_data['wavenumber'] < (bands['window start (cm-1)'][band_index])
                    filter_end = spectral_data['wavenumber'] > (bands['window end (cm-1)'][band_index])
                    absorption = np.array(band_spectral_absorption[gas_list[gas_index]][band_index]) * 100
                    if len(absorption) > 0:
                        y_height = spectral_data[gas_of_interest][filter_start & filter_end].max() * 100
                        if not(math.isnan(y_height)): 
                            axes = plt.gca()
                            axes.set_ylim([0.0,y_height])
                        if np.mean(absorption) > 0.00001:
                            plt.plot(spectral_data['wavenumber'][filter_start & filter_end], absorption)
                            handles_list.append(gas_list[gas_index])
            if len(handles_list) > 0:
                plt.legend(handles_list)
                plt.ylabel('I/Io %')
                plt.xlabel('wavenumber (cm-1)')
                plt.title('Band ' + str(bands['id'][band_index]))
                plt.axvline(x=bands['window start (cm-1)'][band_index],color='g',ls='--',linewidth=0.5)
                plt.axvline(x=bands['window end (cm-1)'][band_index],color='g',ls='--',linewidth=0.5)
                plt.savefig(output_path + 'Absorbance-Band-' + gas_of_interest + "-" + str(bands['id'][band_index]) + '.png')
                plt.clf()
    



    #Plotting band locations for each gas species
    first_band_filter = spectral_data['wavenumber'] < (bands['window start (cm-1)'][0] )
    last_band_filter = spectral_data['wavenumber'] > (bands['window end (cm-1)'][len(bands) - 1] )
    first_band = bands['window start (cm-1)'][0] 
    last_band = bands['window start (cm-1)'][len(bands) - 1] 
    for gas_index in range(len(gas_list)):
        absorption = np.array(spectral_data[gas_list[gas_index]][last_band_filter & first_band_filter]) * 100
        if len(absorption) > 0:
            if np.mean(absorption) > 0.00001:
                plt.plot(spectral_data['wavenumber'][last_band_filter & first_band_filter],absorption)
                plt.ylabel('I/Io %')
                plt.xlabel('wavenumber (cm-1)')
                plt.title(gas_list[gas_index])
                for band_index in range(len(bands)):
                    if bands['wl (cm-1)'][band_index] < np.max(spectral_data['wavenumber']) and \
                    bands['wl (cm-1)'][band_index] > np.min(spectral_data['wavenumber']):
                        plt.axvline(x=bands['wl (cm-1)'][band_index],color='g',ls='--',linewidth=0.5)
                plt.savefig(output_path + 'Band-Locations-' + gas_list[gas_index] + '.png')
                plt.clf()



    #Plotting band histogram absorption plot
    #Only showing top 5 absorbing species to avoid over clutter
    top_few = band_mean_absorption.sum().drop(labels=['band id'])
    top_few = top_few.sort_values(ascending=False)[0:8].keys()
     
    for gas_index in range(len(top_few)):
        if gas_index == 0: 
                plt.bar(band_mean_absorption['band id'].tolist(),band_mean_absorption[top_few[gas_index]])
                accumulation = band_mean_absorption[top_few[gas_index]]
        else: 
                plt.bar(band_mean_absorption['band id'].tolist(),band_mean_absorption[top_few[gas_index]], bottom=accumulation)
                accumulation = accumulation + band_mean_absorption[top_few[gas_index]]
    plt.ylabel('Absorption (cumulative) %')
    plt.xlabel('Band number')
    plt.legend(top_few)
    plt.xticks(rotation=45)
    plt.savefig(output_path + 'Absorbance histogram.png')
    plt.clf()

    
    #Plotting band relative/fraction histrogram
    band_mean_absorption_relative = band_mean_absorption.loc[:, band_mean_absorption.columns != 'band id']
    band_mean_absorption_relative = band_mean_absorption_relative.div(band_mean_absorption_relative.sum(axis=1), axis=0)

    for gas_index in range(len(top_few)):
        if gas_index == 0: 
                plt.bar(band_mean_absorption['band id'].tolist(),band_mean_absorption_relative[top_few[gas_index]])
                accumulation = band_mean_absorption_relative[top_few[gas_index]]
        else: 
                plt.bar(band_mean_absorption['band id'].tolist(),band_mean_absorption_relative[top_few[gas_index]], bottom=accumulation)
                accumulation = accumulation + band_mean_absorption_relative[top_few[gas_index]]
    plt.ylabel('Relative Absorption')
    plt.xlabel('Band number')
    plt.legend(top_few)
    plt.xticks(rotation=45)
    plt.savefig(output_path + 'Absorbance Fractional histogram.png')
    plt.clf()
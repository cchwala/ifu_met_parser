# ----------------------------------------------------------------------------
# Name:         
# Purpose:      
#
# Authors:      
#
# Created:      
# Copyright:    (c) Christian Chwala 2014
# Licence:      The MIT License
# ----------------------------------------------------------------------------

import csv
from datetime import datetime

import numpy as np
import xarray as xr
import tqdm

def read_parsivel_csv_files(fn,
                            spectrum_i_start,
                            spectrum_i_stop,
                            variables_index_dict,
                            show_warnings=True):
    data_dict = {'time': [],
                 'N': [],
                 'N_vec': []}
    for variable in variables_index_dict.keys():
        data_dict[variable] = []

    if not type(fn) == list:
        fn_list = [fn, ]
    else:
        fn_list = fn

    for fn_i in tqdm.tqdm(fn_list):
        with open(fn_i) as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            for row in reader:
                if len(row) == 0:
                    pass
                elif (len(row) == (spectrum_i_start + 1) or
                      len(row) == (spectrum_i_start + 1 + 32 * 32)):
                    date_str = row[0] + ' ' + row[1]
                    data_dict['time'].append(datetime.strptime(date_str, '%d.%m.%Y %H:%M:%S'))
                    for variable in variables_index_dict.keys():
                        data_dict[variable].append(float(row[variables_index_dict[variable]]))
                    if row[spectrum_i_start] == '<SPECTRUM>ZERO</SPECTRUM>':
                        N_matrix = np.zeros((32, 32), dtype=int)
                        N_array = np.zeros(32*32, dtype=int)
                    else:
                        # Convert the (mostly) empty strings of the spectrum
                        # to a list of integers
                        N_list = []
                        for value in row[spectrum_i_start:spectrum_i_stop]:
                            if value == '':
                                N_list.append(0)
                            elif value == '<SPECTRUM>':
                                N_list.append(0)
                            else:
                                N_list.append(int(value))
                        N_matrix = np.array(N_list).reshape((32, 32))

                        N_array = np.array(N_list)

                    data_dict['N'].append(N_matrix)
                    data_dict['N_vec'].append(N_array)
                else:
                    if show_warnings:
                        print('Warning: Skipping row. Length of row should be %d or %d, but is %d instead' % (
                              spectrum_i_start + 1, spectrum_i_start + 1 + 32 * 32, len(row)))
                        print(' Filename: %s' % fn_i)
                        print(' Row: %s' % row)
                        print('')

    ds = xr.Dataset(coords={'time': data_dict['time'],
                            'D': ('D', D),
                            'v': ('v', v),
                            'foo': ('foo', np.arange(32*32))},
                    data_vars={'N': (['time', 'v', 'D'], data_dict['N']),
                               'dD': ('D', dD),
                               'dv': ('v', dv),
                               'N_vec': (['time', 'foo'], data_dict['N_vec'])})
    for variable in variables_index_dict.keys():
        ds[variable] = ('time', data_dict[variable])

    return ds


def read_ifu_tereno_parsivel_csv_files(fn, show_warnings=True):
    variables_index_dict = {'Z': 7,
                            'R': 2}
    return read_parsivel_csv_files(fn,
                                   spectrum_i_start=16,
                                   spectrum_i_stop=-1,
                                   variables_index_dict=variables_index_dict,
                                   show_warnings=show_warnings)


def read_dlr_scalex_parsivel_csv_files(fn, show_warnings=True):
    variables_index_dict = {'Z': 4,
                            'R': 2}
    return read_parsivel_csv_files(fn,
                                   spectrum_i_start=9,
                                   spectrum_i_stop=-1,
                                   variables_index_dict=variables_index_dict,
                                   show_warnings=show_warnings)


# Classification according to volume-equivalent diameter
#
# Mid-value of class [mm]
D = np.array([0.062,  0.187,  0.312,  0.437,  0.562,
              0.687,  0.821,  0.937,  1.062,  1.187,
              1.375,  1.625,  1.875,  2.125,  2.375,
              2.750,  3.250,  3.750,  4.250,  4.750,
              5.500,  6.500,  7.500,  8.500,  9.500,
             11.000, 13.000, 15.000, 17.000, 19.000,
             21.500, 24.500])
# Class spread [mm]
dD = np.array([0.125, 0.125, 0.125, 0.125, 0.125,
               0.125, 0.125, 0.125, 0.125, 0.125,
               0.250, 0.250, 0.250, 0.250, 0.250,
               0.500, 0.500, 0.500, 0.500, 0.500,
               1.000, 1.000, 1.000, 1.000, 1.000,
               2.000, 2.000, 2.000, 2.000, 2.000,
               3.000, 3.000])
# Note: Diameter class 1 and 2 are not evaluated since they are out of the measurement
#       range of the OTT Parsivel2

# Classification according to particle speed
#
# Mid-value of class [m/s]
v = np.array([0.050,  0.150,  0.250,  0.350,  0.450,
              0.550,  0.650,  0.750,  0.850,  0.950,
              1.100,  1.300,  1.500,  1.700,  1.900,
              2.200,  2.600,  3.000,  3.400,  3.800,
              4.400,  5.200,  6.000,  6.800,  7.600,
              8.800, 10.400, 12.000, 13.600, 15.200,
             17.600, 20.800])
# Class spread in [m/s]
dv = np.array([0.100, 0.100, 0.100, 0.100, 0.100,
               0.100, 0.100, 0.100, 0.100, 0.100,
               0.200, 0.200, 0.200, 0.200, 0.200,
               0.400, 0.400, 0.400, 0.400, 0.400,
               0.800, 0.800, 0.800, 0.800, 0.800,
               1.600, 1.600, 1.600, 1.600, 1.600,
               3.200, 3.200])

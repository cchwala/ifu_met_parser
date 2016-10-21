#----------------------------------------------------------------------------
# Name:         
# Purpose:      
#
# Authors:      
#
# Created:      
# Copyright:    (c) Christian Chwala 2014
# Licence:      The MIT License
#----------------------------------------------------------------------------

import ftplib
import os
import tarfile
import tempfile
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
import wradlib as wrl


ftp_server = 'ftp-cdc.dwd.de'
data_dir_recent = '/pub/CDC/grids_germany/hourly/radolan/recent/bin'
data_dir_historic = '/pub/CDC/grids_germany/hourly/radolan/historical/bin'

# prefix and postfix RADOLAN files
filename_timestamp_format = {'radolan_recent': 'raa01-rw_10000-%y%m%d%H%M-dwd---bin.gz',
                             'radolan_historic': '%Y/RW%Y%m.tar.gz',
                             'netcdf': 'RADOLAN_%Y.nc'}


def download_files_from_ftp(t_start, t_stop, local_data_dir, redownload_existing_files=True):
    # Make sure that we work with `datetime` objects
    t_start = pd.to_datetime(t_start)
    t_stop = pd.to_datetime(t_stop)

    # We always set the minutes to 50, since the hourly RADOLAN
    # files are timestamped at 50 minutes.
    # Hence '2016-01' will be extended to '2016-01-01 00:50' which
    # avoids that historic data will be downloaded if only the recent
    # year is given as `t_start`, since '2016-01' would be before the
    # first "recent" file with the timestamp '2016-01-01 00:50'
    if t_start.minute != 50:
        t_start += timedelta(minutes=50)
    if t_stop.minute != 50:
        t_stop += timedelta(minutes=50)

    print('Trying to download RADOLAN files \n from %s \n till %s' % (t_start, t_stop))

    # init ftp session
    ftp_session = ftplib.FTP(ftp_server)
    ftp_session.login()

    # change to dir for recent files to get timestamp
    # of earliest files
    ftp_session.cwd(data_dir_recent)
    fn_list_recent = ftp_session.nlst()
    fn_list_recent.sort()

    # check if we need to look for `recent` and/or `historic` files
    t_recent_min = datetime.strptime(fn_list_recent[0],
                                     filename_timestamp_format['radolan_recent'])
    if t_start < t_recent_min:
        need_to_download_historic_files = True
    else:
        need_to_download_historic_files = False
    if t_stop >= t_recent_min:
        need_to_download_recent_files = True
    else:
        need_to_download_recent_files = False

    print 'We need to download \n recent files: ', need_to_download_recent_files
    print ' historic files:', need_to_download_historic_files

    # Generate the list of files that have to be downloaded
    fn_to_download_list = []
    for t in pd.date_range(t_start, t_stop, freq='H'):
        if t < t_recent_min:
            fn = t.strftime(filename_timestamp_format['radolan_historic'])
            fn_full_path = os.path.join(data_dir_historic, fn)
            if fn_full_path not in fn_to_download_list:
                fn_to_download_list.append(fn_full_path)
        elif t >= t_recent_min:
            fn = t.strftime(filename_timestamp_format['radolan_recent'])
            fn_full_path = os.path.join(data_dir_recent, fn)
            if fn_full_path not in fn_to_download_list:
                fn_to_download_list.append(fn_full_path)
        else:
            raise ValueError('It should be impossible to get here...')

    # Maker sure the local data dir exists
    if not os.path.exists(local_data_dir):
        os.makedirs(local_data_dir)

    # Get a list of files that are already in the local dir so that
    # they have not to be downloaded again
    fn_already_in_local_dir_list = os.listdir(local_data_dir)

    fn_downloaded_list = []
    for fn_full_path in fn_to_download_list:
        fn = os.path.split(fn_full_path)[1]
        fn_local_full_path = os.path.join(local_data_dir, fn)

        if fn not in fn_already_in_local_dir_list:
            print('Downloading %s' % fn_full_path)
            with open(fn_local_full_path, 'wb') as fh:
                ftp_session.retrbinary('RETR %s' % fn_full_path, fh.write)
        else:
            if redownload_existing_files:
                print('Redownloading %s' % fn_full_path)
                with open(fn_local_full_path, 'wb') as fh:
                    ftp_session.retrbinary('RETR %s' % fn_full_path, fh.write)
            else:
                print('Skipping %s' % fn_full_path)
        fn_downloaded_list.append(os.path.join(local_data_dir,fn))

    ftp_session.close()

    return fn_downloaded_list


def read_in_one_bin_file(f):
    data, metadata = wrl.io.read_RADOLAN_composite(f)
    return data, metadata


def read_in_files(fn_list, print_file_names=False):
    data_list = []
    metadata_list = []
    for fn in fn_list:
        if 'tar.gz' in fn:
            if print_file_names:
                print(' Untaring %s' % fn)
            # Unzip and untar RADOLAN-bin files and read them in
            with tarfile.open(fn, "r:gz") as tar:
                for member in tar.getmembers():
                    f = tar.extractfile(member)
                    if f is not None:
                        if print_file_names:
                            print('  Reading in %s' % f.name)
                        data, metadata = read_in_one_bin_file(f)
                        data = clean_radolan_data(data, metadata)
                        data_list.append(data)
                        metadata_list.append(metadata)
        else:
            if print_file_names:
                print(' Reading in %s' % fn)
            data, metadata = read_in_one_bin_file(fn)
            data = clean_radolan_data(data, metadata)
            data_list.append(data)
            metadata_list.append(metadata)

    ds = radolan_to_xarray_dataset(data_list, metadata_list)

    return ds


def clean_radolan_data(data, metadata):
    # mask invalid values
    sec = metadata['secondary']
    data.flat[sec] = -9999
    data[data == -9999] = np.nan
    return data


def radolan_to_xarray_dataset(data, metadata):
    if type(data) != list:
        data = [data,]
    if type(metadata) != list:
        metadata = [metadata,]

    radolan_lat_lon_grids = wrl.georef.get_radolan_grid(900,900, wgs84=True)
    radolan_lons = radolan_lat_lon_grids[:,:,0]
    radolan_lats = radolan_lat_lon_grids[:,:,1]

    ds = xr.Dataset({'precipitation': (['time', 'x', 'y'], data)},
                    coords={'lon': (['x', 'y'], radolan_lons),
                            'lat': (['x', 'y'], radolan_lats),
                            'time': [metadata_i['datetime'] for metadata_i in metadata],
                            'reference_time': pd.Timestamp('1970-01-01')})

    return ds




def _generate_temp_filenames_from_tar_gz_file(fn):
    # Unzip and untar RADOLAN-bin files and read them in
        tar = tarfile.open(fn, "r:gz")
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if f is not None:
                content = f.read()

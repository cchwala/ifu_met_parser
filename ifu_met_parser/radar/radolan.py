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
import gzip
import glob
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
import wradlib as wrl

from tqdm import tqdm

ftp_server = 'ftp-cdc.dwd.de'
data_dir_recent = '/pub/CDC/grids_germany/hourly/radolan/recent/bin'
data_dir_historic = '/pub/CDC/grids_germany/hourly/radolan/historical/bin'

# prefix and postfix RADOLAN files
filename_timestamp_format = {'radolan_recent': 'raa01-rw_10000-%y%m%d%H%M-dwd---bin.gz',
                             'radolan_historic': '%Y/RW%Y%m.tar.gz',
                             'netcdf': 'RADOLAN_%Y.nc'}


############################################
# Functions for downloading raw data files #
############################################


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
    ftp_session = ftplib.FTP(ftp_server, timeout=10)
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
            _download_one_file(ftp_session, fn_full_path, fn_local_full_path)
        else:
            if redownload_existing_files:
                print('Redownloading %s' % fn_full_path)
                _download_one_file(ftp_session, fn_full_path, fn_local_full_path)
            else:
                print('Skipping %s' % fn_full_path)
        fn_downloaded_list.append(os.path.join(local_data_dir,fn))

    ftp_session.close()

    return fn_downloaded_list


def _download_one_file(ftp_session, fn_remote, fn_local):
    try:
        with open(fn_local, 'wb') as fh:
            ftp_session.retrbinary('RETR %s' % fn_remote, fh.write)
    except:
        print(' Could not download %s' % fn_remote)
        print('  Second try...')
        try:
            with open(fn_local, 'wb') as fh:
                ftp_session.retrbinary('RETR %s' % fn_remote, fh.write)
        except:
            print('  Could not download %s' % fn_remote)
            os.remove(fn_local)


########################################
# Functions for parsing raw data files #
########################################

def read_in_one_bin_file(f):
    data, metadata = wrl.io.read_RADOLAN_composite(f)
    return data, metadata


def read_in_one_tar_gz_file(fn, print_filenames=False):
    data_list = []
    metadata_list = []
    # Unzip and untar RADOLAN-bin files and read them in
    with tarfile.open(fn, "r:gz") as tar:
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if f is not None:
                if print_filenames:
                    print('  Reading in %s' % f.name)
                # Unfortunately, the old "historic" files till somewhere
                # in 2014 contain gzipped "bin" files. And since the
                # extraction from the tar-archive returns file handles,
                # assuming that they belong to plain files, the handles
                # have to be converted into GzipFile-handles.
                # Note, that the wradlib parsing function could handle zipped
                # and plain files, but only when file names are provided
                if '.gz' in f.name:
                    with gzip.GzipFile(fileobj=f, mode='rb') as f:
                        data, metadata = read_in_one_bin_file(f)
                else:
                    data, metadata = read_in_one_bin_file(f)
                data = _clean_radolan_data(data, metadata)
                data_list.append(data)
                metadata_list.append(metadata)

    # Sort by datetime
    data_list_sorted, metadata_list_sorted = _sort_by_datetime(data_list,
                                                               metadata_list)
    return data_list_sorted, metadata_list_sorted


def read_in_files(fn_list, print_filenames=False):
    data_list = []
    metadata_list = []
    for fn in fn_list:
        # Distinguish between the two different file types,
        # the "historic" files, which are a tar.gz-archive
        # of the "bin" files and the "recent" files which are
        # also "bin" files.
        if 'tar.gz' in fn:
            if print_filenames:
                print(' Untaring %s' % fn)
            (temp_data_list,
             temp_metadata_list) = read_in_one_tar_gz_file(
                                      fn,
                                      print_filenames=print_filenames)
            data_list += temp_data_list
            metadata_list += temp_metadata_list
        else:
            if print_filenames:
                print(' Reading in %s' % fn)
            data, metadata = read_in_one_bin_file(fn)
            data = _clean_radolan_data(data, metadata)
            data_list.append(data)
            metadata_list.append(metadata)

    # ds = radolan_to_xarray_dataset(data_list, metadata_list)

    return data_list, metadata_list


def _clean_radolan_data(data, metadata):
    # mask invalid values
    sec = metadata['secondary']
    data.flat[sec] = -9999
    data[data == -9999] = np.nan
    return data


def _sort_by_datetime(data_list, metadata_list):
    datetime_array = np.array([metadata_i['datetime'] for metadata_i in metadata_list])
    sorting_index = datetime_array.argsort()

    data_list_sorted = [data_list[i] for i in sorting_index]
    metadata_list_sorted = [metadata_list[i] for i in sorting_index]

    return data_list_sorted, metadata_list_sorted


###################################
# Functions for writing to NetCDF #
###################################

def radolan_to_xarray_dataset(data, metadata):
    if type(data) != list:
        data = [data, ]
    if type(metadata) != list:
        metadata = [metadata, ]

    # Get coordinates
    radolan_xy_grids = wrl.georef.get_radolan_grid(900,900)
    radolan_x = radolan_xy_grids[0, :, 0]
    radolan_y = radolan_xy_grids[:, 0, 1]
    radolan_lat_lon_grids = wrl.georef.get_radolan_grid(900,900, wgs84=True)
    radolan_lons = radolan_lat_lon_grids[:, :, 0]
    radolan_lats = radolan_lat_lon_grids[:, :, 1]

    # Sort data by datetime
    data, metadata = _sort_by_datetime(data, metadata)
    datetime_list_sorted = [metadata_i['datetime'] for metadata_i in metadata]

    # Create a data set
    ds = xr.Dataset({'rainfall_amount': (['time', 'y', 'x'], data)},
                    coords={'x': radolan_x,
                            'y': radolan_y,
                            'lon': (['y', 'x'], radolan_lons),
                            'lat': (['y', 'x'], radolan_lats),
                            'time': datetime_list_sorted,
                            'reference_time': pd.Timestamp('1970-01-01')})

    # Add metadata to data set

    return ds


def append_to_yearly_netcdf(netcdf_file_dir, ds, overwrite_all=False):
    netcdf_fn_list = []
    for year, ds_yearly in ds.groupby('time.year'):
        fn = datetime.strftime(pd.to_datetime(ds_yearly.time.values[0]),
                               format=filename_timestamp_format['netcdf'])
        fn_full_path = os.path.join(netcdf_file_dir, fn)
        netcdf_fn_list.append(fn_full_path)
        if not os.path.isfile(fn_full_path) or overwrite_all:
            netcdf_mode = 'w'
            print('Writing to new file %s' % fn_full_path)
            if not os.path.isdir(os.path.split(fn_full_path)[0]):
                os.makedirs(os.path.split(fn_full_path)[0])
        else:
            netcdf_mode = 'a'
            print('Appending to %s' % fn_full_path)
        ds_yearly.to_netcdf(fn_full_path,
                            mode=netcdf_mode,
                            encoding={'rainfall_amount': {'dtype': 'float32',
                                                          'zlib': True,
                                                          'complevel': 1}})
    return netcdf_fn_list


def create_yearly_netcdfs(raw_data_dir,
                          netcdf_file_dir,
                          clear_netcdf_files=True,
                          print_filenames=False):

    # delete already existing netcdf files if desired
    if clear_netcdf_files:
        for fn in glob.glob(os.path.join(netcdf_file_dir, '*.nc')):
            print('Removing %s' % fn)
            os.remove(fn)

    # Get the list of raw data files
    fn_list_historic = glob.glob(os.path.join(raw_data_dir, 'RW*.tar.gz'))
    fn_list_recent = glob.glob(os.path.join(raw_data_dir, 'raa*bin.gz'))

    # Write NetCDF files for historic files. The parsing and appending
    # is done for each file separately to limit the size of data that
    # has to be kept in memory
    print('Reading and parsing %d historic files' % len(fn_list_historic))
    for fn in tqdm(fn_list_historic):
        data_list_temp, metadata_list_temp = read_in_files(
                                                [fn, ],
                                                print_filenames=print_filenames)
        ds = radolan_to_xarray_dataset(data_list_temp, metadata_list_temp)
        append_to_yearly_netcdf(netcdf_file_dir, ds, overwrite_all=False)

    # Write NetCDF files for recent data, but read in the files
    # chunk-wise to avoid storing everything in memory
    chunk_length = 500
    fn_list_recent_chunks = [fn_list_recent[i:i + chunk_length]
                             for i in xrange(0, len(fn_list_recent), chunk_length)]

    print('Reading and parsing %d recent files in %d chunks' % (len(fn_list_recent),
                                                                len(fn_list_recent_chunks)))
    for fn_list_chunk in tqdm(fn_list_recent_chunks):
        data_list_temp, metadata_list_temp = read_in_files(fn_list_chunk,
                                                           print_filenames=print_filenames)
        ds = radolan_to_xarray_dataset(data_list_temp, metadata_list_temp)
        append_to_yearly_netcdf(netcdf_file_dir, ds, overwrite_all=False)





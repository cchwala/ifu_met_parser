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
from time import sleep

import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
import wradlib as wrl

from tqdm import tqdm

ftp_server = 'ftp-cdc.dwd.de'
data_dir_recent = '/pub/CDC/grids_germany/hourly/radolan/recent/bin'
data_dir_historic = '/pub/CDC/grids_germany/hourly/radolan/historical/bin'

# prefix and postfix RADOLAN files
filename_timestamp_format = {'radolan_recent': 'raa01-rw_10000-%y%m%d%H%M-dwd---bin.gz',
                             'radolan_historic': '%Y/RW%Y%m.tar.gz',
                             'netcdf': 'RADOLAN_RW_%Y.nc'}


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
    try:
        ftp_session = ftplib.FTP(ftp_server, timeout=10)
    except:
        raise IOError('Could not connect to FTP server')

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
            fn_downloaded = _download_one_file(ftp_session, fn_full_path, fn_local_full_path)
        else:
            if redownload_existing_files:
                print('Redownloading %s' % fn_full_path)
                fn_downloaded = _download_one_file(ftp_session, fn_full_path, fn_local_full_path)
            else:
                print('Skipping %s' % fn_full_path)
                fn_downloaded = fn_local_full_path
        if fn_downloaded is not None:
            fn_downloaded_list.append(fn_downloaded)

    ftp_session.close()

    return fn_downloaded_list


def _download_one_file(ftp_session, fn_remote, fn_local):
    try:
        with open(fn_local, 'wb') as fh:
            ftp_session.retrbinary('RETR %s' % fn_remote, fh.write)
            return fn_local
    except:
        print(' Could not download %s' % fn_remote)
        print('  Second try...')
        try:
            with open(fn_local, 'wb') as fh:
                ftp_session.retrbinary('RETR %s' % fn_remote, fh.write)
                return fn_local
        except:
            print('  Could not download %s' % fn_remote)
            os.remove(fn_local)
            return None


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
            try:
                data, metadata = read_in_one_bin_file(fn)
                data = _clean_radolan_data(data, metadata)
                data_list.append(data)
                metadata_list.append(metadata)
            except:
                print('  !!!!!!!! Could not read in file %s !!!!!!!!!!!' % fn)

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

def create_empty_netcdf(fn):
    with netCDF4.Dataset(fn, 'w') as nc_fh:
        # Get RADOLAN coordinates
        radolan_xy_grids = wrl.georef.get_radolan_grid(900, 900)
        radolan_x = radolan_xy_grids[0, :, 0]
        radolan_y = radolan_xy_grids[:, 0, 1]
        radolan_lat_lon_grids = wrl.georef.get_radolan_grid(900, 900, wgs84=True)
        radolan_lons = radolan_lat_lon_grids[:, :, 0]
        radolan_lats = radolan_lat_lon_grids[:, :, 1]

        # create dimensions
        nc_fh.createDimension('x', 900)
        nc_fh.createDimension('y', 900)
        nc_fh.createDimension('time', None)

        # create variables
        nc_fh.createVariable('x', 'f8', ('x'))
        nc_fh.createVariable('y', 'f8', ('y'))
        nc_fh.createVariable('latitudes', 'f8', ('y', 'x'))
        nc_fh.createVariable('longitudes', 'f8', ('y', 'x'))
        nc_fh.createVariable('time', 'f8', ('time'))
        nc_fh.createVariable('rainfall_amount', 'f8', ('time', 'y', 'x'),
                             zlib=True, complevel=5, fill_value=-9999.9)

        # variable attributes
        nc_fh['time'].long_name = 'Time'
        nc_fh['time'].standard_name = 'time'
        nc_fh['time'].units = 'hours since 2000-01-01 00:50:00.0'
        nc_fh['time'].calendar = 'standard'

        nc_fh['x'].long_name = 'RADOLAN Grid x coordinate of projection'
        nc_fh['x'].standard_name = 'projection_x_coordinate'
        nc_fh['x'].units = 'km'

        nc_fh['y'].long_name = 'RADOLAN Grid y coordinate of projection'
        nc_fh['y'].standard_name = 'projection_y_coordinate'
        nc_fh['y'].units = 'km'

        nc_fh['latitudes'].long_name = 'Latitude'
        nc_fh['latitudes'].standard_name = 'latitude'
        nc_fh['latitudes'].units = 'degrees_north'

        nc_fh['longitudes'].long_name = 'Longitude'
        nc_fh['longitudes'].standard_name = 'longitude'
        nc_fh['longitudes'].units = 'degrees_east'

        nc_fh['rainfall_amount'].long_name = 'Hourly rainfall'
        nc_fh['rainfall_amount'].standard_name = 'rainfall_amount'
        nc_fh['rainfall_amount'].units = 'kg s-1'
        nc_fh['rainfall_amount'].coordinates = 'latitudes longitudes'
        nc_fh['rainfall_amount'].grid_mapping = 'RADOLAN_grid'

        # global attributes
        nc_fh.title = 'RADOLAN RW rainfall data'
        nc_fh.source = 'ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/hourly/radolan/'
        nc_fh.institution = 'DWD'
        nc_fh.history = 'Created at ' + str(datetime.utcnow())
        nc_fh.Conventions = 'CF 1.6'

        # Add actual coordinate data
        nc_fh['latitudes'][:, :] = radolan_lats
        nc_fh['longitudes'][:, :] = radolan_lons
        nc_fh['x'][:] = radolan_x
        nc_fh['y'][:] = radolan_y

        # Add projection definition
        nc_fh.createVariable('radolan_grid', 'f8')
        nc_fh['radolan_grid'].long_name = 'RADOLAN Grid'
        nc_fh['radolan_grid'].grid_mapping_name = 'polar_stereographic'
        nc_fh['radolan_grid'].semi_major_axis = 6370040.0
        nc_fh['radolan_grid'].false_easting = 0.0
        nc_fh['radolan_grid'].false_northing = 0.0
        nc_fh['radolan_grid'].standard_parallel = 10
        nc_fh['radolan_grid'].straight_vertical_longitude_from_pole = 30
        nc_fh['radolan_grid'].latitude_of_projection_origin = 90.0

    return None


def append_to_netcdf(fn, data_list, metadata_list):
    if type(data_list) != list:
        data_list = [data_list, ]
    if type(metadata_list) != list:
        metadata_list = [metadata_list, ]
    with netCDF4.Dataset(fn, 'a') as nc_fh:
        current_length = len(nc_fh['time'][:])
        for i, (data, metadata) in enumerate(zip(data_list, metadata_list)):
            i_new = i + current_length
            nc_fh['time'][i_new] = netCDF4.date2num(metadata['datetime'],
                                                    units=nc_fh['time'].units,
                                                    calendar=nc_fh['time'].calendar)
            nc_fh['rainfall_amount'][i_new, :, :] = data


def append_to_yearly_netcdf(netcdf_file_dir,
                            data_list,
                            metadata_list):
    for data, metadata in zip(data_list, metadata_list):
        netcdf_fn = datetime.strftime(metadata['datetime'],
                                      filename_timestamp_format['netcdf'])
        netcdf_fn_full_path = os.path.join(netcdf_file_dir, netcdf_fn)
        if not os.path.isfile(netcdf_fn_full_path):
            print('Creating new file %s' % netcdf_fn_full_path)
            create_empty_netcdf(netcdf_fn_full_path)
        # print('Appending to %s for t=%s' % (netcdf_fn_full_path,
        #                                    str(metadata['datetime'])))
        append_to_netcdf(netcdf_fn_full_path, data, metadata)

    return


def create_yearly_netcdfs(raw_data_dir,
                          netcdf_file_dir,
                          clear_netcdf_files=True,
                          print_filenames=False):

    if not os.path.isdir(netcdf_file_dir):
        os.mkdir(netcdf_file_dir)

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
    # has to be kept in memory. Please note that one historic tar.gz-file
    # still already contains more than 700 bin files and hence the parsed
    # data will be roughly 5 GB.
    print('Reading and parsing %d historic files' % len(fn_list_historic))
    for fn in fn_list_historic:
        print('Working on %s...' % fn)
        data_list_temp, metadata_list_temp = read_in_files(
                                                [fn, ],
                                                print_filenames=print_filenames)
        append_to_yearly_netcdf(netcdf_file_dir,
                                data_list_temp,
                                metadata_list_temp)

    # Write NetCDF files for recent data, but read in the files
    # chunk-wise to avoid storing everything in memory
    chunk_length = 500
    fn_list_recent_chunks = [fn_list_recent[i:i + chunk_length]
                             for i in xrange(0, len(fn_list_recent), chunk_length)]

    print('Reading and parsing %d recent files in %d chunks' %
          (len(fn_list_recent),
           len(fn_list_recent_chunks)))
    for i, fn_list_chunk in enumerate(fn_list_recent_chunks):
        print('Working on chunk %d/%d with %d files...' %
              (i+1,
               len(fn_list_recent_chunks),
               len(fn_list_chunk)))
        data_list_temp, metadata_list_temp = read_in_files(fn_list_chunk,
                                                           print_filenames=print_filenames)
        append_to_yearly_netcdf(netcdf_file_dir,
                                data_list_temp,
                                metadata_list_temp)


################################################################
# Functions for continuously updating NetCDFs with recent data #
################################################################

def get_last_time_stamp_from_netcdfs(netcdf_file_dir, fn_pattern='RADOLAN_*.nc'):
    fn_list = glob.glob(os.path.join(netcdf_file_dir, fn_pattern))
    ds = xr.open_mfdataset(fn_list)
    ts_max = ds.time[:].values.max()
    ds.close()
    return ts_max


def download_latest_files_from_ftp(local_data_dir, netcdf_file_dir):
    netcdf_ts_last = get_last_time_stamp_from_netcdfs(netcdf_file_dir)
    fn_list = download_files_from_ftp(t_start=netcdf_ts_last + np.timedelta64(1, 'h'),
                                      t_stop=datetime.utcnow(),
                                      local_data_dir=local_data_dir,
                                      redownload_existing_files=False)
    return fn_list


def update_recent_netcdf(local_data_dir, netcdf_file_dir, N_retries=10, wait_sec=10):
    retries = 0
    while retries < N_retries:
        try:
            if retries > 0:
                print('Try #%d to download file via ftp' % retries)
            fn_list = download_latest_files_from_ftp(local_data_dir,
                                                     netcdf_file_dir)
            break
        except IOError, e:
            print(e)
            retries += 1
            sleep(wait_sec)

    if retries == N_retries:
        raise IOError('Could not download file, even after %d retries' % retries)


    # Write NetCDF files, but read in the files
    # chunk-wise to avoid storing everything in memory
    chunk_length = 200
    fn_list_chunks = [fn_list[i:i + chunk_length]
                             for i in xrange(0, len(fn_list), chunk_length)]

    print('Reading and parsing %d recent files in %d chunks' %
          (len(fn_list),
           len(fn_list_chunks)))
    for i, fn_list_chunk in enumerate(fn_list_chunks):
        print('Working on chunk %d/%d with %d files...' %
              (i+1,
               len(fn_list_chunks),
               len(fn_list_chunk)))
        data_list_temp, metadata_list_temp = read_in_files(fn_list_chunk,
                                                           print_filenames=True)
        append_to_yearly_netcdf(netcdf_file_dir,
                                data_list_temp,
                                metadata_list_temp)

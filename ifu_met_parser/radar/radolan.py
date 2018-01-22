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


ftp_server = 'ftp-cdc.dwd.de'
data_dir_recent = '/pub/CDC/grids_germany/hourly/radolan/recent/bin'
data_dir_historic = '/pub/CDC/grids_germany/hourly/radolan/historical/bin'

# prefix and postfix RADOLAN files
filename_timestamp_format = {'radolan_recent': 'raa01-rw_10000-%y%m%d%H%M-dwd---bin.gz',
                             'radolan_historic': '%Y/RW%Y%m.tar.gz',
                             'netcdf': 'RADOLAN_RW_%Y.nc'}


# Constants for the RADOLAN rainfall_amount NetCDF-variable
NETCDF_SCALING_FACTOR = 0.1
NETCDF_DATA_TYPE = 'i2'
# Value that represents missing data in the NetCDF file
# for the RADOLAN rainfall_amount
NETCDF_FILL_VALUE_INT = -9999
# Value that has to be used in the data array that will be written
# to the NetCDF file. This value has to be scaled, so that it
# matches what ends up in the file
NETCDF_FILL_VALUE_FLOAT = float(NETCDF_FILL_VALUE_INT)*NETCDF_SCALING_FACTOR


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

    fn_list_recent_cleaned = []
    for fn in fn_list_recent:
        if 'dwd---bin.gz' in fn:
            fn_list_recent_cleaned.append(fn)

    # check if we need to look for `recent` and/or `historic` files
    t_recent_min = datetime.strptime(fn_list_recent_cleaned[0],
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
                try:
                    if '.gz' in f.name:
                        with gzip.GzipFile(fileobj=f, mode='rb') as f:
                            data, metadata = read_in_one_bin_file(f)
                    else:
                        data, metadata = read_in_one_bin_file(f)
                    data = _clean_radolan_data(data, metadata)
                    data_list.append(data)
                    metadata_list.append(metadata)
                except:
                    print('  !!!!!!!! Could not read in file %s !!!!!!!!!!!' % fn)

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
        nc_fh.createVariable('rainfall_amount', NETCDF_DATA_TYPE, ('time', 'y', 'x'),
                             zlib=True, complevel=5,
                             fill_value=NETCDF_FILL_VALUE_INT)
        nc_fh['rainfall_amount'].scale_factor = NETCDF_SCALING_FACTOR
        nc_fh['rainfall_amount'].add_offset = 0

        nc_fh.set_auto_maskandscale(True)

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
        nc_fh['rainfall_amount'].units = 'kg m-2'
        nc_fh['rainfall_amount'].coordinates = 'longitudes latitudes'
        nc_fh['rainfall_amount'].grid_mapping = 'RADOLAN_grid'

        # global attributes
        nc_fh.title = 'RADOLAN RW rainfall data'
        nc_fh.source = 'ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/hourly/radolan/'
        nc_fh.institution = 'DWD'
        nc_fh.history = 'Created at ' + str(datetime.utcnow())
        nc_fh.Conventions = 'CF-1.6'

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
        nc_fh['radolan_grid'].scale_factor_at_projection_origin = 0.9330127019
        nc_fh['radolan_grid'].straight_vertical_longitude_from_pole = 10.0
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
            t_new = netCDF4.date2num(metadata['datetime'],
                                     units=nc_fh['time'].units,
                                     calendar=nc_fh['time'].calendar)
            if t_new not in nc_fh['time'][:]:
                nc_fh['time'][i_new] = t_new
                temp_data = data.copy()
                temp_data[np.isnan(temp_data)] = NETCDF_FILL_VALUE_FLOAT
                nc_fh['rainfall_amount'][i_new, :, :] = temp_data
            else:
                print('Data for %s is already in %s' % (t_new, fn))


def append_to_yearly_netcdf(netcdf_file_dir,
                            data_list,
                            metadata_list):
    for i, (data, metadata) in enumerate(zip(data_list, metadata_list)):
        netcdf_fn = datetime.strftime(metadata['datetime'],
                                      filename_timestamp_format['netcdf'])
        netcdf_fn_full_path = os.path.join(netcdf_file_dir, netcdf_fn)
        if not os.path.isfile(netcdf_fn_full_path):
            print('Creating new file %s' % netcdf_fn_full_path)
            create_empty_netcdf(netcdf_fn_full_path)
        else:
            if i == 0:
                print('Appending to %s' % netcdf_fn_full_path)
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

    fn_list_historic.sort()
    fn_list_recent.sort()

    # Write NetCDF files for historic files.
    read_in_files_and_append_to_yearly_netcdf(fn_list_historic,
                                              netcdf_file_dir,
                                              print_filenames)

    # Write NetCDF files for recent files. The `fn_list` gets
    # divided into chunks within the function to avoid that too
    # the data that has to be temporarily stored in memory gets
    # to large.
    chunk_length = 100
    read_in_files_chunked_and_append_to_yearly_netcdf(fn_list_recent,
                                                      chunk_length,
                                                      netcdf_file_dir,
                                                      print_filenames)


def read_in_files_and_append_to_yearly_netcdf(fn_list,
                                              netcdf_file_dir,
                                              print_filenames=False):
    if type(fn_list) != list:
        fn_list = [fn_list, ]

    data_list = []
    metadata_list = []
    print('Reading and parsing %d files...' % len(fn_list))
    for fn in fn_list:
        data_list_temp, metadata_list_temp = read_in_files(
                                                [fn, ],
                                                print_filenames=print_filenames)
        #data_list.extend(data_list_temp)
        #metadata_list.extend(metadata_list_temp)

        append_to_yearly_netcdf(netcdf_file_dir,
                                data_list_temp,
                                metadata_list_temp)

        data_list_temp = None
        metadata_list_temp = None


def read_in_files_chunked_and_append_to_yearly_netcdf(fn_list,
                                                      chunk_length,
                                                      netcdf_file_dir,
                                                      print_filenames=False):
    if type(fn_list) != list:
        fn_list = [fn_list, ]

    fn_list_chunks = [fn_list[i:i + chunk_length]
                      for i in xrange(0, len(fn_list), chunk_length)]
    print('Reading and parsing %d files in %d chunks' %
          (len(fn_list),
           len(fn_list_chunks)))

    for i, fn_list_chunk in enumerate(fn_list_chunks):
        print(' Working on chunk %d/%d with %d files...' %
              (i+1,
               len(fn_list_chunks),
               len(fn_list_chunk)))
        read_in_files_and_append_to_yearly_netcdf(fn_list_chunk,
                                                  netcdf_file_dir,
                                                  print_filenames=print_filenames)


################################################################
# Functions for continuously updating NetCDFs with recent data #
################################################################

def get_last_time_stamp_from_netcdfs(netcdf_file_dir, fn_pattern='RADOLAN_*.nc'):
    fn_list = glob.glob(os.path.join(netcdf_file_dir, fn_pattern))
    if len(fn_list) > 0:
        ds = xr.open_mfdataset(fn_list)
        ts_max = ds.time[:].values.max()
        ds.close()
    else:
        # If no NetCDFs exist, start with first RADOLAN data from summer 2005
        ts_max = pd.np.datetime64('2005-06-01', 'ns')

    return ts_max


def download_latest_files_from_ftp(local_data_dir, netcdf_file_dir):
    netcdf_ts_last = get_last_time_stamp_from_netcdfs(netcdf_file_dir)
    fn_list = download_files_from_ftp(t_start=netcdf_ts_last + np.timedelta64(1, 'h'),
                                      t_stop=datetime.utcnow(),
                                      local_data_dir=local_data_dir,
                                      redownload_existing_files=False)
    return fn_list


def update_recent_netcdf(local_data_dir, netcdf_file_dir, N_retries=10, wait_sec=10):

    print_filenames = True

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

    # Divide the list into historic tar.gz-files and current .gz-files
    fn_list_historic = []
    fn_list_recent = []
    for fn in fn_list:
        if 'bin.gz' in fn:
            fn_list_recent.append(fn)
        elif 'tar.gz' in fn:
            fn_list_historic.append(fn)
        else:
            print(' !!! File name not recognized !!!')

    fn_list_historic.sort()
    fn_list_recent.sort()

    # Write NetCDF files for historic files.
    read_in_files_and_append_to_yearly_netcdf(fn_list_historic,
                                              netcdf_file_dir,
                                              print_filenames)

    # Write NetCDF files for recent files. The `fn_list` gets
    # divided into chunks within the function to avoid that too
    # the data that has to be temporarily stored in memory gets
    # to large.
    chunk_length = 100
    read_in_files_chunked_and_append_to_yearly_netcdf(fn_list_recent,
                                                      chunk_length,
                                                      netcdf_file_dir,
                                                      print_filenames)


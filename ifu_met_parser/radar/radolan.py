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

import pandas as pd
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


def read_in_one_bin_file(fn):
    data, metadata = wrl.io.read_RADOLAN_composite(fn)
    return data, metadata

def read_in_one_bin_file_from_file_handle(fh):
    with tempfile.NamedTemporaryFile() as tmpfile:
        temp_file_name = tmpfile.name
        f

def read_in_files(fn_list):
    data_list = []
    metadata_list = []
    for fn in fn_list:
        if 'tar.gz' in fn:
            # Unzip and untar RADOLAN-bin files and read them in
            tar = tarfile.open(fn, "r:gz")
            for member in tar.getmembers():
                f = tar.extractfile(member)
                if f is not None:
                    content = f.read()


def _generate_temp_filenames_from_tar_gz_file(fn):
    # Unzip and untar RADOLAN-bin files and read them in
        tar = tarfile.open(fn, "r:gz")
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if f is not None:
                content = f.read()

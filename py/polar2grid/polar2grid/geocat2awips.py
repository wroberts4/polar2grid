#!/usr/bin/env python
# encoding: utf-8
"""Script that uses the `polar2grid` toolbox of modules to take Geocat
hdf4 (.hdf) files and create a properly scaled AWIPS compatible NetCDF file.

:author:       Eva Schiffer (evas)
:contact:      eva.schiffer@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2013 University of Wisconsin SSEC. All rights reserved.
:date:         Oct 2012
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

from polar2grid.core            import Workspace
from polar2grid.core.constants  import *
from polar2grid.core.glue_utils import *
from polar2grid.core.time_utils import utc_now
from .grids.grids               import create_grid_jobs, Cartographer
from polar2grid.modis           import FILE_CONTENTS_GUIDE, find_nav_uid_from_filename
from polar2grid.modis           import Geo_Frontend
from .awips                     import Backend
import remap

import os
import sys
import re
import logging

from multiprocessing import Process
from datetime        import datetime

log       = logging.getLogger(__name__)
GLUE_NAME = "geocat2awips"
LOG_FN    = os.environ.get("GEOCAT2AWIPS_LOG", None) # None interpreted in main

# TODO this is currently being used in a slipshod manner, move to the newer p2g pattern of deleting files
def clean_up_files():
    """Remove as many of the possible files that were created from a previous
    run of this script, including temporary files.

    :note:
        This does not remove the log file because it requires the log file
        to report what's being removed.
    """
    
    list_to_remove = ["latitude*.real4.*",
                      "longitude*.real4.*",
                      "image*.real4.*", 
                      "ll2cr_*.img",
                      "result*.real4.*",
                      "SSEC_AWIPS_GEOCAT*" ]
    
    remove_file_patterns(list_to_remove)

def process_data_sets(filepaths,
                      nav_uid,
                      fornav_D=None, fornav_d=None,
                      grid_configs=None,
                      fornav_m=True,     #
                      forced_grid=None,
                      forced_gpd=None, forced_nc=None,
                      num_procs=1,
                      rescale_config=None,
                      backend_config=None,
                      ) :
    """Process all the files provided from start to finish,
    from filename to AWIPS NC file.
    
    Note: all files provided are expected to share a navigation source.
    """
    
    log.debug("Processing %s navigation set" % (nav_uid,))
    status_to_return = STATUS_SUCCESS
    
    # Handle parameters
    grid_configs = grid_configs or tuple() # needs to be a tuple for use
    
    # create the front and backend objects
    # these calls should load any configuration files needed
    cartographer    = Cartographer(*grid_configs)
    frontend_object = Geo_Frontend()
    backend_object  = Backend(rescale_config=rescale_config, backend_config=backend_config)
    
    # Extract Swaths
    log.info("Extracting swaths...")
    meta_data = { }
    try:
        meta_data = frontend_object.make_swaths(nav_uid, filepaths, cut_bad=True)
    except StandardError:
        log.error("Swath creation failed")
        log.debug("Swath creation error:", exc_info=1)
        status_to_return |= STATUS_FRONTEND_FAIL
    
    # if we weren't able to load any of the swaths... stop now
    if len(meta_data.keys()) <= 0 :
        log.error("Unable to load swaths for any of the bands, quitting...")
        return status_to_return or STATUS_UNKNOWN_FAIL
    
    # for convenience, pull some things out of the meta data
    sat = meta_data["sat"]
    instrument = meta_data["instrument"]
    start_time = meta_data["start_time"]
    band_info = meta_data["bands"]
    flatbinaryfilename_lat = meta_data["fbf_lat"]
    flatbinaryfilename_lon = meta_data["fbf_lon"]
    lat_fill = meta_data.get("lat_fill_value", None)
    lon_fill = meta_data.get("lon_fill_value", None)
    
    log.debug("band_info after prescaling: " + str(band_info.keys()))
    
    # Determine grids
    try:
        log.info("Determining what grids the data fits in...")
        grid_jobs = create_grid_jobs(sat, instrument, nav_uid, band_info,
                                     backend_object, cartographer,
                                     fbf_lat=flatbinaryfilename_lat,
                                     fbf_lon=flatbinaryfilename_lon,
                                     forced_grids=forced_grid,
                                     lon_fill_value=lon_fill,
                                     lat_fill_value=lat_fill,)
    except StandardError:
        log.debug("Grid Determination error:", exc_info=1)
        log.error("Determining data's grids failed")
        status_to_return |= STATUS_GDETER_FAIL
        return status_to_return
    
    ### Remap the data
    try:
        remapped_jobs = remap.remap_bands(sat, instrument, nav_uid,
                                          flatbinaryfilename_lon,
                                          flatbinaryfilename_lat, grid_jobs,
                                          num_procs=num_procs,
                                          fornav_d=fornav_d, fornav_D=fornav_D,
                                          lat_fill_value=lat_fill,
                                          lon_fill_value=lon_fill,
                                          do_single_sample=fornav_m,
                                          )
    except StandardError:
        log.debug("Remapping Error:", exc_info=1)
        log.error("Remapping data failed")
        status_to_return |= STATUS_REMAP_FAIL
        return status_to_return
    
    ### BACKEND ###
    W = Workspace('.')
    # for each grid
    for grid_name, grid_dict in remapped_jobs.items():
        
        # process each band for this grid
        for band_kind, band_id in grid_dict.keys():
            
            band_dict = grid_dict[(band_kind, band_id)]
            
            log.info("Running AWIPS backend for %s%s band grid %s" % (band_kind, band_id, grid_name))
            try:
                # Get the data from the flat binary file
                data = getattr(W, band_dict["fbf_remapped"].split(".")[0]).copy()
                
                # Call the backend
                backend_object.create_product(
                                            sat,
                                            instrument,
                                            nav_uid,
                                            band_kind,
                                            band_id,
                                            band_dict["data_kind"],
                                            data,
                                            start_time=start_time,
                                            grid_name=grid_name,
                                            ncml_template=forced_nc or None,
                                            fill_value=band_dict.get("fill_value", None)
                                            )
            except StandardError:
                log.error("Error in the AWIPS backend for %s%s in grid %s" % (band_kind, band_id, grid_name))
                log.debug("AWIPS backend error:", exc_info=1)
                del remapped_jobs[grid_name][(band_kind, band_id)]
        
        # if all the jobs for a grid failed, warn the user and take that grid off the list
        if len(remapped_jobs[grid_name]) == 0 :
            log.error("All backend jobs for grid %s failed" % (grid_name,))
            del remapped_jobs[grid_name]
    
    # if remapping failed for all grids, warn the user
    if len(remapped_jobs) == 0:
        log.warning("AWIPS backend failed for all grids in this data set")
        status_to_return |= STATUS_BACKEND_FAIL
    
    log.info("Processing of data set is complete")
    
    return status_to_return

def _process_data_sets(*args, **kwargs):
    """Wrapper function around `process_data_sets` so that it can be
    called properly from `run_glue`, where the exitcode is the actual
    returned value from `process_data_sets`.

    This function also checks for exceptions other than the ones already
    checked for in `process_data_sets` so that they are
    recorded properly.
    """
    try:
        stat = process_data_sets(*args, **kwargs)
        sys.exit(stat)
    except MemoryError:
        log.error("geocat2awips ran out of memory, check log file for more info")
        log.debug("Memory error:", exc_info=1)
    except OSError:
        log.error("geocat2awips had a OS error, check log file for more info")
        log.debug("OS error:", exc_info=1)
    except StandardError:
        log.error("geocat2awips had an unexpected error, check log file for more info")
        log.debug("Unexpected/Uncaught error:", exc_info=1)
    except KeyboardInterrupt:
        log.info("geocat2awips was cancelled by a keyboard interrupt")
    
    sys.exit(-1)

def run_glue(filepaths,
             multiprocess=True,
             **kwargs):
    """Convert a Geocat hdf file into a AWIPS NetCDF file.
    
    1.  geocat_guidebook.py       : Info on what's in the files
    2.  geocat_to_swath.py        : Code to load the data
    3.  create_grid_jobs          : Figure out what grids the data fits in
    4.  ll2cr
    5.  fornav
    6.  rescale.py
    7.  awips_netcdf.py
    """
    
    # Rewrite/force parameters to specific format
    filepaths = [ os.path.abspath(os.path.expanduser(x)) for x in sorted(filepaths) ]
    
    # figure out which of our files can be processed and which cannot
    all_used          = set([ ])
    for file_pattern in sorted(FILE_CONTENTS_GUIDE.keys()) :
        all_used.update(set([ x for x in filepaths if re.match(file_pattern, os.path.split(x)[1]) ]))
    
    # now figure out what we used and didn't
    all_provided = set(filepaths)
    not_used = all_provided - all_used
    if len(not_used):
        log.warning("Didn't know what to do with\n%s" % "\n".join(list(not_used)))
    
    # some things that we'll use later for clean up
    processes_to_wait_for = [ ]
    exit_status           = 0
    
    # go through and process our files
    for filepath in all_used:
        log.debug("Processing file %s" % filepath)
        filename = os.path.split(filepath)[1]
        nav_uid  = find_nav_uid_from_filename(filename)
        
        try:
            if multiprocess:
                temp_processes = Process(target=_process_data_sets,
                                         args = ([filepath], nav_uid),
                                         kwargs = kwargs
                                         )
                temp_processes.start()
                processes_to_wait_for.append(temp_processes)
            else:
                stat = _process_data_sets([filepath], nav_uid, **kwargs)
                exit_status = exit_status or stat
        except StandardError:
            log.error("Could not process file %s" % filepath, exc_info=True)
            exit_status = exit_status or 1 # TODO, not right
    
    log.debug("Waiting for subprocesses")
    # look through our processes and wait for any processes we saved to wait for
    for each_process in processes_to_wait_for :
        each_process.join()
        stat = each_process.exitcode
        exit_status = exit_status or stat
    
    return exit_status

def main():
    import argparse
    description = """
    Create Geocat swaths, remap them to a grid, and place that remapped data
    into a AWIPS compatible netcdf file.
    """
    
    parser = argparse.ArgumentParser(description=description)
    
    # Logging related
    parser.add_argument('-v', '--verbose', dest='verbosity', action="count", default=0,
            help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG (default INFO)')
    parser.add_argument('-l', '--log', dest="log_fn", default=None,
            help="""specify the log filename, default
<gluescript>_%%Y%%m%%d_%%H%%M%%S. Date information is provided from data filename
through strftime. Current time if no files.""")
    parser.add_argument('--debug', dest="debug_mode", default=False,
            action='store_true',
            help="Enter debug mode. Keeping intermediate files.")
    
    # general behavior
    #parser.add_argument('-h', '--help', dest='help', default=False, action='store_true',
    #        help="Print help describing how to use this command line tool")
    
    # Remapping and grid related
    parser.add_argument('--fornav-D', dest='fornav_D', default=10,
            help="Specify the -D option for fornav")
    parser.add_argument('--fornav-d', dest='fornav_d', default=1,
            help="Specify the -d option for fornav")
    parser.add_argument('--fornav-m', dest='fornav_m', default=False, action='store_true',
            help="Specify the -m option for fornav")
    parser.add_argument('--grid-configs', dest='grid_configs', nargs="+", default=tuple(),
            help="Specify additional grid configuration files ('grids.conf' for built-ins)")
    parser.add_argument('-g', '--grids', dest='forced_grids', nargs="+", default=["all"],
            help="Force remapping to only some grids, defaults to 'all', use 'all' for determination")
    parser.add_argument('--gpd', dest='forced_gpd', default=None,
            help="Specify a different gpd file to use")
    
    # multiprocess related
    parser.add_argument('--sp', dest='single_process', default=False, action='store_true',
            help="Processing is sequential instead of one process per navigation group")
    parser.add_argument('--num-procs', dest="num_procs", default=1,
            help="Specify number of processes that can be used to run ll2cr/fornav calls in parallel")
    
    # Backend Specific
    parser.add_argument('--nc', dest='forced_nc', default=None,
            help="Specify a different ncml file to use")
    parser.add_argument('--backend-config', dest='backend_config', default=None,
            help="specify alternate backend configuration file")
    parser.add_argument('--rescale-config', dest='rescale_config', default=None,
            help="specify alternate rescale configuration file")
    
    # Input related
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', dest='data_files', nargs="+",
            help="List of one or more hdf files")
    group.add_argument('-d', dest='data_dir', nargs="?",
            help="Data directory to look for input data files")
    group.add_argument('-R', dest='remove_prev', default=False, action='store_true',
            help="Delete any files that may conflict with future processing. Processing is not done with this flag.")
    
    args = parser.parse_args()
    
    # if they asked for help, print that and stop now
    #if args.help :
    #    parser.print_help()
    #    sys.exit(0)
    
    # Figure out what the log should be named
    log_fn = args.log_fn
    if args.remove_prev:
        # They didn't need to specify a filename
        if log_fn is None :
            log_fn = GLUE_NAME + "_removal.log"
        file_start_time = utc_now()
    else:
        # Get input files and the first filename for the logging datetime
        if args.data_files:
            hdf_files = args.data_files[:]
        elif args.data_dir:
            base_dir = os.path.abspath(os.path.expanduser(args.data_dir))
            hdf_files = [ os.path.join(base_dir,x) for x in os.listdir(base_dir) ]
        else:
            # Should never get here because argparse mexc group
            log.error("Wrong number of arguments")
            parser.print_help()
            return -1
    
        # Handle the user using a '~' for their home directory
        hdf_files = [ os.path.realpath(os.path.expanduser(x)) for x in sorted(hdf_files) ]
        for hdf_file in hdf_files:
            if not os.path.exists(hdf_file):
                print "ERROR: File '%s' doesn't exist" % (hdf_file,)
                return -1
        
        # Get the date of the first file if provided
        file_start_time = sorted(Geo_Frontend.parse_datetimes_from_filepaths(hdf_files))[0]
    
    # Determine the log filename
    if log_fn is None :
        log_fn = GLUE_NAME + "_%Y%m%d_%H%M%S.log"
    log_fn = datetime.strftime(file_start_time, log_fn)
    
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    setup_logging(log_filename=log_fn, console_level=levels[min(3, args.verbosity)])
    
    # Don't set this up until after you have setup logging
    sys.excepthook = create_exc_handler(GLUE_NAME)
    
    # Remove previous intermediate and product files
    if args.remove_prev:
        log.info("Removing any possible conflicting files")
        remove_file_patterns(
                Geo_Frontend.removable_file_patterns,
                remap.removable_file_patterns,
                Backend.removable_file_patterns
                )
        return 0
    
    fornav_D = int(args.fornav_D)
    fornav_d = int(args.fornav_d)
    fornav_m = args.fornav_m
    num_procs = int(args.num_procs)
    forced_grids = args.forced_grids
    # Assumes 'all' doesn't appear in the list twice
    if 'all' in forced_grids :
        forced_grids[forced_grids.index('all')] = None
    if args.forced_gpd is not None:
        args.forced_gpd = os.path.realpath(os.path.expanduser(args.forced_gpd))
        if not os.path.exists(args.forced_gpd):
            log.error("Specified gpd file does not exist '%s'" % args.forced_gpd)
            return -1
    if args.forced_nc is not None:
        args.forced_nc = os.path.realpath(os.path.expanduser(args.forced_nc))
        if not os.path.exists(args.forced_nc):
            log.error("Specified nc file does not exist '%s'" % args.forced_nc)
            return -1
    
    stat = run_glue(hdf_files,
                    fornav_D=fornav_D, fornav_d=fornav_d, fornav_m=fornav_m,
                    grid_configs=args.grid_configs, #
                    forced_grid=forced_grids, forced_gpd=args.forced_gpd, forced_nc=args.forced_nc,
                    #create_pseudo=args.create_pseudo, # geocat2awips doesn't make pseudo bands
                    multiprocess=not args.single_process, num_procs=num_procs,
                    rescale_config=args.rescale_config, backend_config=args.backend_config
                    )
    log.debug("Processing returned status code: %d" % stat)
    
    # Remove intermediate files (not the backend)
    if not stat and not args.debug_mode:
        log.info("Removing intermediate products")
        remove_file_patterns(
                Geo_Frontend.removable_file_patterns,
                remap.removable_file_patterns
                )
    
    return stat

if __name__ == "__main__":
    sys.exit(main())


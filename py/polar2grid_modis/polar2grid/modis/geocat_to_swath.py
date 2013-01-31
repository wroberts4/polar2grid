#!/usr/bin/env python
# encoding: utf-8
"""
Read HDF4 Geocat products
Write out product binary files used by ms2gt tools.

:author:       Eva Schiffer (evas)
:contact:      evas@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2013 University of Wisconsin SSEC. All rights reserved.
:date:         Sept 2013
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

import geocat_guidebook
from polar2grid.core.constants import *
from polar2grid.core import roles
from .modis_to_swath import array_appender, file_appender

import numpy
from pyhdf.SD import SD,SDC, SDS, HDF4Error

import os
import re
import sys
import logging

log = logging.getLogger(__name__)

class FileInfoObject (object) :
    """
    An object that will automatically compile some info on our file.
    """
    
    def __init__ (self, file_path, load_file=True) :
        """
        Figure out some per-file information.
        """
        
        self.full_path     = file_path
        self.file_name     = os.path.split(file_path)[1]
        self.matching_re   = _get_matching_re(self.file_name)
        self.datetime      = geocat_guidebook.parse_datetime_from_filename(self.file_name)
        
        self.file_object   = None
        if load_file :
            self.file_object = SD(self.full_path, SDC.READ)
        
    
    def get_geo_file (self) :
        """
        get the geonavigation file as an opened file object, open it if needed
        """
        
        if self.file_object is None :
            self.file_object = SD(self.full_path, SDC.READ)
        
        return self.file_object
    
    def close_files (self) :
        """
        close the various files when we're done with them
        """
        
        if self.file_object is not None :
            self.file_object.close()

def _get_matching_re (file_name) :
    """
    given a file, figure out what regular expression matches it in
    the geocat_guidebook.FILE_CONTENTS_GUIDE
    
    WARNING: if the file somehow matches multiple expressions,
    only the last will be returned
    """
    matched_re = None
    
    for file_expression in geocat_guidebook.FILE_CONTENTS_GUIDE :
        
        if re.match(file_expression, file_name) :
            matched_re = file_expression
    
    return matched_re

def _load_meta_data (file_objects) :
    """
    load meta-data from the given list of FileInfoObject's
    """
    
    if len(file_objects) != 1 :
        raise ValueError("One file was expected for processing in _load_meta_data_and_image_data and more were given.")
    file_object = file_objects[0]
    
    # set up the base dictionaries
    meta_data = {
                 "sat": geocat_guidebook.get_satellite_from_filename(file_object.file_name),
                 "instrument": INST_MODIS, # TODO, this is not a given, once more are wired in this will need to be properly selected
                 "start_time": geocat_guidebook.parse_datetime_from_filename(file_object.file_name),
                 "bands" : { },
                 "rows_per_scan": geocat_guidebook.MODIS_ROWS_PER_SCAN, # TODO, not general?
                 
                 # TO FILL IN LATER
                 "lon_fill_value": None,
                 "lat_fill_value": None,
                 "fbf_lat":        None,
                 "fbf_lon":        None,
                 "swath_rows":     None,
                 "swath_cols":     None,
                 "swath_scans":    None,
                }
    
    # pull information on the data that should be in this file
    file_contents_guide = geocat_guidebook.FILE_CONTENTS_GUIDE[file_object.matching_re]
    
    # based on the list of bands/band IDs that should be in the file, load up the meta data and image data
    for band_kind in file_contents_guide.keys() :
        
        for band_number in file_contents_guide[band_kind] :
            
            data_kind_const = geocat_guidebook.DATA_KINDS[(band_kind, band_number)]
            
            # TODO, when there are multiple files, this will algorithm will need to change
            meta_data["bands"][(band_kind, band_number)] = {
                                                            "data_kind": data_kind_const,
                                                            "remap_data_as": data_kind_const,
                                                            "kind": band_kind,
                                                            "band": band_number,
                                                            "rows_per_scan": geocat_guidebook.MODIS_ROWS_PER_SCAN, # TODO not a long term solution
                                                            
                                                            # TO FILL IN LATER
                                                            "fill_value":    None,
                                                            "fbf_img":       None,
                                                            "swath_rows":    None,
                                                            "swath_cols":    None,
                                                            "swath_scans":   None,
                                                            
                                                            # this is temporary so it will be easier to load the data later
                                                            "file_obj":      file_object,
                                                           }
    
    return meta_data

def _load_geonav_data (meta_data_to_update, file_info_objects, nav_uid=None, cut_bad=False) :
    """
    load the geonav data and save it in flat binary files; update the given meta_data_to_update
    with information on where the files are and what the shape and range of the nav data are
    
    TODO, cut_bad currently does nothing
    FUTURE nav_uid will need to be passed once we are using more types of navigation
    """
    
    list_of_geo_files = [ ]
    for file_info in file_info_objects :
        list_of_geo_files.append(file_info.get_geo_file())
    
    # FUTURE, if the longitude and latitude ever have different variable names, this will need refactoring
    lat_temp_file_name, lat_stats = _load_data_to_flat_file (list_of_geo_files, "lat",
                                                             geocat_guidebook.LATITUDE_VARIABLE_NAME,
                                                             scale_name=geocat_guidebook.DEFAULT_SCALE_FACTOR_NAME,
                                                             offset_name=geocat_guidebook.DEFAULT_ADD_OFFSET_NAME,
                                                             scale_method_name=geocat_guidebook.DEFAULT_SCALE_METHOD_NAME,
                                                             missing_attribute_name=geocat_guidebook.DEFAULT_FILL_VALUE_NAME)
    lon_temp_file_name, lon_stats = _load_data_to_flat_file (list_of_geo_files, "lon",
                                                             geocat_guidebook.LONGITUDE_VARIABLE_NAME,
                                                             scale_name=geocat_guidebook.DEFAULT_SCALE_FACTOR_NAME,
                                                             offset_name=geocat_guidebook.DEFAULT_ADD_OFFSET_NAME,
                                                             scale_method_name=geocat_guidebook.DEFAULT_SCALE_METHOD_NAME,
                                                             missing_attribute_name=geocat_guidebook.DEFAULT_FILL_VALUE_NAME)
    
    # rename the flat file to a more descriptive name
    shape_temp = lat_stats["shape"]
    suffix = '.real4.' + '.'.join(str(x) for x in reversed(shape_temp))
    new_lat_file_name = "latitude_"  + suffix # TODO what to use? + str(nav_uid) + suffix 
    new_lon_file_name = "longitude_" + suffix # TODO what to use? + str(nav_uid) + suffix 
    os.rename(lat_temp_file_name, new_lat_file_name)
    os.rename(lon_temp_file_name, new_lon_file_name)
    
    # based on our statistics, save some meta data to our meta data dictionary
    rows, cols = shape_temp
    meta_data_to_update["lon_fill_value"] = lon_stats["fill_value"]
    meta_data_to_update["lat_fill_value"] = lat_stats["fill_value"]
    meta_data_to_update["fbf_lat"]        = new_lat_file_name
    meta_data_to_update["fbf_lon"]        = new_lon_file_name
    meta_data_to_update["swath_rows"]     = rows
    meta_data_to_update["swath_cols"]     = cols
    meta_data_to_update["swath_scans"]    = rows / geocat_guidebook.MODIS_ROWS_PER_SCAN # TODO, not a long term solution
    meta_data_to_update["nav_set_uid"]    = nav_uid

def _load_data_to_flat_file (file_objects, descriptive_string, variable_name,
                             missing_attribute_name=None, fill_value_default=DEFAULT_FILL_VALUE,
                             scale_name=None, offset_name=None, scale_method_name=None,
                             min_fn=numpy.min, max_fn=numpy.max) :
    """
    given a list of file info objects, load the requested variable and append it into a single flat
    binary file with a temporary name based on the descriptive string
    
    the name of the new flat binary file and dictionary of statistics about the data will be returned
    """
    
    #print("variable_name: " + variable_name)
    
    # a couple of temporaries to hold some stats
    minimum_value = None
    maximum_value = None
    fill_value    = None
    
    # open the file with a temporary name and set up a file appender
    temp_id        = str(os.getpid())
    temp_file_name = temp_id + "." + descriptive_string
    temp_flat_file = file(temp_file_name, 'wb')
    temp_appender  = file_appender(temp_flat_file, dtype=numpy.float32) # set to float32 to keep everything consistent
    
    # append in data from each file
    # TODO, these files aren't currently sorted by date
    for file_object in file_objects :
        
        # get the appropriate file and variable object
        temp_var_object = file_object.select(variable_name)
        
        # extract the variable data
        temp_var_data   = temp_var_object[:].astype(numpy.float32)
        
        # figure out where the missing values are
        temp_fill_value = None
        # if we have an attribute name for the fill value then load it, otherwise use the default
        if missing_attribute_name is not None :
            #print ("attributes: " + str(temp_var_object.attributes()))
            
            temp_fill_value = temp_var_object.attributes()[missing_attribute_name]
        else :
            temp_fill_value = fill_value_default
        # if we already have a fill value and it's not the same as the one we just loaded, fix our data
        if (fill_value is not None) and (temp_fill_value != fill_value) :
            temp_var_data[temp_var_data == temp_fill_value] = fill_value
            temp_fill_value = fill_value
        fill_value      = temp_fill_value
        not_fill_mask   = temp_var_data != fill_value
        
        # if there's scaling information load it
        scale_value  = None
        if scale_name  is not None :
            scale_value  = temp_var_object.attributes()[scale_name]
            scale_value  = float(scale_value)  if scale_value  is not None else scale_value
        offset_value = None
        if offset_name is not None :
            offset_value = temp_var_object.attributes()[offset_name]
            offset_value = float(offset_value) if offset_value is not None else offset_value
        scaling_method = geocat_guidebook.SCALING_METHOD_LINEAR
        if scale_method_name is not None :
            scaling_method = temp_var_object.attributes()[scale_method_name]
            scaling_method = int(scaling_method) if scaling_method is not None else scaling_method
        
        log.debug("Using scale method " + str(scaling_method) + " scale value " + str(scale_value) + " and offset value " + str(offset_value))
        
        # unscale our data if needed
        temp_var_data = geocat_guidebook.unscale_data(temp_var_data, fill_value, scaling_method, scale_factor=scale_value, add_offset=offset_value)
        
        # append the file data to the flat file
        temp_appender.append(temp_var_data)
        
        # at this point we need to calculate some statistics based on the data we're saving
        to_use_temp   = numpy.append(temp_var_data[not_fill_mask], minimum_value) if minimum_value is not None else temp_var_data[not_fill_mask]
        minimum_value = min_fn(to_use_temp) if to_use_temp.size > 0 else minimum_value
        to_use_temp   = numpy.append(temp_var_data[not_fill_mask], maximum_value) if maximum_value is not None else temp_var_data[not_fill_mask]
        maximum_value = max_fn(to_use_temp) if to_use_temp.size > 0 else maximum_value
        
        #print ("variable " + str(variable_name) + " has fill value " + str(fill_value) + " and data range " + str(minimum_value) + " to " + str(maximum_value))
    
    # save some statistics to a dictionary
    stats = {
             "shape": temp_appender.shape,
             "min":         minimum_value,
             "max":         maximum_value,
             "fill_value":     fill_value,
            }
    
    # close the flat binary file object (insuring that all appends are flushed to disk)
    temp_flat_file.close()
    
    return temp_file_name, stats

def _load_image_data (meta_data_to_update, cut_bad=False) :
    """
    load image data into binary flat files based on the meta data provided
    """
    
    # process each of the band kind / id sets
    for band_kind, band_id in meta_data_to_update["bands"] :
        
        # load the data into a flat file
        (scale_name, offset_name, scaling_method) = geocat_guidebook.RESCALING_ATTRS[(band_kind, band_id)]
        temp_image_file_name, image_stats = _load_data_to_flat_file ([meta_data_to_update["bands"][(band_kind, band_id)]["file_obj"].file_object],
                                                                     str(band_kind) + str(band_id),
                                                                     geocat_guidebook.VAR_NAMES[(band_kind, band_id)],
                                                                     missing_attribute_name=geocat_guidebook.FILL_VALUE_ATTR_NAMES[(band_kind, band_id)],
                                                                     scale_method_name=scaling_method, scale_name=scale_name, offset_name=offset_name)
        
        # we don't need this entry with the file object anymore, so remove it
        del meta_data_to_update["bands"][(band_kind, band_id)]["file_obj"]
        
        # rename the file with a more descriptive name
        shape_temp = image_stats["shape"]
        suffix = '.real4.' + '.'.join(str(x) for x in reversed(shape_temp))
        new_img_file_name = "image_" + str(band_kind) + "_" + str(band_id) + suffix
        os.rename(temp_image_file_name, new_img_file_name)
        
        # based on our statistics, save some meta data to our meta data dictionary
        rows, cols = shape_temp
        meta_data_to_update["bands"][(band_kind, band_id)]["fill_value"]  = image_stats["fill_value"]
        meta_data_to_update["bands"][(band_kind, band_id)]["fbf_img"]     = new_img_file_name
        meta_data_to_update["bands"][(band_kind, band_id)]["swath_rows"]  = rows
        meta_data_to_update["bands"][(band_kind, band_id)]["swath_cols"]  = cols
        meta_data_to_update["bands"][(band_kind, band_id)]["swath_scans"] = rows / geocat_guidebook.MODIS_ROWS_PER_SCAN
        
        if rows != meta_data_to_update["swath_rows"] or cols != meta_data_to_update["swath_cols"]:
            msg = ("Expected %d rows and %d cols, but band %s %s had %d rows and %d cols"
                   % (meta_data_to_update["swath_rows"], meta_data_to_update["swath_cols"], band_kind, band_id, rows, cols))
            log.error(msg)
            raise ValueError(msg)

def get_swaths(ifilepaths, cut_bad=False, nav_uid=None):
    """Takes geocat hdf files and creates flat binary files for the information
    required to do further processing.

    :Parameters:
        ifilepaths : list
            List of image data filepaths ('*.hdf') of one kind of band that are
            to be concatenated into a swath. These paths should already be user
            expanded absolute paths.
            TODO, This code does not handle concatination of multiple files.
            TODO, For now only one file will be accepted.
    :Keywords:
        cut_bad : bool
            Specify whether or not to delete/cut out entire scans of data
            when navigation data is bad.  This is done because the ms2gt
            utilities used for remapping can't handle incorrect navigation data
            TODO, for now this doesn't do anything!
    """
    
    # TODO, for now this method only handles one file, eventually it will need to handle more
    if len(ifilepaths) != 1 :
        raise ValueError("One file was expected for processing in get_swaths and more were given.")
    
    # make sure the file exists and get minimal info on it
    assert(os.path.exists(ifilepaths[0]))
    file_info = FileInfoObject(ifilepaths[0])
    
    # get the initial meta data information and raw image data
    log.info("Getting data file info...")
    meta_data = _load_meta_data ([file_info])
    
    # load the geonav data and put it in flat binary files
    log.info("Creating binary files for latitude and longitude data")
    _load_geonav_data (meta_data, [file_info], cut_bad=cut_bad, nav_uid=nav_uid)
    
    # load the image data and put it in flat binary files
    log.info("Creating binary files for image data")
    _load_image_data (meta_data, cut_bad=cut_bad)
    
    return meta_data

class Frontend(roles.FrontendRole):
    def __init__(self):
        pass

    def make_swaths(self, filepaths, **kwargs):
        
        # load up all the meta data
        meta_data = { }
        for temp_filepath in filepaths :
            
            try:
                temp_meta_data = get_swaths([temp_filepath], **kwargs)
                temp_bands     = { } if "bands" not in meta_data else meta_data["bands"]
                meta_data.update(temp_meta_data)
                meta_data["bands"].update(temp_bands)
            except StandardError:
                log.error("Swath creation failed")
                log.debug("Swath creation error:", exc_info=1)
        
        return meta_data

def main():
    import optparse
    usage = """
%prog [options] filename1.h,filename2.h,filename3.h,... struct1,struct2,struct3,...

"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                    action="store_true", default=False, help="run self-tests")
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
                    help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    # parser.add_option('-o', '--output', dest='output',
    #                 help='location to store output')
    # parser.add_option('-I', '--include-path', dest="includes",
    #                 action="append", help="include path to append to GCCXML call")
    (options, args) = parser.parse_args()

    # make options a globally accessible structure, e.g. OPTS.
    global OPTS
    OPTS = options

    if options.self_test:
        import doctest
        doctest.testmod()
        sys.exit(0)

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[min(3, options.verbosity)])

    if not args:
        parser.error( 'incorrect arguments, try -h or --help.' )
        return 9

    import json
    meta_data = make_swaths(args[:])
    print json.dumps(meta_data)
    return 0

if __name__=='__main__':
    sys.exit(main())

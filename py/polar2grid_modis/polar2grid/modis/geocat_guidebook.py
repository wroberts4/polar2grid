#!/usr/bin/env python
# encoding: utf-8
"""
Provide information about Geocat product files for a variety of uses.

:author:       Eva Schiffer (evas)
:contact:      evas@ssec.wisc.edu
:organization: Space Science and Engineering Center (SSEC)
:copyright:    Copyright (c) 2013 University of Wisconsin SSEC. All rights reserved.
:date:         Oct 2013
:license:      GNU GPLv3
:revision:     $Id$
"""
__docformat__ = "restructuredtext en"

from polar2grid.core.constants import *

import sys
import re
import logging
from datetime import datetime

LOG = logging.getLogger(__name__)

# geocat scaling method constants
# 0=no scaling, 1=linear, 2=logarithm, 3=square root
SCALING_METHOD_NO_SCALING    = 0
SCALING_METHOD_LINEAR        = 1
SCALING_METHOD_LOGARITHM     = 2
SCALING_METHOD_SQUARE_ROOT   = 3

# the current assumption is that navigation for all
# geocat variables comes from variables stored in the
# same file as the data you're navigating
LONGITUDE_VARIABLE_NAME      = 'pixel_longitude'
LATITUDE_VARIABLE_NAME       = 'pixel_latitude'

DEFAULT_FILL_VALUE_NAME      = '_FillValue'

DEFAULT_ADD_OFFSET_NAME      = 'add_offset'
DEFAULT_SCALE_FACTOR_NAME    = 'scale_factor'
DEFAULT_SCALE_METHOD_NAME    = 'scaling_method'

# TODO, there are notes in the table suggesting this will need a part specific to the instrument (like the "goesr" part)
# TODO, if so, I should find out how that's formatted and generate the different prefixes programatically
IFR_FOG_PROB_VAR_NAME        = 'goesr_fog_nooptprop_IFR_fog_probability'

# this is true for the 1km data, FUTURE: when we get to other kinds, this will need to be more sophisicated
MODIS_ROWS_PER_SCAN          = 10
# TODO, additional values for other cases, this should cover Aqua and Terra, but we also expect Goes-12, Goes-15, SNPP (VIIRS), and Meteosat-9 (SEVIRI)

# a regular expression that will match geocat files
GEOCAT_FILE_PATTERN            = r'geocatL2\..*?\.\d\d\d\d\d\d\d\.\d\d\d\d\d\d\.hdf'

# a mapping between regular expressions to match files and their band_kind and band_id contents
FILE_CONTENTS_GUIDE = {
                        GEOCAT_FILE_PATTERN:                       {
                                                                     BKIND_IFR:   [BID_FOG],
                                                                    },
                      }

# a mapping between bands and their fill value attribute names
FILL_VALUE_ATTR_NAMES = \
            {
              (BKIND_IFR, BID_FOG):           DEFAULT_FILL_VALUE_NAME,
            }

# a mapping between the bands and their data kinds (in the file)
DATA_KINDS = {
              (BKIND_IFR, BID_FOG):           DKIND_FOG,
             }

# a mapping between the bands and the variable names used in the files to hold them
VAR_NAMES  = {
              (BKIND_IFR, BID_FOG):           IFR_FOG_PROB_VAR_NAME,
             }

# a mapping between bands and the names of their scale and offset attributes
RESCALING_ATTRS = \
             {
              (BKIND_IFR, BID_FOG):           (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
             }

def parse_datetime_from_filename (file_name_string) :
    """parse the given file_name_string and create an appropriate datetime object
    that represents the datetime indicated by the file name; if the file name does
    not represent a pattern that is understood, None will be returned
    """
    
    datetime_to_return = None
    
    # there are at least two file name formats to parse here
    if file_name_string.startswith('geocatL2') :
        temp = file_name_string.split('.')
        datetime_to_return = datetime.strptime(temp[2] + temp[3], "%Y%j%H%M%S")
        # TODO, I need to confirm that this is the right format for the date
    
    return datetime_to_return

def get_satellite_from_filename (data_file_name_string) :
    """given a file name, figure out which satellite it's from
    if the file does not represent a known satellite name
    configuration None will be returned
    """
    
    satellite_to_return = None
    
    if   data_file_name_string.find("Aqua")  >= 0 :
        satellite_to_return = SAT_AQUA
    elif data_file_name_string.find("Terra") >= 0 :
        satellite_to_return = SAT_TERRA
    # TODO, there are other types to process, but I need info on how the names will be structured
    
    return satellite_to_return

# TODO, once this is mature, move it into it's own module so it can be a utility function
def unscale_data (data, fill_value, scaling_method, scale_factor=None, add_offset=None) :
    """
    unscale the given data using the methods defined by geocat
    
    the scaling_method corresponds to the scaling method constants created in geocat:
    0=no scaling, 1=linear, 2=logarithm, 3=square root
    
    currently this method can only handle "no scaling" and "linear" but hopefully it
    will support the others in the FUTURE
    
    data is modified in place and fill values will not be changed
    if a scale factor or add offset is given as None (or not given) it will not be applied
    """
    
    to_return = None
    
    # figure out which scaling method to use
    if   scaling_method == SCALING_METHOD_NO_SCALING :
        
        LOG.debug("No scaling required, using existing data.")
        to_return = data
        
    elif scaling_method == SCALING_METHOD_LINEAR :
        
        LOG.debug("Unscaling Geocat data using linear scaling method.")
        to_return = data
        not_fill_mask = to_return != fill_value
        
        # if we found a scale use it to scale the data
        if scale_factor is not None :
            to_return[not_fill_mask] *= scale_factor
        # if we have an offset use it to offset the data
        if add_offset is not None :
            to_return[not_fill_mask] += add_offset
        
    elif scaling_method == SCALING_METHOD_LOGARITHM :
        
        LOG.warn("Unscaling Geocat data using a logarithm method is not yet supported. Unable to unscale data.")
        
    elif scaling_method == SCALING_METHOD_SQUARE_ROOT :
        
        LOG.warn("Unscaling Geocat data using a square root method is not yet supported. Unable to unscale data.")
    
    return to_return

def main():
    import optparse
    from pprint import pprint
    usage = """
%prog [options] filename1.hdf

"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
            help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_option('-r', '--no-read', dest='read_hdf', action='store_false', default=True,
            help="don't read or look for the hdf file, only analyze the filename")
    (options, args) = parser.parse_args()
    
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[min(3, options.verbosity)])
    
    LOG.info("Currently no command line tests are set up for this module.")
    
    """
    if not args:
        parser.error( 'must specify 1 filename, try -h or --help.' )
        return -1
    
    for fn in args:
        try:
            finfo = generic_info(fn)
        except:
            LOG.error("Failed to get info from filename '%s'" % fn, exc_info=1)
            continue
    
        if options.read_h5:
            generic_read(fn, finfo)
            pprint(finfo)
            if finfo["data_kind"] == K_RADIANCE:
                data_shape = str(finfo["data"].shape)
                print "Got Radiance with shape %s" % data_shape
            elif finfo["data_kind"] == K_REFLECTANCE:
                data_shape = str(finfo["data"].shape)
                print "Got Reflectance with shape %s" % data_shape
            elif finfo["data_kind"] == K_BTEMP:
                data_shape = str(finfo["data"].shape)
                print "Got Brightness Temperature with shape %s" % data_shape
            else:
                data_shape = "Unknown data type"
                print "Got %s" % data_shape
            mask_shape = str(finfo["mask"].shape)
            print "Mask was created with shape %s" % mask_shape
        else:
            pprint(finfo)
    
    """

if __name__ == '__main__':
    sys.exit(main())

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

# FUTURE, these patterns may need to be more flexible depending on future input
# fog variable patterns
IFR_FOG_PROB_VAR_PATTERN     = r'.*?_IFR_fog_probability'
LIFR_FOG_PROB_VAR_PATTERN    = r'.*?_LIFR_fog_probability'
MVFR_FOG_PROB_VAR_PATTERN    = r'.*?_MVFR_fog_probability'
# cloud variable patterns
CLOUD_THICKNESS_VAR_PATTERN  = r'.*?_fog_depth'
CLOUD_PHASE_VAR_PATTERN      = r'.*?_cloud_phase'
# ash variable patterns
ASH_HEIGHT_VAR_PATTERN       = r'.*?_ash_top_height'
ASH_MASS_LOADING_VAR_PATTERN = r'.*?_ash_mass_loading'
ASH_EFF_RADIUS_VAR_PATTERN   = r'.*?_ash_effective_radius'
ASH_BTD_11_12_UM_VAR_PATTERN = r'btd1112'
ASH_11_UM_VAR_PATTERN        = r'channel_14_brightness_temperature'
ASH_VISIBLE_VAR_PATTERN      = r'channel_2_reflectance'
# so2 variable patterns
SO2_LOADING_VAR_PATTERN      = r'.*?_So2_Loading'
SO2_MASK_VAR_PATTERN         = r'.*?_so2_mask'

# this is true for the 1km data, FUTURE: when we get to other kinds, this will need to be more sophisicated
MODIS_ROWS_PER_SCAN          = 10
#AVHRR_ROWS_PER_SCAN          = 1 # FUTURE William confirmed that this is 1, but fornav won't accept less than 2 (so if avhrr input is ever needed, this will need to be reconsidered)
MTSAT_ROWS_PER_SCAN          = 2 # confirmed by William
SEVIRI_ROWS_PER_SCAN         = 2 # confirmed by William
GOES_ROWS_PER_SCAN           = 2 # confirmed by William
# FUTURE, do we need additional values for other cases? this should cover Aqua and Terra, but we also expect Goes-12, Goes-15, SNPP (VIIRS), Meteosat-9 (SEVIRI), and MTSAT-2

# more general rows per scan dictionary
ROWS_PER_SCAN = {
                 (SAT_AQUA,   INST_MODIS)   :      MODIS_ROWS_PER_SCAN,
                 (SAT_TERRA,  INST_MODIS)   :      MODIS_ROWS_PER_SCAN,
                 (SAT_GOES13, INST_GIMAGER) :      GOES_ROWS_PER_SCAN,
                 (SAT_GOES14, INST_GIMAGER) :      GOES_ROWS_PER_SCAN,
                 (SAT_GOES15, INST_GIMAGER) :      GOES_ROWS_PER_SCAN,
                 # FUTURE, add the other satellites that we need to know this information for
                 }

# a regular expression that will match geocat files that are in the general geo_nav group
#GEOCAT_FILE_PATTERN            = r'geocatL2\..*?\.\d\d\d\d\d\d\d\.\d\d\d\d\d\d\.hdf'
GEOCAT_FILE_PATTERN            = r'geocatL2\..*?-(?!(Alaska|CONUS)).*?\.\d\d\d\d\d\d\d\.\d\d\d\d\d\d\.hdf'
# a regular expression that will match Alaska specific files
ALASKA_FILE_PATTERN            = r'geocatL2\..*?(Alaska).*?\.\d\d\d\d\d\d\d\.\d\d\d\d\d\d\.hdf'
# a regular expression that will match CONUS specific files
CONUS_FILE_PATTERN             = r'geocatL2\..*?(CONUS).*?\.\d\d\d\d\d\d\d\.\d\d\d\d\d\d\.hdf'

# not sure if this will work this way in the long run
GEO_FILE_GROUPING = {
                      GEO_NAV_UID:      [GEOCAT_FILE_PATTERN],
                      AK_NAV_UID:       [ALASKA_FILE_PATTERN],
                      CONUS_NAV_UID:    [CONUS_FILE_PATTERN],
                    }

# a mapping between regular expressions to match files and their band_kind and band_id contents
FILE_CONTENTS_GUIDE = {
                        GEOCAT_FILE_PATTERN:                      {
                                                                    BKIND_IFR:   [BID_FOG],
                                                                    BKIND_LIFR:  [BID_FOG],
                                                                    BKIND_MVFR:  [BID_FOG],
                                                                    
                                                                    BKIND_CLDT:  [NOT_APPLICABLE],
                                                                    BKIND_CLDP:  [NOT_APPLICABLE],
                                                                    
                                                                    BKIND_ASHH:  [NOT_APPLICABLE],
                                                                    BKIND_ASHM:  [NOT_APPLICABLE],
                                                                    BKIND_ASHE:  [NOT_APPLICABLE],
                                                                    BKIND_ASHB:  [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASH11: [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASHV:  [NOT_APPLICABLE], # has no attrs
                                                                    
                                                                    BKIND_SO2L:  [NOT_APPLICABLE], # not present in current test files
                                                                    BKIND_SO2M:  [NOT_APPLICABLE], # not present in current test files
                                                                   },

                        ALASKA_FILE_PATTERN:                      {
                                                                    BKIND_IFR:   [BID_FOG],
                                                                    BKIND_LIFR:  [BID_FOG],
                                                                    BKIND_MVFR:  [BID_FOG],

                                                                    BKIND_CLDT:  [NOT_APPLICABLE],
                                                                    BKIND_CLDP:  [NOT_APPLICABLE],

                                                                    BKIND_ASHH:  [NOT_APPLICABLE],
                                                                    BKIND_ASHM:  [NOT_APPLICABLE],
                                                                    BKIND_ASHE:  [NOT_APPLICABLE],
                                                                    BKIND_ASHB:  [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASH11: [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASHV:  [NOT_APPLICABLE], # has no attrs

                                                                    BKIND_SO2L:  [NOT_APPLICABLE], # not present in current test files
                                                                    BKIND_SO2M:  [NOT_APPLICABLE], # not present in current test files
                                                                   },

                        CONUS_FILE_PATTERN:                       {
                                                                    BKIND_IFR:   [BID_FOG],
                                                                    BKIND_LIFR:  [BID_FOG],
                                                                    BKIND_MVFR:  [BID_FOG],

                                                                    BKIND_CLDT:  [NOT_APPLICABLE],
                                                                    BKIND_CLDP:  [NOT_APPLICABLE],

                                                                    BKIND_ASHH:  [NOT_APPLICABLE],
                                                                    BKIND_ASHM:  [NOT_APPLICABLE],
                                                                    BKIND_ASHE:  [NOT_APPLICABLE],
                                                                    BKIND_ASHB:  [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASH11: [NOT_APPLICABLE], # has no attrs
                                                                    BKIND_ASHV:  [NOT_APPLICABLE], # has no attrs

                                                                    BKIND_SO2L:  [NOT_APPLICABLE], # not present in current test files
                                                                    BKIND_SO2M:  [NOT_APPLICABLE], # not present in current test files
                                                                   },
                      }

# a mapping between bands and their fill value attribute names
FILL_VALUE_ATTR_NAMES = \
            {
                (BKIND_IFR,   BID_FOG):           DEFAULT_FILL_VALUE_NAME,
                (BKIND_LIFR,  BID_FOG):           DEFAULT_FILL_VALUE_NAME,
                (BKIND_MVFR,  BID_FOG):           DEFAULT_FILL_VALUE_NAME,
                
                (BKIND_CLDT,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                (BKIND_CLDP,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                
                (BKIND_ASHH,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                (BKIND_ASHM,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                (BKIND_ASHE,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                (BKIND_ASHB,  NOT_APPLICABLE):    None,
                (BKIND_ASH11, NOT_APPLICABLE):    None,
                (BKIND_ASHV,  NOT_APPLICABLE):    None,
                
                (BKIND_SO2L,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                (BKIND_SO2M,  NOT_APPLICABLE):    DEFAULT_FILL_VALUE_NAME,
                
            }

# a mapping between the bands and their data kinds (in the file)
DATA_KINDS = {
                (BKIND_IFR,   BID_FOG):           DKIND_PERCENT,
                (BKIND_LIFR,  BID_FOG):           DKIND_PERCENT,
                (BKIND_MVFR,  BID_FOG):           DKIND_PERCENT,
                
                (BKIND_CLDT,  NOT_APPLICABLE):    DKIND_DISTANCE,
                (BKIND_CLDP,  NOT_APPLICABLE):    DKIND_CATEGORY,
                
                (BKIND_ASHH,  NOT_APPLICABLE):    DKIND_DISTANCE,
                (BKIND_ASHM,  NOT_APPLICABLE):    DKIND_M_LOAD, 
                (BKIND_ASHE,  NOT_APPLICABLE):    DKIND_DISTANCE,
                (BKIND_ASHB,  NOT_APPLICABLE):    DKIND_BTEMP, # this is technically brightness temp difference, is using this type ok?
                (BKIND_ASH11, NOT_APPLICABLE):    DKIND_BTEMP,
                (BKIND_ASHV,  NOT_APPLICABLE):    DKIND_REFLECTANCE,
                
                (BKIND_SO2L,  NOT_APPLICABLE):    DKIND_M_LOAD,
                (BKIND_SO2M,  NOT_APPLICABLE):    DKIND_CATEGORY,
             }

# a mapping between the bands and the variable names used in the files to hold them
VAR_PATTERN  = {
                (BKIND_IFR,   BID_FOG):           IFR_FOG_PROB_VAR_PATTERN,
                (BKIND_LIFR,  BID_FOG):           LIFR_FOG_PROB_VAR_PATTERN,
                (BKIND_MVFR,  BID_FOG):           MVFR_FOG_PROB_VAR_PATTERN,
                
                (BKIND_CLDT,  NOT_APPLICABLE):    CLOUD_THICKNESS_VAR_PATTERN,
                (BKIND_CLDP,  NOT_APPLICABLE):    CLOUD_PHASE_VAR_PATTERN,
                
                (BKIND_ASHH,  NOT_APPLICABLE):    ASH_HEIGHT_VAR_PATTERN,
                (BKIND_ASHM,  NOT_APPLICABLE):    ASH_MASS_LOADING_VAR_PATTERN,
                (BKIND_ASHE,  NOT_APPLICABLE):    ASH_EFF_RADIUS_VAR_PATTERN,
                (BKIND_ASHB,  NOT_APPLICABLE):    ASH_BTD_11_12_UM_VAR_PATTERN,
                (BKIND_ASH11, NOT_APPLICABLE):    ASH_11_UM_VAR_PATTERN,
                (BKIND_ASHV,  NOT_APPLICABLE):    ASH_VISIBLE_VAR_PATTERN,
                
                (BKIND_SO2L,  NOT_APPLICABLE):    SO2_LOADING_VAR_PATTERN,
                (BKIND_SO2M,  NOT_APPLICABLE):    SO2_MASK_VAR_PATTERN,
               }

# a mapping between bands and the names of their scale and offset attributes
RESCALING_ATTRS = \
               {
                (BKIND_IFR,   BID_FOG):           (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_LIFR,  BID_FOG):           (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_MVFR,  BID_FOG):           (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                
                (BKIND_CLDT,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_CLDP,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                
                (BKIND_ASHH,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_ASHM,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_ASHE,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_ASHB,  NOT_APPLICABLE):    (None, None, None),
                (BKIND_ASH11, NOT_APPLICABLE):    (None, None, None),
                (BKIND_ASHV,  NOT_APPLICABLE):    (None, None, None),
                
                (BKIND_SO2L,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
                (BKIND_SO2M,  NOT_APPLICABLE):    (DEFAULT_SCALE_FACTOR_NAME, DEFAULT_ADD_OFFSET_NAME, DEFAULT_SCALE_METHOD_NAME),
               }

def find_nav_uid_from_filename (file_name) :
    """
    Given the file name, figure out which navigation group it goes in and return the
    appropriate nav_uid.

    :param file_name:
    :return:
    """

    matched_uid = None
    for nav_uid in GEO_FILE_GROUPING.keys() :
        for file_pattern in GEO_FILE_GROUPING[nav_uid] :
            if re.match(file_pattern, file_name) is not None :
                matched_uid = nav_uid

    if matched_uid is None :
        LOG.warn("Unable to match file " + file_name + " to a navigation group.")

    return matched_uid

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
        # I confirmed with Corey that this is the correct date format
    
    return datetime_to_return

def get_satellite_from_filename (data_file_name_string) :
    """given a file name, figure out which satellite it's from
    if the file does not represent a known satellite name
    configuration None will be returned
    """
    
    satellite_to_return = None
    instrument_to_return = None
    
    if   data_file_name_string.find("Aqua")  >= 0 :
        satellite_to_return  = SAT_AQUA
        instrument_to_return = INST_MODIS
    elif data_file_name_string.find("Terra") >= 0 :
        satellite_to_return  = SAT_TERRA
        instrument_to_return = INST_MODIS
    elif data_file_name_string.find("GOES-13") >= 0 :
        satellite_to_return  = SAT_GOES13
        instrument_to_return = INST_GIMAGER
    elif data_file_name_string.find("GOES-14") >= 0 :
        satellite_to_return  = SAT_GOES14
        instrument_to_return = INST_GIMAGER
    elif data_file_name_string.find("GOES-15") >= 0 :
        satellite_to_return  = SAT_GOES15
        instrument_to_return = INST_GIMAGER
    elif data_file_name_string.find("SNPP") >= 0 :
        satellite_to_return  = SAT_NPP
        instrument_to_return = INST_VIIRS
    elif data_file_name_string.find("Meteosat-9") >= 0 :
        satellite_to_return = SAT_METEO9
        # TODO, what instrument name to use here?
    elif data_file_name_string.find("MTSAT-2") >= 0 :
        satellite_to_return = SAT_MTSAT2
        # TODO, what instrument name to use here?
    
    return satellite_to_return, instrument_to_return

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
    
    to_return = data
    
    # figure out which scaling method to use
    if   scaling_method == SCALING_METHOD_NO_SCALING :
        
        LOG.debug("No scaling required, using existing data.")
        
    elif scaling_method == SCALING_METHOD_LINEAR :
        
        LOG.debug("Unscaling Geocat data using linear scaling method.")
        not_fill_mask = to_return != fill_value
        
        # if we found a scale use it to scale the data
        if (scale_factor is not None) and (scale_factor is not 1.0) :
            to_return[not_fill_mask] *= scale_factor
        # if we have an offset use it to offset the data
        if (add_offset is   not None) and (add_offset   is not 0.0) :
            to_return[not_fill_mask] += add_offset
        
    elif scaling_method == SCALING_METHOD_LOGARITHM :
        
        LOG.warn("Unscaling Geocat data using a logarithm method is not yet supported. Using raw scaled data.")
        
    elif scaling_method == SCALING_METHOD_SQUARE_ROOT :
        
        LOG.warn("Unscaling Geocat data using a square root method is not yet supported. Using raw scaled data.")
    
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

if __name__ == '__main__':
    sys.exit(main())

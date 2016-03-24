#
# LSST Data Management System
# Copyright 2008-2015 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import lsst.pex.config as pexConfig

import esutil

class IngestIndexedReferenceConfig(pexConfig.Config):
    ref_dataset_name = pexConfig.Field(
            dtype = str,
            default = 'cal_ref_cat',
            doc = 'String to pass to the butler to retrieve persisted files.',
    )
    ra_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = "Name of RA column",
    )
    dec_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = "Name of Dec column",
    )
    mag_column_list = pexConfig.ListField(
            dtype = str,
            default = [],
            doc = """The values in the reference catalog are assumed to be in AB magnitudes.
List of column names to use for photometric information.  At least one entry is required."""
    )
    mag_err_column_map = pexConfig.DictField(
            keytype = str,
            itemtype = str,
            default = {},
            doc = "A map of magnitude column name (key) to magnitude error column (value)."
    )
    is_photometric_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = 'Name of column stating if satisfactory for photometric calibration (optional).'
    )
    is_resolved_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = 'Name of column stating if the object is resolved (optional).'
    )
    is_variable_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = 'Name of column stating if the object is measured to be variable (optional).'
    )
    id_name = pexConfig.Field(
            dtype = str,
            default = '',
            doc = 'Name of column to use as an identifier (optional).'
    )
    extra_col_names = pexConfig.ListField(
            dtype = str,
            default = [],
            doc = 'Extra columns to add to the reference catalog.'
    )
    header_lines = pexConfig.Field(
            dtype = int,
            default = 0,
            doc = 'Number of lines to skip when reading the text reference file.'
    )
    colnames = pexConfig.ListField(
            dtype = str,
            default = [],
            doc = """An ordered list of column names to use in ingesting the catalog.  With an empty
list, column names will be discovered from the first line after the skipped header lines."""
    )
    delimiter = pexConfig.Field(
            dtype = str,
            default = ',',
            doc = 'Delimiter to use when reading text reference files.  Comma is default.'
    )
    level = pexConfig.Field(
            dtype = int,
            default = 8,
            doc = 'Default HTM level.  Level 8 gives ~0.08 sq deg per trixel.',
    )

class HtmIndexer(object):
    def __init__(self, depth=8):
        """!Construct the indexer object

        @param[in] depth  depth of the hierarchy to construct
        """
        self.htm = esutil.htm.HTM(depth)
        # HACK need to call intersect first otherwise it segfaults
        _ = self.htm.intersect(1., 2., 0.00001)

    def get_pixel_ids(self, ctrCoord, radius):
        """!Get all shards that touch a circular aperture

        @param[in] ctrCoord  afwCoord.Coord object of the center of the aperture
        @param[in] radius  afwGeom.Angle object of the aperture radius
        @param[out] Return a list of shards, and a boolean array indicating whether the shard touches the
                    boundary (True) or is fully contained (False).
        """
        pixel_id_list = self.htm.intersect(ctrCoord.getRa().asDegrees(), ctrCoord.getDec().asDegrees(),
                                           radius.asDegrees(), inclusive=True)
        covered_pixel_id_list = self.htm.intersect(ctrCoord.getRa().asDegrees(), ctrCoord.getDec().asDegrees(),
                                                   radius.asDegrees(), inclusive=False)
        is_on_boundary = (pixel_id not in covered_pixel_id_list for pixel_id in pixel_id_list)
        return pixel_id_list, is_on_boundary

    def index_points(self, ra_list, dec_list):
        """!Generate trixel ids for each row in an input file

        @param[in] coord_list  List of coord objects to index
        @param[out] A list of pixel ids
        """
        return self.htm.lookup_id(ra_list, dec_list)

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
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
from lsst.afw.image import fluxFromABMag, fluxErrFromABMagErr

import esutil
import numpy


class IngestIndexedReferenceConfig(pexConfig.Config):
    ref_dataset_name = pexConfig.Field(
        dtype=str,
        default='cal_ref_cat',
        doc='String to pass to the butler to retrieve persisted files.',
    )
    ra_name = pexConfig.Field(
        dtype=str,
        default='',
        doc="Name of RA column",
    )
    dec_name = pexConfig.Field(
        dtype=str,
        default='',
        doc="Name of Dec column",
    )
    mag_column_list = pexConfig.ListField(
        dtype=str,
        default=[],
        doc="""The values in the reference catalog are assumed to be in AB magnitudes.
List of column names to use for photometric information.  At least one entry is required."""
    )
    mag_err_column_map = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        default={},
        doc="A map of magnitude column name (key) to magnitude error column (value)."
    )
    is_photometric_name = pexConfig.Field(
        dtype=str,
        default='',
        doc='Name of column stating if satisfactory for photometric calibration (optional).'
    )
    is_resolved_name = pexConfig.Field(
        dtype=str,
        default='',
        doc='Name of column stating if the object is resolved (optional).'
    )
    is_variable_name = pexConfig.Field(
        dtype=str,
        default='',
        doc='Name of column stating if the object is measured to be variable (optional).'
    )
    id_name = pexConfig.Field(
        dtype=str,
        default='',
        doc='Name of column to use as an identifier (optional).'
    )
    extra_col_names = pexConfig.ListField(
        dtype=str,
        default=[],
        doc='Extra columns to add to the reference catalog.'
    )
    header_lines = pexConfig.Field(
        dtype=int,
        default=0,
        doc='Number of lines to skip when reading the text reference file.'
    )
    colnames = pexConfig.ListField(
        dtype=str,
        default=[],
        doc="""An ordered list of column names to use in ingesting the catalog.  With an empty
list, column names will be discovered from the first line after the skipped header lines."""
    )
    delimiter = pexConfig.Field(
        dtype=str,
        default=',',
        doc='Delimiter to use when reading text reference files.  Comma is default.'
    )
    level = pexConfig.Field(
        dtype=int,
        default=8,
        doc='Default HTM level.  Level 8 gives ~0.08 sq deg per trixel.',
    )

    def validate(self):
        pexConfig.Config.validate(self)
        if not (self.ra_name and self.dec_name and self.mag_column_list):
            raise ValueError("RA column name, Dec column name and at least one magnitude column must be" +
                             " supplied.")
        if len(self.mag_err_column_map) > 0 and not len(self.mag_column_list) == len(self.mag_err_column_map):
            raise ValueError("If magnitude errors are provided, all magnitudes must have an error column")


class IngestReferenceRunner(pipeBase.TaskRunner):
    """!Task runner for the reference catalog ingester
    """
    def run(self, parsedCmd):
        """!Run the task.
        Several arguments need to be collected to send on to the task methods.

        @param[in] parsedCmd  Parsed command including command line arguments.
        @param[out] Struct containing the result of the indexing.
        """
        files = parsedCmd.files
        butler = parsedCmd.butler
        task = self.TaskClass(config=self.config, log=self.log, butler=butler)
        task.writeConfig(parsedCmd.butler, clobber=self.clobberConfig, doBackup=self.doBackup)

        result = task.create_indexed_catalog(files)
        if self.doReturnResults:
            return pipeBase.Struct(
                result = result,
            )


class IngestIndexedReferenceTask(pipeBase.CmdLineTask):
    """!Class for both producing indexed reference catalogs and for loading them.

    This implements an indexing scheme based on hierarchical triangular mesh (HTM).
    The term index really means breaking the catalog into localized chunks called
    shards.  In this case each shard contains the entries from the catalog in a single
    HTM trixel
    """
    ConfigClass = IngestIndexedReferenceConfig
    RunnerClass = IngestReferenceRunner
    _DefaultName = 'IngestIndexedReferenceTask'

    _flags = ['photometric', 'resolved', 'variable']

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_argument("files", nargs="+", help="Names of files to index")
        return parser

    def __init__(self, *args, **kwargs):
        """!Constructor for the HTM indexing engine

        @param[in] butler  dafPersistence.Butler object for reading and writing catalogs
        """
        self.butler = kwargs.pop('butler')
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.indexer = HtmIndexer(self.config.level)

    def create_indexed_catalog(self, files):
        """!Index a set of files comprising a reference catalog.  Outputs are persisted in the
        data repository.

        @param[in] files  A list of file names to read.
        """
        rec_num = 0
        first = True
        for filename in files:
            names = True
            if self.config.colnames:
                names = self.config.colnames
            arr = numpy.genfromtxt(filename, dtype=None, skip_header=self.config.header_lines,
                                   delimiter=self.config.delimiter,
                                   names=names)
            index_list = self.indexer.index_points(arr[self.config.ra_name], arr[self.config.dec_name])
            if first:
                schema = self.make_schema(arr.dtype)
                # persist empty catalog to hold the master schema
                dataId = self.make_data_id('master_schema')
                self.butler.put(self.get_catalog(dataId, schema), self.config.ref_dataset_name,
                                dataId=dataId)
                first = False
            pixel_ids = set(index_list)
            for pixel_id in pixel_ids:
                dataId = self.make_data_id(pixel_id)
                catalog = self.get_catalog(dataId, schema)
                els = numpy.where(index_list == pixel_id)
                # Just in case someone has only one line in the file.
                for row in numpy.atleast_1d(arr)[els]:
                    record = catalog.addNew()
                    rec_num = self._fill_record(record, row, rec_num)
                self.butler.put(catalog, self.config.ref_dataset_name, dataId=dataId)

    @staticmethod
    def make_data_id(pixel_id):
        """!Make a data id.  Meant to be overridden.
        @param[in] pixel_id  An identifier for the pixel in question.
        @param[out] dataId (dictionary)
        """
        return {'pixel_id': pixel_id}

    @staticmethod
    def compute_coord(row, ra_name, dec_name):
        """!Create a afwCoord object from a numpy.array row
        @param[in] row  dict like object with ra/dec info in degrees
        @param[in] ra_name  name of RA key
        @param[in] dec_name  name of Dec key
        @param[out] IcrsCoord object constructed from the RA/Dec values
        """
        return afwCoord.IcrsCoord(row[ra_name]*afwGeom.degrees,
                                  row[dec_name]*afwGeom.degrees)

    def _set_flags(self, record, row):
        """!Set the flags for a record.  Relies on the _flags class attribute
        @param[in,out] record  SourceCatalog record to modify
        @param[in] row  dict like object containing flag info
        """
        names = record.schema.getNames()
        for flag in self._flags:
            if flag in names:
                attr_name = 'is_{}_name'.format(flag)
                record.set(flag, bool(row[getattr(self.config, attr_name)]))

    def _set_mags(self, record, row):
        """!Set the flux records from the input magnitudes
        @param[in,out] record  SourceCatalog record to modify
        @param[in] row  dict like object containing magnitude values
        """
        for item in self.config.mag_column_list:
            record.set(item+'_flux', fluxFromABMag(row[item]))
        if len(self.config.mag_err_column_map) > 0:
            for err_key in self.config.mag_err_column_map.keys():
                error_col_name = self.config.mag_err_column_map[err_key]
                record.set(err_key+'_fluxSigma', fluxErrFromABMagErr(row[error_col_name], row[err_key]))

    def _set_extra(self, record, row):
        """!Copy the extra column information to the record
        @param[in,out] record  SourceCatalog record to modify
        @param[in] row  dict like object containing the column values
        """
        for extra_col in self.config.extra_col_names:
            record.set(extra_col, row[extra_col])

    def _fill_record(self, record, row, rec_num):
        """!Fill a record to put in the persisted indexed catalogs

        @param[in,out] record  afwTable.SourceRecord in a reference catalog to fill.
        @param[in] row  A row from a numpy array constructed from the input catalogs.
        """
        record.setCoord(self.compute_coord(row, self.config.ra_name, self.config.dec_name))
        if self.config.id_name:
            record.setId(row[self.config.id_name])
        else:
            rec_num += 1
            record.setId(rec_num)
        # No parents
        record.setParent(-1)

        self._set_flags(record, row)
        self._set_mags(record, row)
        self._set_extra(record, row)
        return rec_num

    def get_catalog(self, dataId, schema):
        """!Get a catalog from the butler or create it if it doesn't exist

        @param[in] dataId  Identifier for catalog to retrieve
        @param[in] schema  Schema to use in catalog creation if the butler can't get it
        @param[out] afwTable.SourceCatalog for the specified identifier
        """
        if self.butler.datasetExists(self.config.ref_dataset_name, dataId=dataId):
            return self.butler.get(self.config.ref_dataset_name, dataId=dataId)
        return afwTable.SourceCatalog(schema)

    def make_schema(self, dtype):
        """!Make the schema to use in constructing the persisted catalogs.

        @param[in] dtype  A numpy.dtype to use in constructing the schema
        @param[out] The schema for the output source catalog.
        """
        mag_column_list = self.config.mag_column_list
        mag_err_column_map = self.config.mag_err_column_map
        if len(mag_err_column_map) > 0 and (
            not len(mag_column_list) == len(mag_err_column_map) or
                not sorted(mag_column_list) == sorted(mag_err_column_map.keys())):
            raise ValueError("Every magnitude column must have a corresponding error column")
        # makes a schema with a coord, id and parent_id
        schema = afwTable.SourceTable.makeMinimalSchema()

        def add_field(name):
            if dtype[name].kind == 'S':
                # dealing with a string like thing.  Need to get type and size.
                at_type = afwTable.aliases[str]
                at_size = dtype[name].itemsize
                return schema.addField(name, type=at_type, size=at_size)
            else:
                at_type = afwTable.aliases[dtype[name].type]
                return schema.addField(name, at_type)

        for item in mag_column_list:
            schema.addField(item+'_flux', float)
        if len(mag_err_column_map) > 0:
            for err_item in mag_err_column_map.keys():
                schema.addField(err_item+'_fluxSigma', float)
        for flag in self._flags:
            attr_name = 'is_{}_name'.format(flag)
            if getattr(self.config, attr_name):
                schema.addField(flag, 'Flag')
        for col in self.config.extra_col_names:
            add_field(col)
        return schema


class HtmIndexer(object):
    def __init__(self, depth=8):
        """!Construct the indexer object

        @param[in] depth  depth of the hierarchy to construct
        """
        self.htm = esutil.htm.HTM(depth)
        # HACK need to call intersect first otherwise it segfaults
        self.htm.intersect(1., 2., 0.00001)

    def get_pixel_ids(self, ctrCoord, radius):
        """!Get all shards that touch a circular aperture

        @param[in] ctrCoord  afwCoord.Coord object of the center of the aperture
        @param[in] radius  afwGeom.Angle object of the aperture radius
        @param[out] A pipeBase.Struct with the list of shards, shards, and a boolean arry, boundary_mask,
                    indicating whether the shard touches the boundary (True) or is fully contained (False).
        """
        pixel_id_list = self.htm.intersect(ctrCoord.getRa().asDegrees(), ctrCoord.getDec().asDegrees(),
                                           radius.asDegrees(), inclusive=True)
        covered_pixel_id_list = self.htm.intersect(ctrCoord.getRa().asDegrees(),
                                                   ctrCoord.getDec().asDegrees(),
                                                   radius.asDegrees(), inclusive=False)
        is_on_boundary = (pixel_id not in covered_pixel_id_list for pixel_id in pixel_id_list)
        return pixel_id_list, is_on_boundary

    def index_points(self, ra_list, dec_list):
        """!Generate trixel ids for each row in an input file

        @param[in] coord_list  List of coord objects to index
        @param[out] A list of pixel ids
        """
        return self.htm.lookup_id(ra_list, dec_list)

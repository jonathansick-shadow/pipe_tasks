#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2016 LSST/AURA
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy
from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Task
from lsst.afw.geom import Box2D


class SetPrimaryFlagsConfig(Config):
    nChildKeyName = Field(dtype=str, default="deblend_nChild",
                          doc="Name of field in schema with number of deblended children")
    pseudoFilterList = ListField(dtype=str, default=['sky'],
                                 doc="Names of filters which should never be primary")


class SetPrimaryFlagsTask(Task):
    ConfigClass = SetPrimaryFlagsConfig

    def __init__(self, schema, **kwargs):
        Task.__init__(self, **kwargs)
        self.schema = schema
        self.isPatchInnerKey = self.schema.addField(
            "detect_isPatchInner", type="Flag",
            doc="true if source is in the inner region of a coadd patch",
        )
        self.isTractInnerKey = self.schema.addField(
            "detect_isTractInner", type="Flag",
            doc="true if source is in the inner region of a coadd tract",
        )
        self.isPrimaryKey = self.schema.addField(
            "detect_isPrimary", type="Flag",
            doc="true if source has no children and is in the inner region of a coadd patch "
                + "and is in the inner region of a coadd tract "
                  "and is not \"detected\" in a pseudo-filter (see config.pseudoFilterList)",
        )

    def run(self, sources, skyMap, tractInfo, patchInfo, includeDeblend=True):
        """Set is-primary and related flags on sources

        @param[in,out] sources   a SourceTable
            - reads centroid fields and an nChild field
            - writes is-patch-inner, is-tract-inner and is-primary flags
        @param[in] skyMap   sky tessellation object (subclass of lsst.skymap.BaseSkyMap)
        @param[in] tractInfo   tract object (subclass of lsst.skymap.TractInfo)
        @param[in] patchInfo   patch object (subclass of lsst.skymap.PatchInfo)
        @param[in] includeDeblend   include deblend information in isPrimary?
        """
        nChildKey = None
        if includeDeblend:
            nChildKey = self.schema.find(self.config.nChildKeyName).key

        # set inner flags for each source and set primary flags for sources with no children
        # (or all sources if deblend info not available)
        innerFloatBBox = Box2D(patchInfo.getInnerBBox())

        # When the centroider fails, we can still fall back to the peak, but we don't trust
        # that quite as much - so we use a slightly smaller box for the patch comparison.
        # That's trickier for the tract comparison, so we just use the peak without extra
        # care there.
        shrunkInnerFloatBBox = Box2D(innerFloatBBox)
        shrunkInnerFloatBBox.grow(-1)

        pseudoFilterKeys = []
        for filt in self.config.pseudoFilterList:
            try:
                pseudoFilterKeys.append(self.schema.find("merge_peak_%s" % filt).getKey())
            except Exception:
                self.log.warn("merge_peak is not set for pseudo-filter %s" % filt)

        tractId = tractInfo.getId()
        for source in sources:
            centroidPos = source.getCentroid()
            if numpy.any(numpy.isnan(centroidPos)):
                continue
            if source.getCentroidFlag():
                # Use a slightly smaller box to guard against bad centroids (see above)
                isPatchInner = shrunkInnerFloatBBox.contains(centroidPos)
            else:
                isPatchInner = innerFloatBBox.contains(centroidPos)
            source.setFlag(self.isPatchInnerKey, isPatchInner)

            skyPos = source.getCoord()
            sourceInnerTractId = skyMap.findTract(skyPos).getId()
            isTractInner = sourceInnerTractId == tractId
            source.setFlag(self.isTractInnerKey, isTractInner)

            if nChildKey is None or source.get(nChildKey) == 0:
                for pseudoFilterKey in pseudoFilterKeys:
                    if source.get(pseudoFilterKey):
                        isPseudo = True
                        break
                else:
                    isPseudo = False

                source.setFlag(self.isPrimaryKey, isPatchInner and isTractInner and not isPseudo)


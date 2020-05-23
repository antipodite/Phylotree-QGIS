# -*- coding: utf-8 -*-

"""
/***************************************************************************
 PhyloTree
                                 A QGIS plugin
 Create, draw and link a phylogenetic tree to vector features
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2020-05-08
        copyright            : (C) 2020 by Isaac Stead
        email                : isaac.stead@protonmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Isaac Stead'
__date__ = '2020-05-08'
__copyright__ = '(C) 2020 by Isaac Stead'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis import processing
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterFile,
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsPoint,
                       QgsPointXY,
                       QgsLineString,
                       QgsGeometry,
                       QgsWkbTypes)
from os.path import splitext
from phylo_tree.trees import drawtree

class PhyloTreeAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT = 'OUTPUT'
    INPUTLAYER = 'INPUTLAYER'
    INPUTTREE = 'INPUTTREE'
    OUT_FIELDS = {
        'id':     QVariant.Int,
        'label':  QVariant.String,
    }
    SCALE_X = 6.0
    SCALE_Y = 8.0

    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUTLAYER,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )
        # Add input tree file source
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUTTREE,
                self.tr('Tree file')
            )
        )
        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output layer')
            )
        )
    
    def processAlgorithm(self, parameters, context, feedback):
        
        layer = self.parameterAsSource(parameters, self.INPUTLAYER, context)
        fname = self.parameterAsFile(parameters, self.INPUTTREE, context)

        # Set up fields for the output layer
        out_fields = QgsFields()
        for name, typ in self.OUT_FIELDS.items():
            out_fields.append(QgsField(name, typ))
        
        (sink, dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT, context, out_fields,
            QgsWkbTypes.LineString, layer.sourceCrs()
        )

        feedback.pushConsoleInfo(str(layer.wkbType()))

        tree = drawtree.buildtree(fname)
        center = (115, -33)
        tree.scale(self.SCALE_X, self.SCALE_Y)
        tree.translate(center)

        # Draw the tree on the map
        polylines = self.create_line_tree(tree, out_fields)
        sink.addFeatures(polylines, QgsFeatureSink.FastInsert)

        # Link tree to input layer features
        features = layer.getFeatures()
        links = self.link_leaves(tree, features, 'Language')
        sink.addFeatures(links, QgsFeatureSink.FastInsert)

        return {self.OUTPUT: dest_id}

    def position_tree(self, tree, inputlayer):
        """
        Using the convex hull of the points in the input layer,
        determine the best place to draw the tree.
        """
        pass

    def link_leaves(self, tree, feats, fieldname):
        """
        Create Polylines linking leaves of tree to input layer
        features
        """
        # Match leaf nodes up with features in input layer
        matches = []
        leaf_table = {leaf.name: leaf for leaf in tree.leaves()}
        for f in feats:
            try:
                name = f[fieldname]
                leaf = leaf_table[name]
                matches.append( (leaf, f) )
            except KeyError:
                pass
        # Create lines linking each pair
        out = []
        for leaf, feat in matches:
            start = QgsPointXY(leaf.x, leaf.y)
            end   = feat.geometry().asPoint()
            line  = QgsGeometry.fromPolylineXY([start, end])
            feat  = QgsFeature()
            feat.setGeometry(line)
            out.append(feat)
        return out

    def create_square_tree(self, tree, fields):
        linedata = tree.construct_squaretree()
        polylines = []
        for line in linedata:
            name, start, end = line
            x1, y1 = start
            x2, y2 = end
            startp = QgsPointXY(x1, y1)
            endp   = QgsPointXY(x2, y2)
            geom   = QgsGeometry.fromPolylineXY([startp, endp])
            feat   = QgsFeature(fields)
            feat.setGeometry(geom)
            feat['label'] = name
            polylines.append(feat)
        return polylines

    def create_line_tree(self, tree, fields):
        out = []
        for i, node in enumerate(tree.walk()):
            if node.parent:
                # Tree geometry
                start = QgsPointXY(node.parent.x, node.parent.y)
                end   = QgsPointXY(node.x, node.y)
                line  = QgsGeometry.fromPolylineXY([start, end])
                feat  = QgsFeature(fields)
                feat.setGeometry(line)
                # Tree labels
                feat['id'], feat['label'] = i, node.name
                out.append(feat)
        return out
    
    def create_point_tree(self, tree):
        out = []
        for i, node in enumerate(tree.walk()):
            feat = QgsFeature(out_fields)
            feat['id'] = i
            feat['label'] = node.name
            if node.parent:
                feat['parent'] = node.parent.name
            feat['startx'] = node.x
            feat['starty'] = node.y
            point = QgsPointXY(node.x, node.y)
            points.append(point) # Save refs to points so can use in next step
            geom = QgsGeometry.fromPointXY(point)
            feat.setGeometry(geom)
            out.append(feat)
        return out

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Create and link phylogenetic tree'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return ''

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return PhyloTreeAlgorithm()

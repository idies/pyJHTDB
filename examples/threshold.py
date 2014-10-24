########################################################################
#
#  Copyright 2014 Johns Hopkins University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Contact: turbulence@pha.jhu.edu
# Website: http://turbulence.pha.jhu.edu/
#
########################################################################

## modified version of lorenz_ui.py taken from
## http://docs.enthought.com/mayavi/mayavi/auto/example_lorenz_ui.html#example-lorenz-ui

import numpy as np

from traits.api import HasTraits, Range, Instance, on_trait_change, Array, Tuple, Str, Button
from traitsui.api import View, Item, HSplit, Group

from mayavi import mlab
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor

import pyJHTDB
import pyJHTDB.dbinfo
info = pyJHTDB.dbinfo.isotropic1024coarse

lTDB = pyJHTDB.libJHTDB()

################################################################################
# `plotter` class.
################################################################################
class plotter(HasTraits):

    # The parameters for the point.
    ox = Range(
            0, 1023, 0,
            desc='first corner x',
            enter_set=True,
            auto_set=False)
    oy = Range(
            0, 1023, 0,
            desc='first corner y',
            enter_set=True,
            auto_set=False)
    oz = Range(
            0, 1023, 0,
            desc='first corner z',
            enter_set=True,
            auto_set=False)

    nx = Range(
            4, 16, 4,
            desc='box size along x',
            enter_set=True,
            auto_set=False)
    ny = Range(
            4, 16, 4,
            desc='box size along y',
            enter_set=True,
            auto_set=False)
    nz = Range(
            4, 16, 4,
            desc='box size along z',
            enter_set=True,
            auto_set=False)
    threshold = Range(
            .0, 200., 15.,
            desc='threshold value',
            enter_set=True,
            auto_set=False)
    time = Range(
            float(info['time'][0]), float(info['time'][-1]), .0,
            desc='time',
            enter_set=True,
            auto_set=False)

    compute = Button(label = 'Compute')

    pcoords = np.zeros((1, 3), dtype = np.float32)

    # The mayavi(mlab) scene.
    scene = Instance(MlabSceneModel, args=())

    # The coordinate system
    xaxis = Instance(HasTraits)
    yaxis = Instance(HasTraits)
    zaxis = Instance(HasTraits)

    # The vectors
    vectors = Instance(HasTraits)

    ########################################
    # The UI view to show the user.
    view = View(HSplit(
                    Group(
                        Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                             height=500, width=500, show_label=False)),
                    Group(
                        Item('ox'),
                        Item('oy'),
                        Item('oz'),
                        Item('nx'),
                        Item('ny'),
                        Item('nz'),
                        Item('threshold'),
                        Item('time'),
                        Item('compute'))
                    ),
                resizable=True
                )

    ######################################################################
    # Trait handlers.
    ######################################################################

    # Note that in the `on_trait_change` call below we listen for the
    # `scene.activated` trait.  This conveniently ensures that the flow
    # is generated as soon as the mlab `scene` is activated (which
    # happens when the configure/edit_traits method is called).  This
    # eliminates the need to manually call the `update_flow` method etc.
    @on_trait_change('compute, scene.activated')
    def update_vectors(self):
        lTDB.initialize()
        bla = lTDB.getThreshold(
                threshold = self.threshold,
                time = self.time,
                field = 'vorticity',
                data_set = info['name'],
                cx = self.ox, cy = self.oy, cz = self.oz,
                nx = self.nx, ny = self.ny, nz = self.nz)
        points = np.zeros((bla.shape[0], 3), np.float32)
        points[:, 0] = bla['x'] * info['dx']
        points[:, 1] = bla['y'] * info['dy']
        points[:, 2] = bla['z'] * info['dz']
        gradu = lTDB.getData(
                self.time,
                points,
                sinterp = 60,
                tinterp = 0,
                data_set = info['name'],
                getFunction = 'getVelocityGradient')
        lTDB.finalize()
        vorticity = points.copy()
        vorticity[:, 0] = gradu[:, 7] - gradu[:, 5]
        vorticity[:, 1] = gradu[:, 2] - gradu[:, 6]
        vorticity[:, 2] = gradu[:, 3] - gradu[:, 1]
        self.vectors.mlab_source.reset(
                x = points[:, 0],
                y = points[:, 1],
                z = points[:, 2],
                u = vorticity[:, 0],
                v = vorticity[:, 1],
                w = vorticity[:, 2])
        self.scene.mlab.view(focalpoint='auto')
        return None

    @on_trait_change('ox, nx, scene.activated')
    def update_xbox(self):
        self.xaxis.mlab_source.set(x = np.array([self.ox, self.ox + self.nx]) * info['dx'])
        self.yaxis.mlab_source.set(x = np.array([self.ox, self.ox]) * info['dx'])
        self.zaxis.mlab_source.set(x = np.array([self.ox, self.ox]) * info['dx'])
        self.scene.mlab.view(focalpoint='auto')
        return None
    @on_trait_change('oy, ny, scene.activated')
    def update_ybox(self):
        self.yaxis.mlab_source.set(y = np.array([self.oy, self.oy + self.ny]) * info['dy'])
        self.xaxis.mlab_source.set(y = np.array([self.oy, self.oy]) * info['dy'])
        self.zaxis.mlab_source.set(y = np.array([self.oy, self.oy]) * info['dy'])
        self.scene.mlab.view(focalpoint='auto')
        return None
    @on_trait_change('oz, nz, scene.activated')
    def update_zbox(self):
        self.zaxis.mlab_source.set(z = np.array([self.oz, self.oz + self.nz]) * info['dz'])
        self.xaxis.mlab_source.set(z = np.array([self.oz, self.oz]) * info['dz'])
        self.yaxis.mlab_source.set(z = np.array([self.oz, self.oz]) * info['dz'])
        self.scene.mlab.view(focalpoint='auto')
        return None


    ######################################################################
    # Private interface.
    ######################################################################

    def _xaxis_default(self):
        return self.scene.mlab.plot3d(
                np.array([self.ox, self.ox + self.nx]) * info['dx'],
                np.array([self.oy, self.oy          ]) * info['dy'],
                np.array([self.oz, self.oz          ]) * info['dz'],
                color = (1, 0, 0),
                tube_radius = None)
    def _yaxis_default(self):
        return self.scene.mlab.plot3d(
                np.array([self.ox, self.ox          ]) * info['dx'],
                np.array([self.oy, self.oy + self.ny]) * info['dy'],
                np.array([self.oz, self.oz          ]) * info['dz'],
                color = (0, 1, 0),
                tube_radius = None)
    def _zaxis_default(self):
        return self.scene.mlab.plot3d(
                np.array([self.ox, self.ox          ]) * info['dx'],
                np.array([self.oy, self.oy          ]) * info['dy'],
                np.array([self.oz, self.oz + self.nz]) * info['dz'],
                color = (0, 0, 1),
                tube_radius = None)
    def _vectors_default(self):
        lTDB.initialize()
        bla = lTDB.getThreshold(
                threshold = self.threshold,
                time = self.time,
                field = 'vorticity',
                data_set = info['name'],
                cx = self.ox, cy = self.oy, cz = self.oz,
                nx = self.nx, ny = self.ny, nz = self.nz)
        points = np.zeros((bla.shape[0], 3), np.float32)
        points[:, 0] = bla['x'] * info['dx']
        points[:, 1] = bla['y'] * info['dy']
        points[:, 2] = bla['z'] * info['dz']
        gradu = lTDB.getData(
                self.time,
                points,
                sinterp = 60,
                tinterp = 0,
                data_set = info['name'],
                getFunction = 'getVelocityGradient')
        lTDB.finalize()
        vorticity = points.copy()
        vorticity[:, 0] = gradu[:, 7] - gradu[:, 5]
        vorticity[:, 1] = gradu[:, 2] - gradu[:, 6]
        vorticity[:, 2] = gradu[:, 3] - gradu[:, 1]
        arrows = self.scene.mlab.quiver3d(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                vorticity[:, 0],
                vorticity[:, 1],
                vorticity[:, 2])
        return arrows


if __name__ == '__main__':
    # Instantiate the class and configure its traits.
    lor = plotter()
    lor.configure_traits()


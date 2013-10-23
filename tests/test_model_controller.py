import os
import shutil
import unittest
import random
import json
import numpy
import urllib
from pytest import raises
from datetime import datetime, timedelta
from paegan.transport.models.transport import Transport
from paegan.transport.models.behavior import LarvaBehavior
from paegan.transport.particles.particle import Particle
from paegan.location4d import Location4D
from paegan.transport.exceptions import ModelError, DataControllerError
from paegan.utils.asarandom import AsaRandom
from paegan.transport.model_controller import ModelController
from shapely.geometry import Point, Polygon
import os
import pytz
import logging
from paegan.logger.easy_logger import EasyLogger


class ModelControllerTest(unittest.TestCase):

    def setUp(self):
        self.start_lat = 60.75
        self.start_lon = -147
        self.start_depth = 0
        self.num_particles = 4
        self.time_step = 3600
        self.num_steps = 10
        self.start_time = datetime(2012, 8, 1, 00)
        self.transport = Transport(horizDisp=0.05, vertDisp=0.0003)

        self.log = EasyLogger('testlog.txt', level=logging.PROGRESS)

        self.output_path = "/data/lm/output"
        self.cache_path = "/data/lm/cache"
        self.bathy_file = "/data/lm/bathy/ETOPO1_Bed_g_gmt4.grd"
        self.shoreline_path = "/data/lm/shore"

    def tearDown(self):
        self.log.close()

    def test_run_from_point(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_run_from_point")

        models = [self.transport]

        p = Point(self.start_lon, self.start_lat, self.start_depth)

        model = ModelController(geometry=p, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=False, use_shoreline=True,
            time_chunk=10, horiz_chunk=4)

        cache_path = os.path.join(self.cache_path, "test_run_from_point.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=cache_path)

    def test_run_from_point_with_wfs(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_run_from_point_with_wfs")

        models = [self.transport]

        p = Point(self.start_lon, self.start_lat, self.start_depth)

        model = ModelController(geometry=p, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=False, use_shoreline=True,
            time_chunk=10, horiz_chunk=4, shoreline_path='http://geo.asascience.com/geoserver/shorelines/ows', shoreline_feature='shorelines:10m_land_polygons')

        cache_path = os.path.join(self.cache_path, "test_run_from_point.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=cache_path)

    def test_run_from_polygon(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_run_from_polygon")

        models = [self.transport]

        poly = Point(self.start_lon, self.start_lat, self.start_depth).buffer(0.001)

        model = ModelController(geometry=poly, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=False, use_shoreline=True,
            time_chunk=10, horiz_chunk=4)

        cache_path = os.path.join(os.path.dirname(__file__), "..", "paegan/transport/_cache/test_run_from_polygon.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=cache_path)

    def test_interp(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_interp")

        models = [self.transport]

        num_steps = 100

        output_path = os.path.join(self.output_path, "test_interp")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=self.num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=10, horiz_chunk=4)

        cache_path = os.path.join(self.cache_path, "test_interp.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

    def test_nearest(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_nearest")

        models = [self.transport]
        
        num_steps = 100

        output_path = os.path.join(self.output_path, "test_nearest")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=10, horiz_chunk=4, time_method='nearest')

        cache_path = os.path.join(self.cache_path, "test_nearest.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

    def test_start_on_land(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_start_on_land")

        # Set the start position and time for the models
        start_lat = 60.15551950079041
        start_lon = -148.1999130249019

        models = [self.transport]

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=False, use_shoreline=True,
            time_chunk=10, horiz_chunk=4, time_method='nearest')

        cache_path = os.path.join(self.cache_path, "test_start_on_land.nc")

        with raises(ModelError):
            model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", cache=cache_path)

    def test_bad_dataset(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_bad_dataset")

        models = [self.transport]

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=self.num_steps, npart=self.num_particles, models=models, use_bathymetry=False, use_shoreline=True,
            time_chunk=10, horiz_chunk=4, time_method='nearest')

        cache_path = os.path.join(self.cache_path, "test_bad_dataset.nc")
        
        with raises(DataControllerError):
            model.run("http://asascience.com/thisisnotadataset.nc", cache=cache_path)

    """
    def test_behavior_growth_and_settlement(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_behavior_growth_and_settlement")

        # 6 days
        num_steps = 144

        num_particles = 2

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/behavior_for_run_testing.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=4, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_behavior_growth_and_settlement")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline','Pickle']

        cache_path = os.path.join(self.cache_path, "test_behavior_growth_and_settlement.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
    """

    def test_quick_settlement(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_quick_settlement")

        num_steps = 24

        num_particles = 4

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/behavior_quick_settle.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=12, horiz_chunk=2, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_quick_settlement")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "test_quick_settlement.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

    def test_timechunk_greater_than_timestep(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_timechunk_greater_than_timestep")

        # 6 days
        num_steps = 10

        num_particles = 2

        models = [self.transport]

        model = ModelController(latitude=self.start_lat, longitude=self.start_lon, depth=self.start_depth, start=self.start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=2)

        cache_path = os.path.join(self.cache_path, "test_timechunk_greater_than_timestep.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path)

    """
    def test_kayak_island(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: Kayak Island")

        # 6 days
        num_steps = 1632

        num_particles = 100

        time_step = 3600

        behavior_config = json.loads(urllib.urlopen("http://behaviors.larvamap.asascience.com/library/50ef1bb1cc7b61000700001d.json").read())
        lb = LarvaBehavior(data=behavior_config[u'results'][0])

        #behavior_config = json.loads(open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/nick.json"))).read())
        #lb = LarvaBehavior(data=behavior_config[u'results'][0])

        models = [Transport(horizDisp=0.01, vertDisp=0.001)]
        models.append(lb)

        start_time = datetime(2011, 5, 2, 00, tzinfo=pytz.utc)

        start_lat = 59.93517413488866
        start_lon = -144.496213677788
        depth = -1

        shoreline_path = os.path.join(self.shoreline_path, "westcoast", "New_Land_Clean.shp")

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=depth, start=start_time, step=time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=5, time_method='interp', shoreline_path=shoreline_path)

        output_path = os.path.join(self.output_path, "kayak_island")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']
        
        cache_path = os.path.join(self.cache_path, "kayak_island.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L1_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
        
    def test_sheep_bay(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: Sheep Bay")

        # 6 days
        num_steps = 1632

        num_particles = 100

        time_step = 3600

        behavior_config = json.loads(urllib.urlopen("http://behaviors.larvamap.asascience.com/library/50ef1bb1cc7b61000700001d.json").read())
        lb = LarvaBehavior(data=behavior_config[u'results'][0])

        models = [Transport(horizDisp=0.01, vertDisp=0.001)]
        models.append(lb)

        start_time = datetime(2011, 5, 2, 00, tzinfo=pytz.utc)

        start_lat = 60.60899655733162
        start_lon = -145.97402533055956
        depth = -1

        shoreline_path = os.path.join(self.shoreline_path, "westcoast", "New_Land_Clean.shp")

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=depth, start=start_time, step=time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=5, time_method='interp', shoreline_path=shoreline_path)

        output_path = os.path.join(self.output_path, "sheep_bay")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']
        
        cache_path = os.path.join(self.cache_path, "sheep_bay.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L2_FCST.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)
    """

    def test_diel_migration(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: test_diel_migration")

        num_steps = 168

        num_particles = 4

        start_time = datetime(2013,4,1,0, tzinfo=pytz.utc)

        # Behavior
        behavior_config = open(os.path.normpath(os.path.join(os.path.dirname(__file__),"./resources/files/diel_suncycles.json"))).read()
        lb = LarvaBehavior(json=behavior_config)

        models = [self.transport]
        models.append(lb)

        model = ModelController(latitude=60.68, longitude=-146.42, depth=self.start_depth, start=start_time, step=self.time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=24, horiz_chunk=2, time_method='nearest')

        output_path = os.path.join(self.output_path, "test_diel_migration")
        shutil.rmtree(output_path, ignore_errors=True)
        os.makedirs(output_path)
        output_formats = ['Shapefile','NetCDF','Trackline']

        cache_path = os.path.join(self.cache_path, "test_diel_migration.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_DAS.nc", bathy=self.bathy_file, cache=cache_path, output_path=output_path, output_formats=output_formats)

    def test_PWS_L0(self):
        self.log.logger.info("**************************************")
        self.log.logger.info("Running: PWS_L0")

        # 1 days
        num_steps = 24

        num_particles = 4

        time_step = 3600

        models = [Transport(horizDisp=1.0, vertDisp=0.5)]

        start_time = datetime(2011, 2, 1, 06, tzinfo=pytz.utc)

        start_lat = 50.989393258377746
        start_lon = -149.75202941894338
        depth = -20

        model = ModelController(latitude=start_lat, longitude=start_lon, depth=depth, start=start_time, step=time_step, nstep=num_steps, npart=num_particles, models=models, use_bathymetry=True, use_shoreline=True,
            time_chunk=10, horiz_chunk=4, time_method='interp')

        cache_path = os.path.join(self.cache_path, "pwsl0.nc")
        model.run("http://thredds.axiomalaska.com/thredds/dodsC/PWS_L0_FCST.nc", bathy=self.bathy_file, cache=cache_path)
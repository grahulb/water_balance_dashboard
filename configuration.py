PLUGIN_MODE = 'TEST-SUITE'	#	Possible values: 'DEBUG', 'TEST-SUITE', 'REAL'

#	Set the following for debugging or testing mode
DEBUG_BASE_FOLDER_PATH = 'D:\MTP related\Datasets\Sample'
TEST_SUITE_BASE_FOLDER_PATH = 'D:\PoCRA\deliverables\one-time_cluster_shapefiles\Prabhani'
DEBUG_OR_TEST_CROPS = ['bajri']	#	Possible values: subset of dict_crop.keys()
DEBUG_OR_TEST_RABI_CROPS = ['harbhara']
DEBUG_OR_TEST_GRADUATED_RENDERING_INTERVAL_POINTS = [0, 25, 50, 75]

#	Input-Output protocol constants
RAINFALL_CSV_FILENAME = 'Rainfall.csv'
ET0_CSV_FILENAME = 'ET0_file.csv'
POINTWISE_OUTPUT_CSV_FILENAME = 'kharif_model_pointwise_output.csv'
ALL_CUSTERS_SHAPEFILE = 'D:\PoCRA\deliverables\one-time_cluster_shapefiles\Clusters.shp'
#	Optional inputs for debugging/testing
OVERRIDE_FILECROPS_BY_DEBUG_OR_TEST_CROPS = True
CROPS_FILENAME = 'crops.csv'

#	Computation Settings
STEP = 1000.0
DEFAULT_SOWING_THRESHOLD = 30
START_DATE_INDEX = 0
MONSOON_END_DATE_INDEX = 132
END_DATE_INDEX = 364
END_DATE_INDICES = [50, 125, 200]
CADASTRAL_VULNERABILITY_DISPLAY_COLOUR_INTERVALS_COUNT = 4
CADASTRAL_VULNERABILITY_DISPLAY_COLOURS_DICT = {0: [200, 200, 255], 1: [255, 0, 255], 2: [255, 165, 0], 3: [255, 0, 0]}
from os import path

#Color Constants------------------------------------------------------------------------------------
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
YELLOW = (255, 255, 0)
DARKGREY = (40, 40, 40)
LIGHTGREY = (100, 100, 100)

#Game options---------------------------------------------------------------------------------------
FPS = 500
TITLE = "NEURAL SNAKE - FPS: "+str(FPS)
RESOLUTION = (1024, 768)
HEIGHT = RESOLUTION[1]
WIDTH = RESOLUTION[0]
BGCOLOR = BLACK

TILESIZE = 32
GRIDWIDTH = WIDTH/TILESIZE
GRIDHEIGHT = HEIGHT/TILESIZE

PLAYERVELOCITY = 1
GAME_FOLDER = path.dirname(__file__)

#GAME RUNNING OPTIONS-------------------------------------------------------------------------------

RUNTIME_NUM = 100000

#DISABLE TRAINING_MODE_ENABLED IN ORDER TO TEST OUT A SPECIFIC CHROMOSOME
TRAINING_MODE_ENABLED = True
TEST_GENOME_NUM = 0
TEST_CHROMOMSOME_NUM = 0
FITNESS_GRAPH_TYPE = "CHROMOSOME" #"GENOME"

#NEURAL NETWORK SETTINGS. (rates are percentages)
MUTATION_RATE       = 2
CROSSOVER_RATE      = 50
GENOME_SIZE         = 25
SELECTION_PERCENTAGE= 0.2
CLONING_ENABLED     = False
#activation function being used: ( tanh or sigmoid )
ACTIVATION_FUNCTION = "TANH" #"SIGMOID"
SELECTION_METHOD    = "ELITISM" #"RCHANCE"
NORMALIZE_INPUTS    = False
LEARNING_RATE       = 1.0 #range = [0, 1]

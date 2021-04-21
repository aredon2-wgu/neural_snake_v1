#Color Constants
from os import path

BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
YELLOW = (255, 255, 0)
DARKGREY = (40, 40, 40)
LIGHTGREY = (100, 100, 100)

#Game options
TITLE = "Snake - Rchance Selection"
RESOLUTION = (1024, 768)
HEIGHT = RESOLUTION[1]
WIDTH = RESOLUTION[0]
FPS = 400
BGCOLOR = BLACK

TILESIZE = 32
GRIDWIDTH = WIDTH/TILESIZE
GRIDHEIGHT = HEIGHT/TILESIZE

PLAYERVELOCITY = 1
GAME_FOLDER = path.dirname(__file__)

#disables first inital random new genome
FIRST_TIME_TRAINING = True

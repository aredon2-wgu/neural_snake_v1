#Helper Functions

import pygame as pg
from settings import *

def drawText(surface, string, size, posX, posY, color=WHITE):
    #Setup Basic Font
    basicFont = pg.font.SysFont(None, size)
    # set up the text
    text = basicFont.render(string, True, color)
    textRect = text.get_rect()
    textRect.centerx = posX
    textRect.centery = posY
    #textRect.centerx = canvas.get_rect().centerx
    #textRect.centery = canvas.get_rect().centery
    # draw the text onto the surface
    surface.blit(text, textRect)

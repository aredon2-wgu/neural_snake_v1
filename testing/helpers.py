#Helper Functions
import pygame as pg
from settings import *
import os
import numpy as np

#TODO: write all the documentation comments for all helper Functions

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

def format_time(t):

    #seconds, fractions = divmod(t, 10)

    minutes, seconds = divmod(t, 60)
    hours, minutes = divmod(minutes, 60)
    #return "%d:%d:%02d.%2d" % (hours, minutes, seconds, fractions)
    return "%d:%d:%02d" % (hours, minutes, seconds)



def sigmoid(x): #TODO change this to its appropriate name
    return 1/(1+np.exp(-x))

def tanhActivation(x):
    return np.tanh(x)

def average(num_list):
    #takes the average values of a list
    lsums = 0
    for num in num_list:
        lsums += num

    if len(num_list) == 0:
        return 0
    else:
        return lsums/len(num_list)

def standard_deviation(num_list):
    #Takes a list of numbers and finds the standard deviation.
    lst = num_list
    mean = average(lst)
    sigma = 0
    for num in lst:
        sigma += (num-mean)**2
    variance = sigma/len(lst)
    return variance**0.5

def file_is_empty(directory):
    if os.stat(directory).st_size == 0:
        return True
    else: return False

def percent_difference(numA_New, numB_Old):
    #returns the percent difference of two integer numbers.
    if numB_Old != 0:
        return (abs(numA_New-numB_Old)/numB_Old)*100
    else:
        return 0

def distanceOf(pointA, pointB):
    #returns the distance from pointA to pointB
    return ( (pointB[0] - pointA[0])**2 + (pointB[1] - pointA[1])**2 )**0.5

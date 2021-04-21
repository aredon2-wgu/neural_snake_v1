#sprites file

import pygame as pg
import random
from settings import *

class Square(pg.sprite.Sprite):
    def __init__(self, game, x, y, groups, color=WHITE):
        self.groups = groups
        pg.sprite.Sprite.__init__(self, self.groups)
        self.game = game
        self.color = color
        self.image = pg.Surface((TILESIZE, TILESIZE))
        self.image.fill(self.color)
        self.rect = self.image.get_rect()
        self.x = x
        self.y = y

    def update(self):
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

class Player(Square):
    #the game attr is the game object
    def __init__(self, game, x, y):
        self.groups = game.all_sprites_group, game.players_group
        self.game = game
        self.color = WHITE
        self.x = x
        self.y = y
        Square.__init__(self, self.game, self.x, self.y, self.groups, self.color)
        self.dx, self.dy = 1, 0
        self.vel = PLAYERVELOCITY
        self.trail = []
        self.trailSquaresCount = 0
        self.length = 0

    #dx and dy are the direction in x and y
    def move(self):
        if not self.collide_with_walls(self.dx, self.dy) and not self.collide_with_tail(self.dx, self.dy):
            if self.eat_food(self.dx, self.dy):
                #append current position to trail and add one to length
                dpos = [self.x - self.dx, self.y - self.dy]
                self.trail.append(dpos)
                self.length += 1
                self.game.food.spawnRand()

            if self.length > 0:
                #add current position to the TRAIL and pop the last trail element
                self.trail.insert(0, [self.x, self.y])
                self.trail.pop(len(self.trail)-1)

            self.x += self.dx*self.vel
            self.y += self.dy*self.vel
        else: self.die()

    #traverse the self.trail list and create a square for each element in the list
    def update_trail(self):
        trailGroups = self.game.all_sprites_group, self.game.trails_group

        #if length increases then add a new square
        if self.trailSquaresCount < self.length:
            Square(self.game, self.trail[0][0], self.trail[0][1], trailGroups, BLUE)
            self.trailSquaresCount += 1
        #if not then update the positions of the rest of the squares
        else:
            for i, square in enumerate(self.game.trails_group):
                square.x = self.trail[i][0]
                square.y = self.trail[i][1]

    def die(self):
        self.game.quit()

    def eat_food(self, dx, dy):
        for food in self.game.foods_group:
            if food.x == self.x and food.y == self.y:
                return True
        return False

    def collide_with_tail(self, dx, dy):
        for tail in self.game.trails_group:
            if tail.x == self.x and tail.y == self.y:
                return True
        return False

    def collide_with_walls(self, dx, dy):
        for wall in self.game.walls_group:
            if wall.x == self.x and wall.y == self.y:
                return True
        return False

    def update(self):
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

        self.update_trail()
        self.move()

class Wall(Square):
    def __init__(self, game, x, y):
        self.groups = game.all_sprites_group, game.walls_group
        self.game = game
        self.color = GREEN
        self.x = x
        self.y = y
        Square.__init__( self, self.game, self.x, self.y, self.groups, self.color )
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

class Food(Square):
    def __init__( self, game, x, y ):
        self.groups = game.all_sprites_group, game.foods_group
        self.game = game
        self.color = RED
        self.x = x
        self.y = y
        Square.__init__( self, self.game, self.x, self.y, self.groups, self.color )

    def spawnRand(self):
        #if the player goes over me(food) then change my x,y coordinates to somewhere else
        randX = random.randint( 1, int( ( WIDTH )/TILESIZE )-15 )+13
        randY = random.randint( 1, int( ( HEIGHT )/TILESIZE )-15 )+13
        for tail in self.game.trails_group:
            if ( tail.x == randX ) and ( tail.y == randY ):
                randX = random.randint( 1, int( ( WIDTH )/TILESIZE )-13 )+13
                randY = random.randint( 1, int( ( HEIGHT )/TILESIZE )-13 )+13
        self.x = randX
        self.y = randY

    def update(self):
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

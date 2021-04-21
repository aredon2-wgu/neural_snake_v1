#sprites file

import pygame as pg
import random#, time
from settings import *
from helpers import *

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
        self.old_dx, self.old_dy = 1, 0
        self.dx, self.dy = 1, 0
        self.vel = PLAYERVELOCITY
        self.trail = []
        self.trailSquaresCount = 0
        self.length = 0
        self.vision_inputs = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    #dx and dy are the direction in x and y
    def move(self):
        if (self.x + self.dx == self.x - self.old_dx or self.y + self.dy == self.y - self.old_dy):
            self.dx = self.old_dx
            self.dy = self.old_dy
        if not self.collide_with_walls(self.dx, self.dy) and not self.collide_with_tail(self.dx, self.dy) \
           and self.game.timeleft > 0: 
            if self.eat_food(self.dx, self.dy):
                #append current position to trail and add one to length
                dpos = [self.x - self.dx, self.y - self.dy]
                self.trail.append(dpos)
                self.length += 1
                if self.game.timeleft < 30:
                    self.game.timeleft += 10
                self.game.food.spawnRand()

            if self.length > 0:
                #add current position to the TRAIL and pop the last trail element
                self.trail.insert(0, [self.x, self.y])
                self.trail.pop(len(self.trail)-1)

            self.old_dx = self.dx
            self.old_dy = self.dy
            self.x += self.dx*self.vel
            self.y += self.dy*self.vel
            
        else:
            self.die()


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
        self.game.new()
        #self.game.tf = time.time()
        #self.game.time_survived = self.game.tf - self.game.t0
        #self.game.quit()

    def lookinDirection(self, direction):
        vision_in_direction = [0,0,0]
        self_pos = [self.x, self.y]
        position = self_pos
        distance = 0

        food_found = False
        tail_found = False

        position[0] += direction[0]
        position[1] += direction[1]
        distance +=1
        while (not (position[0] < 8 or position[1] < 1 or position[0] >= 31 or position[1] >= 23)):

            #check for food in the position
            if (not food_found and position[0] == self.game.food.x and position[1] == self.game.food.y):
                vision_in_direction[0] = distance
                food_found = True

            #check for tail at the position
            for tail in self.game.trails_group:
                if (not tail_found and position[0] == tail.x and position[1] == tail.y):
                    vision_in_direction[1] = distance
                    tail_found = True

            #look further in the direction
            position[0] += direction[0]
            position[1] += direction[1]
            distance += 1
        vision_in_direction[2] = distance

        return vision_in_direction


    def vision(self):
        inputs = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        #up
        temp = self.lookinDirection([0,-1])
        inputs[0] = temp[0]
        inputs[1] = temp[1]
        inputs[2] = temp[2]
        #upright
        temp = self.lookinDirection([1,-1])
        inputs[3] = temp[0]
        inputs[4] = temp[1]
        inputs[5] = temp[2]
        #right
        temp = self.lookinDirection([1,0])
        inputs[6] = temp[0]
        inputs[7] = temp[1]
        inputs[8] = temp[2]
        #downright
        temp = self.lookinDirection([1,1])
        inputs[9] = temp[0]
        inputs[10] = temp[1]
        inputs[11] = temp[2]
        #down
        temp = self.lookinDirection([0,1])
        inputs[12] = temp[0]
        inputs[13] = temp[1]
        inputs[14] = temp[2]
        #downleft
        temp = self.lookinDirection([-1,1])
        inputs[15] = temp[0]
        inputs[16] = temp[1]
        inputs[17] = temp[2]
        #left
        temp = self.lookinDirection([-1,0])
        inputs[18] = temp[0]
        inputs[19] = temp[1]
        inputs[20] = temp[2]
        #upleft
        temp = self.lookinDirection([-1,-1])
        inputs[21] = temp[0]
        inputs[22] = temp[1]
        inputs[23] = temp[2]

        return inputs

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
        self.vision_inputs = self.vision()

        self.update_trail()
        self.move()

class Wall(Square):
    def __init__(self, game, x, y):
        self.groups = game.all_sprites_group, game.walls_group
        self.game = game
        self.color = LIGHTGREY
        self.x = x
        self.y = y
        Square.__init__( self, self.game, self.x, self.y, self.groups, self.color )
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

class Food(Square):
    def __init__( self, game):
        self.groups = game.all_sprites_group, game.foods_group
        self.game = game
        self.color = RED
        #self.x = x
        #self.y = y
        self.spawnRand()
        Square.__init__( self, self.game, self.x, self.y, self.groups, self.color )

    def spawnRand(self):
        #if the player goes over me(food) then change my x,y coordinates to somewhere else
        randX = random.randint( 1, int( ( WIDTH )/TILESIZE )-12 )+10
        randY = random.randint( 1, int( ( HEIGHT )/TILESIZE )-4 )+2
        for tail in self.game.trails_group:
            if ( tail.x == randX ) and ( tail.y == randY ):
                randX = random.randint( 1, int( ( WIDTH )/TILESIZE )-12 )+10
                randY = random.randint( 1, int( ( HEIGHT )/TILESIZE )-12 )+10
        self.x = randX
        self.y = randY

    def update(self):
        self.rect.x = self.x * TILESIZE
        self.rect.y = self.y * TILESIZE

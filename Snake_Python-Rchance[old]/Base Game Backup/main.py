#Python 3.6.4
#pygame for 3.6
#pip -m install pygame

#TODO: print point text on the screen
#TODO: trigger game over screen which shows total points and the word game over in the player.die() method
#TODO: Make it so that you can die if you dont get points in a reasonable amount of time
#TODO: create better variable names for everything
#TODO: make it so that the highest score is shown (read from a file)

import pygame as pg
import random, sys
from os import path
from settings import *
from helpers import *
from sprites import *


class Game:
    def __init__(self):
        # initialize game window, etc
        pg.init()
        self.screen = pg.display.set_mode(RESOLUTION)
        pg.display.set_caption(TITLE)
        self.clock = pg.time.Clock()
        self.running = True
        pg.key.set_repeat(500, 100)
        self.load_data()
        self.points = 0

    def load_data(self):
        #Load map data from map.txt and store it in self.map_data
        game_folder = path.dirname(__file__)
        self.map_data = []
        with open(path.join(game_folder, "map.txt"), "rt") as f:
            for line in f:
                self.map_data.append(line)

    def new(self):
        # start a new game
        self.all_sprites_group = pg.sprite.Group()
        self.players_group = pg.sprite.Group()
        self.walls_group = pg.sprite.Group()
        self.foods_group = pg.sprite.Group()
        self.trails_group = pg.sprite.Group()
        for col, tiles in enumerate(self.map_data):
            for row, tile in enumerate(tiles):
                if tile == '1':
                    Wall(self, row, col)
                elif tile == 'P':
                    self.player = Player(self, row, col)

        self.food = Food(self, 0, 0)
        self.food.spawnRand()

    def run(self):
        # Game Loop
        self.playing = True
        while self.playing:
            self.dt = self.clock.tick(FPS) / 1000
            self.events()
            self.update()
            self.draw()

    def quit(self):
        print("Quitting Game...")
        pg.quit()
        sys.exit()

    def update(self):
        # Game Loop - Update
        self.all_sprites_group.update()
        self.points = self.player.length * 10

    def events(self):
        # Game Loop - events
        for event in pg.event.get():
            # check for closing window
            if event.type == pg.KEYDOWN:
                #Arrow Keys
                if event.key == pg.K_UP:
                    #self.player.move(dy=-1)
                    self.player.dy = -1
                    self.player.dx = 0
                elif event.key == pg.K_DOWN:
                    #self.player.move(dy=1)
                    self.player.dy = 1
                    self.player.dx = 0
                if event.key == pg.K_LEFT:
                    #self.player.move(dx=-1)
                    self.player.dx = -1
                    self.player.dy = 0
                elif event.key == pg.K_RIGHT:
                    #self.player.move(dx=1)
                    self.player.dx = 1
                    self.player.dy = 0

                if event.key == pg.K_q:
                    self.quit()

            if event.type == pg.QUIT:
                self.quit()

    def draw_grid(self):
        for x in range(0, WIDTH, TILESIZE):
            pg.draw.line(self.screen, LIGHTGREY, (x, 0), (x, HEIGHT))
        for y in range(0, HEIGHT, TILESIZE):
            pg.draw.line(self.screen, LIGHTGREY, (0, y), (WIDTH, y))

    def draw(self):
        # Game Loop - draw
        self.screen.fill(BGCOLOR)
        self.draw_grid()
        self.all_sprites_group.draw(self.screen)
        drawText(self.screen, "Points: " + str(self.points), 80, 140, 100, color=WHITE)

        # *after* drawing everything, flip the display
        pg.display.flip()

    def show_start_screen(self):
        # game splash/start screen
        pass

    def show_go_screen(self):
        # game over/continue
        pass

g = Game()
g.show_start_screen()
while g.running:
    g.new()
    g.run()
    g.show_go_screen()

pg.quit()

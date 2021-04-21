import pygame as pg
import random, sys, os, json, time
from os import path
from settings import *
from helpers import *
from sprites import *
from NeuralNet import *


class Game:
    def __init__(self):
        # initialize game window, etc
        self.t0 = time.time()
        pg.init()
        self.screen = pg.display.set_mode(RESOLUTION)
        pg.display.set_caption(TITLE)
        self.clock = pg.time.Clock()
        self.running = True
        #pg.key.set_repeat(500, 100) #Purposely disabled to not interfere with the Neural Network
        self.load_data()
        self.points = 0
        #self.t0 = 0
        #self.tf = 0
        self.time_survived = 0
        self.timeleft = 15
        self.print_runtime = 0
        #Neural Network stuff
        self.NN = Neural_network()

        self.first_time_training = FIRST_TIME_TRAINING
        if self.first_time_training:
            self.NN.initialize_ref_file()
            self.NN.new_genome()

        self.current_species = self.NN.get_current_species()
        self.current_genome = self.NN.get_current_genome()
        self.current_chromosome = self.NN.get_current_chromosome(flag=1)
        self.current_genome_PD = 0
        self.convergence = 0
        self.scores = []
        self.high_fitness_score = 0

        self.number_of_run_times = 1,000,000
        self.times_ran = 0
        self.num_of_chromosomes = 0
        self.outputs = [0,0,0,0]

        self.training_mode = True
        self.play_Test_Chromosome = 0

    def load_data(self):
        #Load map data from map.txt and store it in self.map_data
        self.map_data = []
        with open(path.join(GAME_FOLDER, "map.txt"), "rt") as f:
            for line in f:
                self.map_data.append(line)
            f.close()

    def new(self):
        # start a new game
        #self.t0 = time.time()
        if self.times_ran == self.number_of_run_times:
            print("completed run")
            self.quit()

        self.times_ran += 1

        self.timeleft = 15
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

        self.food = Food(self)

        #Neural Network stuff
        if self.training_mode == True:
            self.first_time_training = False
            #print("chromosome_"+str(self.current_chromosome))

            if self.current_chromosome < self.NN.genome_size:
                #evaluate fitness score and add it to the scores list
                #print(self.current_chromosome)

                t = round(self.time_survived, 3)
                #print(t)
                #print()
                x = self.current_chromosome
                fitness_score = round( ( t + (self.points**2)), 3)
                if fitness_score > self.high_fitness_score:
                    self.high_fitness_score = fitness_score

                #TODO::write fitness score to txt file for graphing
                write_data = str(self.num_of_chromosomes) +","+ str(fitness_score)
                self.writeFitnessData(write_data)

                self.scores.append(fitness_score)
                self.time_survived = 0
                if self.current_chromosome != self.NN.genome_size - 1:
                    self.current_chromosome += 1
                    self.num_of_chromosomes += 1
                    self.NN.update_reference_file("last_chromosome", 1)
                else:
                    self.current_chromosome = 0
                    self.NN.update_reference_file("last_chromosome", 0, flag=1)
                    #remove the first score entry from the list.
                    #self.scores.pop(0)

                    #Check for speciation and create a new species if there is
                    if self.current_genome > 0:
                        if self.NN.speciation_check() == True:
                            self.NN.new_species()
                    #write fitness scores to file
                    self.NN.new_score_file(self.scores)

                    #selection
                    selected = self.NN.selection_Rchance()

                    #crossover
                    offspring = self.NN.crossover(selected)

                    #mutation
                    mutated = self.NN.mutation(offspring)

                    #new genome
                    self.NN.new_genome(mutated)

                    #reset the score list for the next genome
                    self.scores = []

        else:
            #play_Test_Chromosome is used for manual replaying of one specific
            #chromosome when not in training mode.
            self.current_chromosome = self.play_Test_Chromosome

    def run(self):
        # Game Loop
        self.playing = True
        while self.playing:
            self.dt = self.clock.tick(FPS)/1000
            self.events()
            self.update()
            self.draw()

    def writeFitnessData(self, data):
        file_name = "Rchance_fitness_data.txt"
        if not path.isfile(file_name):
            with open(file_name, 'w') as f:
                f.write(data)
                f.close()
        else:
            with open(file_name, 'a') as f:
                f.write("\n" + data)
                f.close()

    def quit(self):
        print("Quitting Game...")
        self.tf = time.time()
        self.runTime = self.tf - self.t0
        print("RunTime: " + str(self.runTime) + " secs")
        print("RunTime: " + str(self.runTime/60) + " mins")
        pg.quit()
        sys.exit()

    def update(self):
        # Game Loop - Update
        self.all_sprites_group.update()
        self.points = self.player.length
        self.timeleft -= self.dt
        self.time_survived += self.dt
        self.print_runtime += self.dt
        #Neural Network stuff
        self.current_species = self.NN.get_current_species()
        self.current_genome = self.NN.get_current_genome()
        '''
        try:
            self.current_genome_PD = self.NN.speciation_check(flag=1)
        except:
            self.current_genome_PD = self.current_genome_PD

        try:
            self.convergence = self.NN.pop_convergence(flag=1)
        except:
            self.convergence = self.convergence
        '''


        current_chromosome_genes = self.NN.get_current_chromosome()
        self.inputs = self.player.vision_inputs
        self.outputs = self.NN.forwardProp(self.inputs, current_chromosome_genes)

    def events(self):
        # Game Loop - events
        for event in pg.event.get():
            # check for closing window
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_q:
                    self.quit()
            #Arrow Keys controlled by the Neural Network

            if event.type == pg.QUIT:
                self.quit()

        if self.outputs[0] == 1: #event.key == pg.K_UP:
            self.player.dy = -1
            self.player.dx = 0
        elif self.outputs[1] == 1: #event.key == pg.K_DOWN:
            self.player.dy = 1
            self.player.dx = 0
        elif self.outputs[2] == 1: #event.key == pg.K_LEFT:
            self.player.dx = -1
            self.player.dy = 0
        elif self.outputs[3] == 1: #event.key == pg.K_RIGHT:
            self.player.dx = 1
            self.player.dy = 0

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
        drawText(self.screen, "Points: " + str(self.points), 60, 140, 60, color=WHITE)
        drawText(self.screen, "Highest Fitness: " + str(self.high_fitness_score), 30, 140, 90, color=WHITE)
        drawText(self.screen, "Time Left: " + str(int(self.timeleft)), 60, 140, 120, color=WHITE)

        #Neural Network stats text
        drawText(self.screen, "Population Size: " + str(int(self.NN.readFile(self.NN.reference_file)["num_of_chromosomes"])), 30, 140, 150, color=WHITE)
        drawText(self.screen, "Current Species: " + str(int(self.current_species)), 30, 140, 200, color=WHITE)
        drawText(self.screen, "Current Genome: " + str(int(self.current_genome)), 30, 140, 250, color=WHITE)
        drawText(self.screen, "Current Chromosome: " + str(int(self.current_chromosome)), 30, 140, 300, color=WHITE)
        drawText(self.screen, "Mutation Rate: " + str(int(self.NN.mutation_rate))+"%", 30, 140, 325, color=WHITE)
        drawText(self.screen, "Outputs: " + str(self.outputs), 30, 140, 350, color=WHITE)
        drawText(self.screen, "Inputs:", 30, 140, 400, color=WHITE)
        drawText(self.screen, str(self.inputs[0:8]), 30, 140, 440, color=WHITE)
        drawText(self.screen, str(self.inputs[8:16]), 30, 140, 480, color=WHITE)
        drawText(self.screen, str(self.inputs[16:24]), 30, 140, 520, color=WHITE)
        drawText(self.screen, "RunTime: "+str(format_time(self.print_runtime))+" sec", 30, 140, 550, color=WHITE)
        #drawText(self.screen, "Current Genome PD: " + str(int(self.current_genome_PD)), 30, 140, 400, color=WHITE)
        #drawText(self.screen, "Convergence: " + str(int(self.convergence)), 30, 140, 450, color=WHITE)


        # *after* drawing everything, flip the display
        pg.display.flip()

    '''
    def show_start_screen(self):
        # game splash/start screen
        pass

    def show_go_screen(self):
        # game over/continue
        pass
    '''
g = Game()
#g.show_start_screen()
while g.running:
    g.new()
    g.run()
    #g.show_go_screen()

pg.quit()

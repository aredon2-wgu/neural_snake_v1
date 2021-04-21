import pygame as pg
import random, sys, os, json, time, math
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
        self.time_survived = 0
        self.timeleft = 15
        self.print_runtime = 0
        #Neural Network stuff
        self.NN = Neural_network()

        self.training_mode = TRAINING_MODE_ENABLED #disable in order to play a specific chromosome

        if self.training_mode:
            self.first_time_training = True
        else:
            #disables first inital random new genome
            self.first_time_training = False

        if self.first_time_training:
            self.NN.initialize_ref_file()
            self.NN.new_genome()

        self.current_species = self.NN.get_current_species()
        self.current_genome = self.NN.get_current_genome()
        self.current_chromosome = self.NN.get_current_chromosome(flag=1)
        self.current_genome_PD = 0
        self.convergence = 0
        self.scores = []
        self.averageTTF = 0
        self.TTFs = []
        self.high_fitness_score = 0
        self.dt = 0

        self.number_of_run_times = RUNTIME_NUM
        self.times_ran = 0
        self.num_of_chromosomes = 0
        self.outputs = [0,0,0,0]

        self.play_test_chromosome = TEST_CHROMOMSOME_NUM
        self.play_test_genome = TEST_GENOME_NUM
        self.playing = False

    def load_data(self):
        #Load map data from map.txt and store it in self.map_data
        self.map_data = []
        with open(path.join(GAME_FOLDER, "map.txt"), "rt") as f:
            for line in f:
                self.map_data.append(line)
            f.close()

    def new(self):
        #Neural Network stuff
        self.playing = False
        if self.training_mode == True:
            if self.current_chromosome < self.NN.genome_size:
                t = round(self.time_survived, 3)

                #Evaluate Fitness score
                #fitness_score = round( ( t + self.points**2), 3)
                fitness_score = 0
                if self.points == 0 and t > 14:
                    fitness_score = round((self.points+1)**2 + math.log((16-t)**2), 3)
                else:
                    fitness_score = round((self.points+1)**2 + math.log(t + 1), 3)

                if fitness_score > self.high_fitness_score:
                    self.high_fitness_score = fitness_score

                if not self.first_time_training:
                    self.scores.append(fitness_score)

                if FITNESS_GRAPH_TYPE == "CHROMOSOME":
                    write_data = str(self.times_ran) +","+ str(fitness_score)
                    self.writeFitnessData(write_data)

                #Checks conditions for the current chromosome to increment or reset
                if self.current_chromosome == 0 and self.first_time_training:
                    print("first run")
                elif self.current_chromosome < self.NN.genome_size - 1 and not self.first_time_training:
                    self.current_chromosome += 1
                    self.NN.update_reference_file("last_chromosome", 1)
                elif self.current_chromosome >= self.NN.genome_size - 1  and not self.first_time_training:
                    #write max score to graphing file

                    if FITNESS_GRAPH_TYPE == "GENOME":
                        write_data = str(self.current_genome) +","+ str(max(self.scores))
                        self.writeFitnessData(write_data)

                    self.current_chromosome = 0
                    self.NN.update_reference_file("last_chromosome", 0, flag=1)

                    #Check for speciation and create a new species if there is
                    '''
                    if self.current_genome > 0:
                        if self.NN.speciation_check() == True:
                            self.NN.new_species()
                    '''
                    #write fitness scores to file
                    self.NN.new_score_file(self.scores)

                    #GENETIC ALGORITHM CODE
                    #selection
                    selected = None
                    #if SELECTION_METHOD == "RCHANCE":
                    #selected = self.NN.selection_Rchance()
                    #elif SELECTION_METHOD == "ELITISM":
                    selected = self.NN.selection_elitism()
                    #crossover and mutation
                    cross_mutated = self.NN.crossover(selected)

                    #mutation
                    #mutated = self.NN.mutation(offspring)

                    #new genome
                    self.NN.new_genome(cross_mutated)

                    #reset the score list for the next genome
                    self.scores = []
                    #Restart the clock to avoid recording lagtime in next generation
                    self.clock = pg.time.Clock()
            self.first_time_training = False
        else:
            '''play_Test_Chromosome is used for manual replaying of one specific
            chromosome when not in training mode.'''
            self.current_chromosome = self.play_test_chromosome
            self.current_genome = self.play_test_genome

            #evaluate the score for each play through
            t = round(self.time_survived, 3)
            fitness_score = round( ( t + (self.points**2)), 3)
            #fitness_score = round( self.points, 3)
            if fitness_score > self.high_fitness_score:
                self.high_fitness_score = fitness_score

        #START A NEW GAME: Reset All variables for new game
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
        self.ControlButtons = ControlButtons(self, 3, 12)
        self.SnakeVision = SnakeVision(self, 3, 17)
        self.times_ran += 1
        self.timeleft = 15
        #self.time_survived = 0
        self.num_of_chromosomes += 1

        #check if the game has ran for the specified amount of times
        if self.times_ran > self.number_of_run_times:
            print("completed run")
            self.quit()

        self.playing = True
        self.time_survived = 0

    def run(self):
        # Game Loop
        while self.playing:
            self.dt = self.clock.tick(FPS)/1000
            multiprocessing(multithreading(self.events(), workers=5), workers=10)
            multiprocessing(multithreading(self.update(), workers=5), workers=10)
            multiprocessing(multithreading(self.draw(), workers=5), workers=10)

    def writeFitnessData(self, data):
        #Writes fitness data to file for graphing
        file_name = "Elitism_fitness_data.txt"
        if not path.isfile(file_name):
            with open(file_name, 'w') as f:
                f.write(data)
                f.close()
        else:
            with open(file_name, 'a') as f:
                f.write("\n" + data)
                f.close()

    def quit(self):
        #Completely Exit the Game and close the pygame window
        print("Quitting Game...")
        self.tf = time.time()
        self.runTime = self.tf - self.t0
        print("RunTime: " + str(self.runTime) + " secs")
        print("RunTime: " + str(self.runTime/60) + " mins")
        pg.quit()
        sys.exit()

    def update(self):
        #Update time, player length, NN stuff, and sprites
        self.all_sprites_group.update()
        self.ControlButtons.update(self.max_prob_dir)

        self.points = self.player.length
        self.timeleft -= self.dt
        self.time_survived += self.dt
        self.print_runtime += self.dt

        #Neural Network stuff
        if self.training_mode == True:
            self.current_species = self.NN.get_current_species()
            self.current_genome = self.NN.get_current_genome()
            self.current_chromosome_genes = self.NN.get_current_chromosome()
        elif self.training_mode == False:
            self.current_chromosome_genes = self.NN.get_current_chromosome(flag=2, g_idx=self.play_test_genome, last_chromosome=self.play_test_chromosome)
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

        self.inputs = self.player.vision_inputs
        multiprocessing(multithreading(self.SnakeVision.update(self.inputs), workers=3), workers=4)

        self.outputs = self.NN.forwardProp(self.inputs, self.current_chromosome_genes)

    def events(self):
        #Event Handling mainly for keypresses and player IO
        for event in pg.event.get():
            # check for closing window
            if event.type == pg.KEYDOWN:
                if event.key == pg.K_q:
                    self.quit()
            if event.type == pg.QUIT:
                self.quit()

        #look at all of the inputs see which one has the highest probability
        self.max_prob_dir = [0, 0]
        greater_num = 0
        if ACTIVATION_FUNCTION == "SIGMOID":
            greater_num = 0.5
        elif ACTIVATION_FUNCTION == "TANH":
            greater_num = 0

        for i, val in enumerate(self.outputs):
            if val > self.max_prob_dir[1] and (val > greater_num):
                self.max_prob_dir = [i, val]

        #Arrow Keys controlled by the Neural Network
        if self.max_prob_dir[0] == 0: #event.key == pg.K_UP:
            self.player.dy = -1
            self.player.dx = 0
        elif self.max_prob_dir[0] == 1: #event.key == pg.K_DOWN:
            self.player.dy = 1
            self.player.dx = 0
        elif self.max_prob_dir[0] == 2: #event.key == pg.K_LEFT:
            self.player.dx = -1
            self.player.dy = 0
        elif self.max_prob_dir[0] == 3: #event.key == pg.K_RIGHT:
            self.player.dx = 1
            self.player.dy = 0

    def draw_grid(self):
        #Draw tile grid in background
        for x in range(0, WIDTH, TILESIZE):
            pg.draw.line(self.screen, LIGHTGREY, (x, 0), (x, HEIGHT))
        for y in range(0, HEIGHT, TILESIZE):
            pg.draw.line(self.screen, LIGHTGREY, (0, y), (WIDTH, y))

    def draw(self):
        #Draw sprites and text
        self.screen.fill(BGCOLOR)
        multiprocessing(multithreading(self.draw_grid(), workers=3), workers=4)

        multiprocessing(multithreading(self.all_sprites_group.draw(self.screen), workers=3), workers=4)

        drawText(self.screen, "Points: " + str(self.points), 50, 140, 60, color=WHITE)
        drawText(self.screen, "Time Left: " + str(int(self.timeleft)), 50, 140, 90, color=WHITE)

        #Neural Network stats text
        drawText(self.screen, "Highest Fitness: " + str(self.high_fitness_score), 30, 140, 130, color=WHITE)
        drawText(self.screen, "Population Size: " + str(int(self.NN.readFile(self.NN.reference_file)["num_of_chromosomes"])), 30, 140, 160, color=WHITE)

        drawText(self.screen, "Species: " + str(int(self.current_species)), 30, 140, 190, color=WHITE)
        drawText(self.screen, "Generation: " + str(int(self.current_genome)), 30, 140, 210, color=WHITE)
        drawText(self.screen, "Snake: " + str(int(self.current_chromosome)), 30, 140, 230, color=WHITE)

        drawText(self.screen, "Mutation Rate: " + str(int(MUTATION_RATE))+"%", 30, 140, 250, color=WHITE)
        drawText(self.screen, "Crossover Rate: " + str(int(CROSSOVER_RATE))+"%", 30, 140, 270, color=WHITE)
        drawText(self.screen, "Cloning Enabled: " + str(CLONING_ENABLED), 30, 140, 290, color=WHITE)
        drawText(self.screen, "Activation: " + ACTIVATION_FUNCTION, 30, 140, 320, color=WHITE)
        drawText(self.screen, "Selection: " + SELECTION_METHOD, 30, 140, 350, color=WHITE)

        multiprocessing(multithreading(self.ControlButtons.draw(), workers=3), workers=4)
        multiprocessing(multithreading(self.SnakeVision.draw(), workers=3), workers=4)

        drawText(self.screen, "RunTime: "+str(format_time(self.print_runtime))+" sec", 30, 140, 700, color=WHITE)
        #drawText(self.screen, "Current Genome PD: " + str(int(self.current_genome_PD)), 30, 140, 400, color=WHITE)
        #drawText(self.screen, "Convergence: " + str(int(self.convergence)), 30, 140, 450, color=WHITE)

        # *after* drawing everything, flip the display
        pg.display.flip()

g = Game()
#g.show_start_screen()
multiprocessing(g.new(), workers=10)
while g.running:
    #multiprocessing(g.run(), workers=20)
    g.run()
    #g.show_go_screen()
pg.quit()

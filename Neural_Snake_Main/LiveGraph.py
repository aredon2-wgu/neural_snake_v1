import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
from settings import *
from helpers import *
import numpy as np

#style.use('fivethirtyeight')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

file_name = "Elitism_fitness_data.txt"
def animate(i):
	try:
		graph_data = open(file_name, 'r').read()
		lines = graph_data.split('\n')
		xs = []
		ys = []
		for line in lines:
			if len(line) > 1:
				x, y = line.split(',')
				xs.append(float(x))
				ys.append(float(y))
		ax.clear()
		#new_data = np.squeeze(xs, ys)
		ax.plot(xs, ys)
	except:
		print("no file to read")

	plt.ylabel("Max Fitness")
	plt.xlabel(FITNESS_GRAPH_TYPE)
	plt.title("FITNESS SCORE OVER TIME\n"+ACTIVATION_FUNCTION+" | "+SELECTION_METHOD+" | Ni:"+str(NORMALIZE_INPUTS)+" | MR:"+str(MUTATION_RATE)+"% | CR:"+\
	str(CROSSOVER_RATE)+"% | GS:"+str(GENOME_SIZE)+" | LR:"+str(LEARNING_RATE))

ani = animation.FuncAnimation(fig, animate, interval=10)

plt.show()

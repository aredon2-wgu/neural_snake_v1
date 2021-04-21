import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import numpy as np

#style.use('fivethirtyeight')
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

file_name = "Rchance_fitness_data.txt"
def animate(i):
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
    plt.ylabel("Fitness")
    plt.xlabel("Chromosome")
    plt.title("Rchance Fitness Score Over time")
ani = animation.FuncAnimation(fig, animate, interval=20)

plt.show()

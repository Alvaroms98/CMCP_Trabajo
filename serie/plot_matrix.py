import numpy as np
import matplotlib.pyplot as plt
import re

with open('output.txt') as f:
    lines = f.readlines()

tiempo = re.findall("\d+.\d+", lines[1])[0]
titulo = lines[0]

matrix = np.loadtxt('matrix_poisson_origin.txt')

plt.title(titulo)
plt.xlabel("Tiempo de cálculo: " + str(tiempo) + " segundos")
plt.imshow(matrix)
plt.savefig('plot.png')
plt.show()



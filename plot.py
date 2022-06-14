import matplotlib.pyplot as plt
import numpy as np

xs = np.arange(0, 4.25, 0.25)
ys = []
names = ['startEPerMagneses', 'endEPerMagneses', 'startOrders', 'endOrders']
name = names[1]

with open(name + '.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        ys.append(float(line.strip('\n')))

print(xs)
print(ys)

plt.scatter(xs, ys)
plt.xlabel('TkB')
plt.ylabel(name)
plt.show()
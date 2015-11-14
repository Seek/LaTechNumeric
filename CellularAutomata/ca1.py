import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def count_neighbors(nda):
    n = np.zeros(nda.shape,  dtype= int)
    n += np.roll(nda, nda.shape[0])
    n += np.roll(nda, -nda.shape[0])
    n += np.roll(nda, 1, 1)
    n += np.roll(nda, -1, 1)
    return n
width = 8
height = 8
ca = np.random.randint(0, 2, (width, height))
print(ca)
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.matshow(ca, interpolation='nearest', cmap="viridis")
print(count_neighbors(ca))
im = ax2.matshow(count_neighbors(ca), interpolation='nearest', cmap="viridis")
fig.colorbar(im, ax=ax2)
plt.show()
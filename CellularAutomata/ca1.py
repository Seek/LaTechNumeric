import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

width = 8
height = 8
ca = np.random.randint(0, 2, (width, height))
print(ca)
fig, ax = plt.subplots()
ax.imshow(ca, interpolation='nearest', cmap="Greys")
plt.show()
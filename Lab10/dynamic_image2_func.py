import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()

def f(x, y):
    return np.sin(x) + np.cos(y)

def updatefig(i):
    global x,y
    x += np.pi / 15.
    y += np.pi / 20.    
    im.set_array(f(x,y))
    plt.title(str(i))
    return im,


x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)
im = plt.imshow(f(x, y), animated=True)

ani=animation.FuncAnimation(fig, updatefig,np.arange(60), interval=50)

plt.show()

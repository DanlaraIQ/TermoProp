import matplotlib.pyplot as plt

x = [1, 2, 3, 4, 5]
y = [6, 7, 3, 6, 9]
xx = [7, 5, 6, 3, 5]
yy = [8, 2, 4, 7, 3]
plt.plot(x, y)
plt.plot(xx, yy)
plt.legend(["Grafica 1", "Grafica 2"])
plt.show()

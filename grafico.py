import matplotlib.pyplot as plt
import numpy as np

xp, yp, xs, ys = [], [], [], []

for i in range(10):     
    f = open("calidad/30-5000-" + str(i+1) + "s.txt",'r')

    xi, yi = [], []
    for row in f:
        row = row.split(' ')
        xi.append(float(row[0]))
        yi.append(int(row[1]))

    xs.append(xi)
    ys.append(yi)


for i in range(10):     
    f = open("calidad/30-5000-" + str(i+1) + "p.txt",'r')

    xi, yi = [], []
    for row in f:
        row = row.split(' ')
        xi.append(float(row[0]))
        yi.append(int(row[1]))

    xp.append(xi)
    yp.append(yi)


mean_x_axis = [i/10.0 for i in range(1, 100)]

ys_interp = [np.interp(mean_x_axis, xs[i], ys[i]) for i in range(len(xs))]
mean_y_axis = np.mean(ys_interp, axis=0)

plt.plot(mean_x_axis, mean_y_axis)

ys_interp = [np.interp(mean_x_axis, xp[i], yp[i]) for i in range(len(xp))]
mean_y_axis = np.mean(ys_interp, axis=0)

plt.plot(mean_x_axis, mean_y_axis)


plt.legend(["Secuencial", "Paralelo"], loc ="upper right")
plt.grid()
plt.xlim([0, 10])
plt.xticks(np.arange(0, 11, step=1))

plt.xlabel('Segundos', fontsize = 12)
plt.ylabel('Valor solución', fontsize = 12)
  
plt.title('Mejor solución en el tiempo: 30-5000 ', fontsize = 20)
plt.show()
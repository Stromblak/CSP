import matplotlib.pyplot as plt
import numpy as np

instancias = 5
inst = "30-5000"
instancia = "calidad/" + inst + "-"

# normal
xs, ys = [], []
for i in range(5):     
	f = open(instancia + str(i+1) + ".txt",'r')

	xi, yi = [], []
	for row in f:
		row = row.split(' ')
		xi.append(float(row[0]))
		yi.append(int(row[1]))

	xs.append(xi)
	ys.append(yi)


# CUDA
xc, yc = [], []
for i in range(5):     
	f = open(instancia + str(i+1) + "c.txt",'r')

	xi, yi = [], []
	for row in f:
		row = row.split(' ')
		xi.append(float(row[0]))
		yi.append(int(row[1]))

	xc.append(xi)
	yc.append(yi)


mean_x_axis = [i/10.0 for i in range(100)]

# secuencial
ys_interp = [np.interp(mean_x_axis, xs[i], ys[i]) for i in range(len(xs))]
mean_y_axis = np.mean(ys_interp, axis=0)
plt.plot(mean_x_axis, mean_y_axis, label= 'Secuencial')

# CUDA
ys_interp = [np.interp(mean_x_axis, xc[i], yc[i]) for i in range(len(xc))]
mean_y_axis = np.mean(ys_interp, axis=0)
plt.plot(mean_x_axis, mean_y_axis, label= 'CUDA')


plt.grid()
plt.xlim([0, 10])
plt.xticks(np.arange(0, 11, step=1))

plt.xlabel('Segundos', fontsize = 12)
plt.ylabel('Valor solución', fontsize = 12)
  
plt.suptitle('Mejor solución en el tiempo', fontsize = 20)
plt.title(inst)

plt.legend(loc ="upper right")
plt.show()
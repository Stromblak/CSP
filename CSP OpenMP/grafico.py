import matplotlib.pyplot as plt
import numpy as np

instancias = 5

instancia = "calidad/10-5000-"

# secuencial
xs, ys = [], []
for i in range(5):     
    f = open(instancia + str(i+1) + "s.txt",'r')

    xi, yi = [], []
    for row in f:
        row = row.split(' ')
        xi.append(float(row[0]))
        yi.append(int(row[1]))

    xs.append(xi)
    ys.append(yi)


# paralelo
hilos = [2, 8, 32]
xp, yp = [], []
for h in hilos:
	xi2, yi2 = [], []
	for i in range(5):
		f = open(instancia + str(i+1) + "p" + str(h) + ".txt",'r')

		xi, yi = [], []
		for row in f:
			row = row.split(' ')
			xi.append(float(row[0]))
			yi.append(int(row[1]))

		xi2.append(xi)
		yi2.append(yi)
	
	xp.append(xi2)
	yp.append(yi2)



mean_x_axis = [i/10.0 for i in range(100)]

# secuencial
ys_interp = [np.interp(mean_x_axis, xs[i], ys[i]) for i in range(len(xs))]
mean_y_axis = np.mean(ys_interp, axis=0)
plt.plot(mean_x_axis, mean_y_axis, label= 'Secuencial')


# paralelo
for j in range( len(hilos) ):
	ys_interp = [np.interp(mean_x_axis, xp[j][i], yp[j][i]) for i in range(len(xp[j]))]
	mean_y_axis = np.mean(ys_interp, axis=0)
	plt.plot(mean_x_axis, mean_y_axis, label= 'Paralelo ' + str(hilos[j]) + " hilos")


plt.grid()
plt.xlim([0, 10])
plt.xticks(np.arange(0, 11, step=1))

plt.xlabel('Segundos', fontsize = 12)
plt.ylabel('Valor solución', fontsize = 12)
  
plt.suptitle('Mejor solución en el tiempo', fontsize = 20)
plt.title('10-5000')

plt.legend(loc ="upper right")
plt.show()




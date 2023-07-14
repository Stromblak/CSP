CSP1.cpp corresponde al avance 1
CSP2.cu  corresponde al avance 2
CSP3.cpp corresponde al avance final, implementacion de OpenMP solo
CSP3.cu  corresponde al avance final, OpenMP + CUDA


Los archivos .cpp se pueden compilar usando g++ <archivo> -o <nombre> -fopenmp

El archivo CSP2.cu se compila con: nvcc CSP2.cu -o <nombre>

El archivo CSP3.cu se compila con: nvcc CSP3.cu -o <nombre> -Xcompiler "-openmp"
	- Esto me funciono a mi, sin embargo, no puedo asegurar que se compile en otras partes.


Finalmente, para ejecutar: <nombre>	

Sin embargo, hay algunas flags utiles:
	-i <direccion de instancia>: para ejecutar el algoritmo con otra instancia, por default es instancias/30-10000-1.txt
	-cuda: para usar CUDA, esto es para los archivos .cu
	-h <hilos>: para indicar el numero de hilos de OpenMP, por default es secuencial, no tiene utilidad en CSP2.cu
	-t <segundos>: por si se quiere cambiar el tiempo limite del algoritmo




PD: Esta la presentacion del avance 2 aqui, debido a que no hubo tarea para entregarla.

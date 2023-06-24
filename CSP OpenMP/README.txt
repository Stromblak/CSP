Lo importante en relacion al codigo:
	La clase Dataset contiene la funcion improve_solution(), la cual es la busqueda local implementada en el paper.
	La clase Sim contiene la funcion mutacion_busquedaLocal_evaluacion(), la cual fue paralelizada


Para compilar: g++ CSP.cpp -o CSP -fopenmp
Para ejecutar: CSP -i <direccion de la instancias> -h <numero de hilos> -t <tiempo limite: segundos>
	ej: -i instancias/30-5000-1.txt, este es la instancia por default
	-h <>: con 0 es secuencial, valor default 4
	-t <>: valor default 10


Esas son las flags mas relevantes segun yo.
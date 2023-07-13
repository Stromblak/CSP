# pragma omp parallel num_threads(hilos) if(hilos)
{
	while(tActual < tiempoMaximo){	
		// Crear antibiotico
		# pragma omp single
			crearAntibiotico();

		// Mutacion y busqueda local
		# pragma omp for
			for(int i=0; i<poblacion; i++){
				if(prob(pm)) bacterias[i].mutar();
				if(prob(bl)) bacterias[i].busquedaLocal();
				bacterias[i].actualizarFitness();
			}

		# pragma omp barrier // Sincronizacion

		// Evaluacion y clasificacion de bacterias
		# pragma omp single
		{
			evaluacion();
			clasificacion();
		}

		# pragma omp barrier // Sincronizacion

		// Elegir bacteria para actualizar feromonas
		# pragma omp single
		{
			if(!donadoras.empty()) bacteriaFeromonas = &bacterias[donadoras[rng()%donadoras.size()]];
			else bacteriaFeromonas = &bacterias[rng()%poblacion];
		}

		# pragma omp barrier // Sincronizacion

		// Actualizar feromonas
		# pragma omp for
			for(int i=0; i<dataset->m; i++){
				char base = bacteriaFeromonas->solucion[i];
				int aux = dataset->feromonas[i][base];
				dataset->feromonas[i][base] = aux*(1.0-rho) + (dataset->m - bacteriaFeromonas->fitness);
			}

		// Crear descendencia
		if(donadoras.size() > 1)
			# pragma omp for 
				for(int i=0; i<receptoras.size(); i++){
					int r1 = rng()%donadoras.size();
					int r2 = (donadoras[r1] + 1 + rng()%(donadoras.size()-1))%donadoras.size();

					bacterias[receptoras[i]].hijo(&bacterias[donadoras[r1]], &bacterias[donadoras[r2]]);
				}

		// Calcular tiempo y sincronizar
		# pragma omp single
			tActual = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count()/1000.0;
		
		// Sincronizacion
		# pragma omp barrier
	}
}
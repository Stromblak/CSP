#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <queue>
#include <chrono>
#include <string>
#include <cooperative_groups.h>
#include <omp.h>
using namespace std;
using namespace chrono;
using namespace cooperative_groups;

string instancia = "instancias/30-10000-1.txt";
int tiempoMaximo = 10, printMejorIter = 0, printNuevaMejor = 1;
double fer0 = 1000;
int seed = 140398212;
minstd_rand rng;
vector<pair<double, int>> calidadTiempo;

int usarCUDA = 0, hilos = 0, hilos2 = hilos;

int poblacion = 50;			// Poblacion
int torneos = 1;			// Numero de torneos para crear el antibiotico
double ph = 0.05;			// Probabilidad de que el individuo use 
int N = 2;					// Limite de busquedas sin mejoras para la heuristica

double pm = 0.001;			// Probabilidad mutacion
double bl = 0.2;			// Probabilidad busqueda local
double rho = 0.1;			// Factor de evaporacion de feromonas


//kernel
__global__ void kernelHamming(char *s1, char *s2, int *dist){
	extern __shared__ int temp[];
	int offset = blockIdx.x*blockDim.x;
	thread_group g = this_thread_block();
	
	int tid = g.thread_rank();

	if(s1[tid + offset] != s2[tid + offset]) temp[tid] = 1;
	else temp[tid] = 0;
	g.sync();

	for(int i=512; i>0; i/=2){
		if(tid < i && i < blockDim.x) temp[tid] += temp[tid + i];
		g.sync();
	}

	if(threadIdx.x == 0) atomicAdd(dist, temp[0]);
}

// funcion de probabilidad
bool prob(double p){
	if(rng()%1000 < p*1000) return true;
	else return false;
}


class DatasetCUDA{
	private:
		char *solDst;
		vector<int*> dist;
		int bloques, hebras, sharedBytes, n, m, id;
		vector<char*> datasetDst;
		vector<DatasetCUDA> datasetCUDA;
		vector<cudaStream_t> stream;

	public:
		DatasetCUDA(int n2, int m2, vector<string> dataset, int identificador){
			n = n2;
			m = m2;
			id = identificador;
			cudaSetDevice(id);

			datasetDst = vector<char*>(n);
			dist = vector<int*>(n);


			cudaMalloc(&solDst, m*sizeof(char));
			for(int i=0; i<n; i++){
				cudaMalloc(&datasetDst[i], m*sizeof(char));
				cudaMallocHost(&dist[i], sizeof(int));
			}

			stream = vector<cudaStream_t>(n);
			for(int i=0; i<n; i++) cudaStreamCreate(&stream[i]);

			bloques = m/1000;
			hebras = 1000;
			sharedBytes = hebras * sizeof(int);

			for(int i=0; i<n; i++) 
				cudaMemcpyAsync(datasetDst[i], &dataset[i][0], m*sizeof(char), cudaMemcpyHostToDevice, stream[i]);
		}

		vector<int> calidadCUDA(string sol){
			cudaSetDevice(id);
			cudaMemcpyAsync(solDst, &sol[0], m*sizeof(char), cudaMemcpyHostToDevice);
			for(int i=0; i<n; i++) *dist[i] = 0;

			for(int i=0; i<n; i++) kernelHamming<<<bloques, hebras, sharedBytes, stream[i]>>>(solDst, datasetDst[i], dist[i]);
			cudaDeviceSynchronize();

			vector<int> distHamming(n);
			for(int i=0; i<n; i++) distHamming[i] = *dist[i];

			return distHamming;
		}
};



// Procesamiento del dataset, algoritmo del paper, calculo distancia hamming
class Dataset{
	private:
		vector<string> dataset;
		vector<DatasetCUDA> datasetCUDA;

		void procesarDatos(){
			bases 		= vector<char>{'A', 'C', 'G', 'T'};
			feromonas 	= vector<unordered_map<char, double>>(m);
			
			for(int col=0; col<m; col++) for(char c: bases) feromonas[col][c] = fer0;

			for(int i=0; i<poblacion; i++) datasetCUDA.push_back( DatasetCUDA(n, m, dataset, i) );
		}

	public:
		int n, m;
		vector<char> bases;
		vector<unordered_map<char, double>> feromonas;

		Dataset(string instancia){	
			ifstream archivo(instancia);
			string gen;
 
			while(archivo >> gen) dataset.push_back(gen);
			archivo.close();

			n = dataset.size();
			m = dataset[0].size();
			procesarDatos();
		}

		auto heuristica(int id){
			return improve_solution(solRandom(id).first, id);
		}

		pair<string, int> improve_solution(string s, int id){
			vector<int> d, d2(n);

			if(usarCUDA) d = calidadCUDA(s, id);
			else d = calidadCPU(s);

			for(int i=0; i<n; i++) d2[i] = d[i];

			int dc = *max_element(d.begin(), d.end());
			int OldMax = dc;
			int NoImprov = 0;
			while (NoImprov <= N){
				int b = rng()%dataset.size();
			
				for(int j=0; j<m; j++){
					if( dataset[b][j] != s[j]){
						int max = -1;
						
						for(int k=0; k<n; k++){
							if(k != b){
								if(s[j] == dataset[k][j] && dataset[b][j] != dataset[k][j]) d[k] += 1;
								else if(s[j] != dataset[k][j] && dataset[b][j] == dataset[k][j]) d[k] -= 1;

								if(max < d[k]) max = d[k];
							}
						}

						if(dc >= max){
							dc = max;
							s[j] = dataset[b][j];

							for(int k=0; k<n; k++) d2[k] = d[k];
						}else for(int k=0; k<n; k++) d[k] = d2[k];
					}
				}

				if(OldMax != dc){
					OldMax = dc;
					NoImprov = 0;
				}else NoImprov += 1;
			}
			
			return pair<string, int>{s, dc};
		}

		pair<string, int> solRandom(int id){
			string sol(m, ' ');
			for(char &c: sol) c = bases[rng()%4];
			return pair<string, int>(sol, calidad(sol, id));
		}

		vector<int> calidadCPU(string sol){
			vector<int> distHamming(n);
			for(int i=0; i<n; i++) for(int j=0; j<m; j++) if( sol[j] != dataset[i][j] ) distHamming[i]++;

			return distHamming;
		}

		vector<int> calidadCUDA(string sol, int bacteria){
			return datasetCUDA[bacteria].calidadCUDA(sol);
		}

		int calidad(string sol, int bacteria){
			vector<int> distHamming;

			if(usarCUDA) distHamming = calidadCUDA(sol, bacteria);
			else distHamming = calidadCPU(sol);

			return *max_element(distHamming.begin(), distHamming.end());
		}
		
		char seguirFeromonas(int col){
			double total = 0;
			for(char base: bases) total += feromonas[col][base];

			uniform_real_distribution<> dis(0.0, total);

			double r = dis(rng);
			double suma = 0;

			char b = rng()%bases.size();
			for(char base: bases){
				suma += feromonas[col][base];
				if(r <= suma){
					b = base;
					break;
				}
			}

			return b;
		} 
};

// Clase que guarda informacion de cada individuo
class Bacteria{
	private:
		Dataset *dataset;
		bool actualizado;
		int id;

	public:
		string solucion;
		int fitness;

		Bacteria(Dataset *d, int identificador){
			dataset = d;
			id = identificador;

			pair<string, int> p;
			if(prob(ph)) p = dataset->heuristica(id);
			else p = dataset->solRandom(id);

			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void mutar(){
			for(int i=0; i<solucion.size(); i++) if(prob(pm)) solucion[i] = dataset->seguirFeromonas(i);
			actualizado = false;
		}

		void busquedaLocal(){
			auto p = dataset->improve_solution(solucion, id);

			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void actualizarFitness(){
			if(actualizado) return;

			fitness = dataset->calidad(solucion, id);
			actualizado = true;
		}

		void hijo(Bacteria *b1, Bacteria *b2){
			for(int i=0; i<solucion.size(); i++){
				if(prob(0.5)) solucion[i] = b1->solucion[i];
				else solucion[i] = b2->solucion[i];
			}

			mutar();
			actualizarFitness();
		}
};

// Simulador
class Sim{
	private:
		Dataset *dataset;
		int antibiotico, mejor, tiempoMejor;
		high_resolution_clock::time_point ti;
		vector<Bacteria> bacterias;
		vector<int> donadoras, receptoras;
		string mejorSol;

		void crearPoblacion(){
			for(int i=0; i<poblacion; i++) bacterias.push_back( Bacteria(dataset, i) );

			mejor = dataset->m+1;
			for(Bacteria &b: bacterias) if(b.fitness < mejor) mejor = b.fitness;
		}

		void crearAntibiotico(){
			antibiotico = 0;

			for(int i=0; i<torneos; i++){
				int r1 = rng()%poblacion;
				int r2 = (r1 + 1 + rng()%(poblacion-1))%poblacion;
				
				antibiotico = max(antibiotico, max(bacterias[r1].fitness, bacterias[r2].fitness));
			}
		}

		void evaluacion(){
			int mejorIter = dataset->m + 1, mejorBacteriaIter;

			for(int i=0; i<poblacion; i++){
				if(bacterias[i].fitness < mejorIter){
					mejorIter = bacterias[i].fitness;
					mejorBacteriaIter = i;
				}
			}

			if(mejorIter < mejor){
				mejor = mejorIter;
				mejorSol = bacterias[mejorBacteriaIter].solucion;

				tiempoMejor = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count()/1000.0;
				if(printNuevaMejor) cout << "Nueva solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			}

			if(printMejorIter) cout << mejorIter << " " << antibiotico << endl;
			
			float ms = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count();
			calidadTiempo.push_back( {ms/1000.0, mejor} );


		}

		void clasificacion(){
			receptoras.clear();
			donadoras.clear();
			
			for(int i=0; i<poblacion; i++){
				if(bacterias[i].fitness < antibiotico) donadoras.push_back(i);
				else receptoras.push_back(i);
			}
		}
		
		void ciclos(){
			Bacteria* bacteriaFeromonas;
			int tActual = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count()/1000.0;

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

					// Evaluacion, clasificacion de bacterias y eleccion de bacteria para actualizar feromonas
					# pragma omp single
					{
						evaluacion();
						clasificacion();
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

					// Calcular tiempo
					# pragma omp single
						tActual = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count()/1000.0;
					
					# pragma omp barrier // Sincronizacion
				}
			}
		}

	public:
		Sim(string instancia){
			ti = high_resolution_clock::now();
			dataset = new Dataset(instancia);
			mejor = dataset->m+1;
		}

		auto iniciar(){
			crearPoblacion();
			calidadTiempo.push_back( {0.0, mejor} );

			ciclos();

			fstream file;
			file.open("4.txt", ios_base::out);
		
			for(auto par: calidadTiempo) 
				file << par.first << " " << par.second << endl;

			file.close();


			cout << "Mejor solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			return pair<int, int>(mejor, tiempoMejor);
		}
};


void analisisExperimentalCUDA(){
	printNuevaMejor = 0;
	vector<string> n = {"30-", "50-"};
	vector<string> m = {"5000-", "10000-"};
	vector<int> h = {2, 4, 8};
	int instancias = 5;
	fstream file;
	rng.seed(seed);


	for(auto ni: n){
		for(auto mi: m){
			for(int i=0; i<instancias; i++){
				string instancia = "instancias/" + ni + mi + to_string(i+1) + ".txt";
				string salida = "calidad/" + ni + mi + to_string(i+1);
	
				//  CUDA
				usarCUDA = 1;
				for(int hi: h){
					hilos = hi;
					calidadTiempo = vector<pair<double, int>>();	
					rng.seed(seed);

					Sim(instancia).iniciar();
					
					file.open(salida + "-" + to_string(hi) + ".txt", ios_base::out);
					for(auto par: calidadTiempo) file << par.first << " " << par.second << endl;
					file.close();
				}


				// Sin CUDA
				usarCUDA = 0;
				hilos = 0;
				calidadTiempo = vector<pair<double, int>>();
				rng.seed(seed);

				Sim(instancia).iniciar();
				
				file.open(salida + ".txt", ios_base::out);
				for(auto par: calidadTiempo) file << par.first << " " << par.second << endl;
				file.close();
			}
		}
	}
}



int main(int argc, char *argv[]){
	rng.seed(time(NULL));
	int a = 0;

	for(int i=0; i<argc; i++){
		if( !strcmp(argv[i], "-i" ) ) instancia = argv[i+1];
		if( !strcmp(argv[i], "-t" ) ) tiempoMaximo = atof(argv[i+1]);
		if( !strcmp(argv[i], "-h" ) ){
			hilos = atoi(argv[i+1]);
			hilos2 = hilos;
		}
		if( !strcmp(argv[i], "-p" ) ) poblacion = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-tor" ) ) torneos = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-pg" ) ) ph = atof(argv[i+1]);
		if( !strcmp(argv[i], "-nb" ) ) N = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-pm" ) ) pm = atof(argv[i+1]);
		if( !strcmp(argv[i], "-bl" ) ) bl = atof(argv[i+1]);
		if( !strcmp(argv[i], "-r" ) ) rho = atof(argv[i+1]);
		if( !strcmp(argv[i], "-a" ) ) a = 1;

		if( !strcmp(argv[i], "-cuda" ) ) usarCUDA = 1;
	}

	if(a) analisisExperimentalCUDA();
	else{
		Sim s(instancia);
		s.iniciar();
	}
	
	cout << "fin" << endl;
	return 0;
}
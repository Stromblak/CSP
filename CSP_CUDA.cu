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
using namespace std;
using namespace chrono;
using namespace cooperative_groups;

string instancia = "instancias/30-5000-1.txt";
int tiempoMaximo = 10, printMejorIter = 0, printNuevaMejor = 1, hilos = 0, hilos2 = 0;
double fer0 = 1000;
int seed = 140398212;
minstd_rand rng;
vector<pair<double, int>> calidadTiempo;

int usarCUDA = 0;

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

// Procesamiento del dataset, algoritmo del paper, calculo distancia hamming
class Dataset{
	private:
		vector<string> dataset;
		char *solDst;
		vector<char*> datasetDst;
		vector<int*> dist;
		vector<cudaStream_t> stream;
		int bloques, hebras, sharedBytes;

		void procesarDatos(){
			bases 		= vector<char>{'A', 'C', 'G', 'T'};
			feromonas 	= vector<unordered_map<char, double>>(m);
			
			for(int col=0; col<m; col++) for(char c: bases) feromonas[col][c] = fer0;
		}

		void procesarCUDA(){
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
			procesarCUDA();
		}

		auto heuristica(){
			return improve_solution(solRandom().first);
		}

		pair<string, int> improve_solution(string s){
			vector<int> d, d2(n);

			if(usarCUDA) d = calidadCUDA(s);
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

		pair<string, int> solRandom(){
			string sol(m, ' ');
			for(char &c: sol) c = bases[rng()%4];
			return pair<string, int>(sol, calidad(sol));
		}

		vector<int> calidadCPU(string sol){
			vector<int> distHamming(n);
			for(int i=0; i<n; i++) for(int j=0; j<m; j++) if( sol[j] != dataset[i][j] ) distHamming[i]++;

			return distHamming;
		}

		vector<int> calidadCUDA(string sol){
			cudaMemcpyAsync(solDst, &sol[0], m*sizeof(char), cudaMemcpyHostToDevice);
			for(int i=0; i<n; i++) *dist[i] = 0;

			for(int i=0; i<n; i++) kernelHamming<<<bloques, hebras, sharedBytes, stream[i]>>>(solDst, datasetDst[i], dist[i]);
			cudaDeviceSynchronize();

			vector<int> distHamming(n);
			for(int i=0; i<n; i++) distHamming[i] = *dist[i];

			return distHamming;
		}

		int calidad(string sol){
			vector<int> distHamming;

			if(usarCUDA) distHamming = calidadCUDA(sol);
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

	public:
		string solucion;
		int fitness;

		Bacteria(Dataset *d){
			dataset = d;

			pair<string, int> p;
			if(prob(ph)) p = dataset->heuristica();
			else p = dataset->solRandom();

			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void mutar(){
			for(int i=0; i<solucion.size(); i++) if(prob(pm)) solucion[i] = dataset->seguirFeromonas(i);
			actualizado = false;
		}

		void busquedaLocal(){
			auto p = dataset->improve_solution(solucion);

			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void actualizarFitness(){
			if(actualizado) return;

			fitness = dataset->calidad(solucion);
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
			bacterias = vector<Bacteria>(poblacion, Bacteria(dataset));
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

		void clasificacion(){
			receptoras.clear();
			donadoras.clear();
			
			for(int i=0; i<poblacion; i++){
				if(bacterias[i].fitness < antibiotico) donadoras.push_back(i);
				else receptoras.push_back(i);
			}
		}

		void mutacion_busquedaLocal_evaluacion(){
			for(Bacteria &b: bacterias){
				if(prob(pm)) b.mutar();
				if(prob(bl)) b.busquedaLocal();
				b.actualizarFitness();
			}

			int b, mejorIter = dataset->m + 1;
			for(int i=0; i<poblacion; i++){
				if(bacterias[i].fitness < mejorIter){
					mejorIter = bacterias[i].fitness;
					b = i;
				}
			}

			if(mejorIter < mejor){
				mejor = mejorIter;
				mejorSol = bacterias[b].solucion;

				tiempoMejor = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count()/1000.0;
				if(printNuevaMejor) cout << "Nueva solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			}

			if(printMejorIter) cout << mejorIter << " " << antibiotico << endl;
		}

		void administrarAntibiotico(){		
			if(donadoras.size() <= 1) return;

			for(int i: receptoras){
				int r1 = rng()%donadoras.size();
				int r2 = (donadoras[r1] + 1 + rng()%(donadoras.size()-1))%donadoras.size();

				bacterias[i].hijo(&bacterias[donadoras[r1]], &bacterias[donadoras[r2]]);
			}
		}

		void actualizarFeromonas(){
			int b;
			if(!donadoras.empty()) b = donadoras[rng()%donadoras.size()];
			else b = rng()%poblacion;

			for(int i=0; i<dataset->m; i++){
				int aux = dataset->feromonas[i][bacterias[b].solucion[i]];
				dataset->feromonas[i][bacterias[b].solucion[i]] = aux*(1.0-rho) + (dataset->m - bacterias[b].fitness);
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

			while( duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count() < tiempoMaximo*1000){
				
				crearAntibiotico();
				mutacion_busquedaLocal_evaluacion();

				float ms = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count();
				calidadTiempo.push_back( {ms/1000.0, mejor} );

				clasificacion();
				actualizarFeromonas();
				administrarAntibiotico();
			}
			
			cout << "Mejor solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			return pair<int, int>(mejor, tiempoMejor);
		}
};


void analisisExperimentalCUDA(){
	printNuevaMejor = 0;
	vector<string> n = {"10-", "30-"};
	vector<string> m = {"5000-"};
	int instancias = 5;
	fstream file;

	for(auto ni: n){
		for(auto mi: m){
			for(int i=0; i<instancias; i++){
				string instancia = "instancias/" + ni + mi + to_string(i+1) + ".txt";

				//  CUDA
				usarCUDA = 1;
				calidadTiempo = vector<pair<double, int>>();	
				rng.seed(seed);

				Sim(instancia).iniciar();
				
				file.open("calidad/" + ni + mi + to_string(i+1)+ "c.txt", ios_base::out);
				for(auto par: calidadTiempo) file << par.first << " " << par.second << endl;
				file.close();

				// Sin CUDA
				usarCUDA = 0;
				calidadTiempo = vector<pair<double, int>>();
				rng.seed(seed);

				Sim(instancia).iniciar();
				
				file.open("calidad/" + ni + mi + to_string(i+1) + ".txt", ios_base::out);
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
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
#include <omp.h>
using namespace std;
using namespace chrono;

string instancia = "instancias/30-5000-1.txt";
int tiempoMaximo = 10, printMejorIter = 0, printNuevaMejor = 1, hilos = 4;
double fer0 = 1000;
minstd_rand rng;
vector<pair<float, int>> calidad;

int poblacion = 50;			// Poblacion
int torneos = 1;			// Numero de torneos para crear el antibiotico
double ph = 0.05;			// Probabilidad de que el individuo use 
int N = 2;					// Limite de busquedas sin mejoras para la heuristica

double pm = 0.001;			// Probabilidad mutacion
double bl = 0.2;			// Probabilidad busqueda local
double rho = 0.1;			// Factor de evaporacion de feromonas

int limite;


bool prob(double p){
	if(rng()%1000 < p*1000) return true;
	else return false;
}


class Dataset{
	private:
		vector<string> dataset;
		vector<unordered_map<char, vector<int>>> posiciones;

		void procesarDatos(){
			bases 		= vector<char>{'A', 'C', 'G', 'T'};
			posiciones 	= vector<unordered_map<char, vector<int>>>(m);
			feromonas 	= vector<unordered_map<char, double>>(m);
			
			for(int col=0; col<m; col++){
				for(int j=0; j<n; j++) posiciones[col][dataset[j][col]].push_back(j);
				for(char c: bases) feromonas[col][c] = fer0;
			}
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

		auto heuristica(){
			return improve_solution(solRandom().first);
		}

		pair<string, int> improve_solution(string s){
			vector<short> d(n, m), d2(n);

			for(int i=0; i<m; i++) for(int j: posiciones[i][s[i]]) d[j]--;
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

		int calidad(string sol){
			vector<short> hamming(n, m);

			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) hamming[j]--;
			return *max_element(hamming.begin(), hamming.end());
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
			# pragma omp parallel for if(hilos)
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

			while( duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count() < tiempoMaximo*1000){
				
				crearAntibiotico();
				mutacion_busquedaLocal_evaluacion();

				float ms = duration_cast<milliseconds>(high_resolution_clock::now() -  ti).count();
				calidad.push_back( { ms/1000.0, mejor} );

				clasificacion();
				actualizarFeromonas();
				administrarAntibiotico();
			}
			
			cout << "Mejor solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			return pair<int, int>(mejor, tiempoMejor);
		}
};


void analisis2(){
	printNuevaMejor = 0;
	vector<string> nm = {"10-1000-", "10-5000-", "30-1000-", "30-5000-"};
	vector<int> lim = {595, 2960, 680, 3450};

	int i = 0;
	for(int i=0; i<4; i++){
		cout << nm[i] << endl;
		limite = lim[i];

		string instancia = "instancias/" + nm[i] + to_string(i+1) + ".txt";

		auto start = high_resolution_clock::now();
		for(int i=0; i<10; i++){
			Sim s(instancia);
			s.iniciar();
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start).count();

		cout << duration/(10.0 * 1000.0) << endl;
	}
}

void analisis(){
	printNuevaMejor = 0;
	vector<string> n = {"30-"};
	vector<string> m = {"5000-"};

	for(auto ni: n){
		for(auto mi: m){
			cout << ni + mi << endl;

			hilos = 8;
			for(int i=0; i<10; i++){
				calidad = vector<pair<float, int>>();

				string instancia = "instancias/" + ni + mi + to_string(i+1) + ".txt";

				Sim s(instancia);
				s.iniciar();
				
				fstream file;
				file.open("calidad/" + ni + mi + to_string(i+1) + "p.txt", ios_base::out);
			
				for(int i=0; i<calidad.size(); i++) 
					file << to_string(calidad[i].first) + " " + to_string(calidad[i].second) << endl;

				file.close();
			}

			hilos = 0;
			for(int i=0; i<10; i++){
				calidad = vector<pair<float, int>>();

				string instancia = "instancias/" + ni + mi + to_string(i+1) + ".txt";

				Sim s(instancia);
				s.iniciar();
				
				fstream file;
				file.open("calidad/" + ni + mi + to_string(i+1) + "s.txt", ios_base::out);
			
				for(int i=0; i<calidad.size(); i++) 
					file << to_string(calidad[i].first) + " " + to_string(calidad[i].second) << endl;

				file.close();
			}
		}
	}
}


int main(int argc, char *argv[]){
	rng.seed(time(NULL));
	for(int i=0; i<argc; i++){
		if( !strcmp(argv[i], "-i" ) ) instancia = argv[i+1];
		if( !strcmp(argv[i], "-t" ) ) tiempoMaximo = atof(argv[i+1]);
		if( !strcmp(argv[i], "-h" ) ) hilos = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-p" ) ) poblacion = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-tor" ) ) torneos = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-pg" ) ) ph = atof(argv[i+1]);
		if( !strcmp(argv[i], "-nb" ) ) N = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-pm" ) ) pm = atof(argv[i+1]);
		if( !strcmp(argv[i], "-bl" ) ) bl = atof(argv[i+1]);
		if( !strcmp(argv[i], "-r" ) ) rho = atof(argv[i+1]);
	}
	omp_set_num_threads(hilos);
	
	analisis();
	//Sim s(instancia);
	//s.iniciar();

	cout << "fin" << endl;
	return 0;
}
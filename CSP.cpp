#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <queue>
using namespace std;

#include <chrono>
using namespace chrono;

double total, suma;

// Parametros
string instancia = "instancias/10-1000.txt";
int tiempoMaximo = 30;
int tuning = 0, update = 0;
minstd_rand rng;
		
double fer0 = 1000;

// Variables
double determinismo = 0.9;
int poblacion = 500;		// Poblacion
int torneos = 17;			// Numero de torneos para crear el antibiotico
double pg = 0.05;			// Probabilidad de que el individuo use 
int nb = 23;				// Limite de busquedas sin mejoras

double pm = 0.001;			// Probabilidad mutacion
double bl = 0.01;			// Probabilidad busqueda local
double rho = 0.5;			// Factor de evaporacion


bool prob(double p){
	if(rng()%1000 < p*1000) return true;
	else return false;
}


class Dataset{
	private:
		vector<int> indices;
		vector<string> dataset;
		unordered_map<char, int> calidadBase;
		vector<unordered_map<char, int>> contador;
		vector<unordered_map<char, vector<int>>> posiciones;
		double th, det;
		int th_m;

		void procesarDatos(){
			bases 		= vector<char>{'A', 'C', 'G', 'T'};
			contador 	= vector<unordered_map<char, int>>(m);
			posiciones 	= vector<unordered_map<char, vector<int>>>(m);
			feromonas 	= vector<unordered_map<char, double>>(m);

			// posiciones y contador de cada base en cada columna
			for(int col=0; col<m; col++){
				for(int j=0; j<n; j++) posiciones[col][dataset[j][col]].push_back(j);
				for(char c: bases){
					contador[col][c] = posiciones[col][c].size();
					feromonas[col][c] = fer0;
				}
			}

			// crear indice que contiene: columa con base con mayor repeticion, ..., columna con base con menor repeticion
			vector<pair<int, int>> indicesAux(m);
			for(int col=0; col<m; col++){
				int mayor = -1;
				for(auto par: contador[col]) mayor = max(par.second, mayor);
				indicesAux[col] = {mayor, col};
			}
			sort(indicesAux.begin(), indicesAux.end(), greater<pair<int, int>>());

			indices = vector<int>(m);
			for(int col=0; col<m; col++) indices[col] = indicesAux[col].second;
		}
		
		char encontrarMejorBase(int col){
			int calidadBaseMin = m+1;
			vector<char> minimos;

			for(auto par: calidadBase) calidadBaseMin = min(calidadBaseMin, par.second);
			for(auto par: calidadBase) if(par.second == calidadBaseMin) minimos.push_back(par.first);

			if(minimos.size() == 1) return minimos[0];

			int repMin = n+1;
			vector<char> minRepeticion;
			for(char c: minimos) repMin = min(repMin, contador[col][c]);
			for(char c: minimos) if(contador[col][c] == repMin) minRepeticion.push_back(c);

			return minRepeticion[ rng() % minRepeticion.size() ];
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
			det = determinismo;
			procesarDatos();
		}

		auto heuristica(){
			string mejor;
			int cal = m+1;

			for(int i=0; i<n; i++){
				int calidadNueva = calidad(dataset[i]);
				if(calidadNueva < cal ){
					mejor = dataset[i];
					cal = calidadNueva;
				}
			}

			return mejorar_solucion(mejor, nb);
		}

		pair<string, int> mejorar_solucion(string sol, short N){
			vector<short> d(n, m), d2(n, m);
			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) d[j]--;
			for(int i=0; i<m; i++) d2[i] = d[i];

			int dc = *max_element(d.begin(), d.end());
			int oldMax = dc;
			int sinMejora = 0;
			

			while (sinMejora <= N){
				int b = 0;
				
				for(int i=0; i<n; i++){
					if(d[i] == dc){
						b = i;
						break;
					}
				}

				for(int j=0; j<m; j++){
					if( dataset[b][j] != sol[j]){
						int max = -1;
						
						for(int k=0; k<n; k++){
							if(k != b){
								if(sol[j] == dataset[k][j] && dataset[b][j] != dataset[k][j]) d[k] += 1;
								else if(sol[j] != dataset[k][j] && dataset[b][j] == dataset[k][j]) d[k] -= 1;

								if(max < d[k]) max = d[k];
							}
						}

						if(dc >= max){
							dc = max;
							sol[j] = dataset[b][j];

							for(int k=0; k<n; k++) d2[k] = d[k];
						}else for(int k=0; k<n; k++) d[k] = d2[k];
					}
				}

				if(oldMax != dc){
					oldMax = dc;
					sinMejora = 0;
				}else sinMejora += 1;
			}
			
			cout << sol << endl;
			return pair<string, int>{sol, calidad(sol)};
		}

		auto solRandom(){
			string sol(m, ' ');
			for(char &c: sol) c = bases[rng()%4];
			return pair<string, int>(sol, calidad(sol));
		}

		auto busquedaLocal(string sol){
			vector<int> hamming(n, m);
			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) hamming[j]--;

			for(int col: indices){
				if( !prob(bl) ) continue;

				for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]--;

				for(char base: bases){
					calidadBase[base] = 0;

					for(int i=0; i<n; i++){
						int dif = 0;
						if(base != dataset[i][col]) dif = 1;

						if(hamming[i] + dif > calidadBase[base]) calidadBase[base] = hamming[i] + dif;
					}
				}
				sol[col] = encontrarMejorBase(col);
				
				for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]++;
			}

			int cal = *max_element(hamming.begin(), hamming.end());
			return pair<string, int>{sol, cal};
		}

		int calidad(string sol){
			vector<short> hamming(n, m);

			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) hamming[j]--;
			return *max_element(hamming.begin(), hamming.end());
		}
				
		char seguirFeromonas(int col){
			double total = 0;
			for(char base: bases) total += feromonas[col][base];

			int r = rng()%(int)total;
			double suma = 0;

			char b;
			for(char base: bases){
				suma += feromonas[col][base];
				if(r < suma){
					b = base;
					break;
				}
			}

			return b;
		} 

		char baseDif(char c){
			int i = rng()%4;
			char a = bases[i];
			if(a == c) a = bases[(i+1+rng()%3)%4];
			return a;
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
			if(rng()%100 < pg*100) p = dataset->heuristica();
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
			auto p = dataset->busquedaLocal(solucion);
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
		int antibiotico, mejor, tiempoMejor, ti;
		vector<Bacteria> bacterias;
		vector<int> donadoras, receptoras;
		string mejorSol;
		int mejorIter;

		void crearPoblacion(){
			bacterias = vector<Bacteria>(poblacion, Bacteria(dataset));
		}

		void evaluarBacterias(){
			int b;
			mejorIter = dataset->m + 1;

			for(int i=0; i<poblacion; i++){
				bacterias[i].actualizarFitness();
				if(bacterias[i].fitness < mejorIter){
					mejorIter = bacterias[i].fitness;
					b = i;
				}
			}

			if(mejorIter < mejor){
				mejor = mejorIter;
				mejorSol = bacterias[b].solucion;

				tiempoMejor = time(NULL) - ti;
				if(!tuning) cout << "Nueva solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			}

			if(!tuning && update) cout << mejorIter << " " << antibiotico << endl;
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

		void mutacion(){
			for(Bacteria &b: bacterias) if(prob(pm)) b.mutar();
		}

		void busquedaLocal(){
			for(Bacteria &b: bacterias) if(prob(bl)) b.busquedaLocal();
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

			agregarFeromonas(bacterias[b].solucion, bacterias[b].fitness);	
		}

		void agregarFeromonas(string sol, int cal){
			for(int i=0; i<sol.size(); i++) dataset->feromonas[i][sol[i]] += dataset->m - cal;
		}

	public:
		Sim(string instancia){
			ti = time(NULL);
			dataset = new Dataset(instancia);
			antibiotico = 0;
			mejor = dataset->m;
		}

		auto iniciar(){
			crearPoblacion();

			while(time(NULL) - ti <= tiempoMaximo){
				crearAntibiotico();

				mutacion();
				busquedaLocal();
				evaluarBacterias();

				clasificacion();
				actualizarFeromonas();
				administrarAntibiotico();
			}
			
			if(!tuning) cout << "Mejor solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			else cout << -mejor << endl;

			return pair<int, int>(mejor, tiempoMejor);
		}
};


int main(int argc, char *argv[]){
	rng.seed(time(NULL));
	for(int i=0; i<argc; i++){
		// Parametros
		if( !strcmp(argv[i], "-i" ) ) instancia = argv[i+1];
		if( !strcmp(argv[i], "-t" ) ) tiempoMaximo = atof(argv[i+1]);
		if( !strcmp(argv[i], "-tuning" ) ) tuning = atof(argv[i+1]);

		// Variables
		if( !strcmp(argv[i], "-d" ) ) determinismo = atof(argv[i+1]);

		if( !strcmp(argv[i], "-p" ) ) poblacion = atoi(argv[i+1]);
		if( !strcmp(argv[i], "-tor" ) ) torneos = atoi(argv[i+1]);

		if( !strcmp(argv[i], "-pg" ) ) pg = atof(argv[i+1]);
		if( !strcmp(argv[i], "-nb" ) ) nb = atoi(argv[i+1]);

		if( !strcmp(argv[i], "-pm" ) ) pm = atof(argv[i+1]);
		if( !strcmp(argv[i], "-bl" ) ) bl = atof(argv[i+1]);

		if( !strcmp(argv[i], "-r" ) ) rho = atof(argv[i+1]);
	}

	//Dataset d(instancia);
	//d.heuristica();

	Sim s(instancia);
	s.iniciar();

	cout << "fin" << endl;
	return 0;
}
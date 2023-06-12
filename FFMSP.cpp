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
string instancia = "instancias/100-800-001.txt";
double threshold = 0.80;
int tiempoMaximo = 90;
int tuning = 0, update = 0;
minstd_rand rng;

double fer0 = 1000;

// Variables
double determinismo = 0.9;
int poblacion = 500;
int torneos = 17;
double pg = 0.05;
int nb = 23;

double pm = 0.03;
double bl = 0.01;
double rho = 0.1;


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

			for(int col=0; col<m; col++){
				for(int j=0; j<n; j++) posiciones[col][dataset[j][col]].push_back(j);
				for(char c: bases){
					contador[col][c] = posiciones[col][c].size();
					feromonas[col][c] = fer0;
				}
			}

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
			int calidadBaseMax = 0;
			vector<char> maximos;
			for(auto par: calidadBase) calidadBaseMax = max(calidadBaseMax, par.second);
			for(auto par: calidadBase) if(par.second == calidadBaseMax) maximos.push_back(par.first);

			if(maximos.size() == 1) return maximos[0];

			int repMin = n;
			vector<char> minRepeticion;
			for(char c: maximos) repMin = min(repMin, contador[col][c]);
			for(char c: maximos) if(contador[col][c] == repMin) minRepeticion.push_back(c);

			if(minRepeticion.size() == 1) return minRepeticion[0];
			else return minRepeticion[ rng() % minRepeticion.size() ];
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
			th = threshold;
			det = determinismo;
			th_m = th*m;
			procesarDatos();
		}

		auto miniGrasp(int busquedas){
			vector<int> hamming(n);
			string sol(m, ' ');
			
			int columnasListas = 1;
			for(int col: indices){
				if(rng()%100 < det*100){
					for(char base: bases){
						calidadBase[base] = 0;
						for(int i=0; i<n; i++){
							int dif = 0;
							if(dataset[i][col] != base) dif = 1;
							if( hamming[i] + dif >= (int)(th*columnasListas) ) calidadBase[base]++;
						}
					}
					sol[col] = encontrarMejorBase(col);
					
				}else sol[col] = bases[rng()%4];
				
				for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]++;
				columnasListas++;
			}

			string mejorSol = sol;
			int mejorCal = calidadBase[ sol[indices.back()] ];

			// busqueda local
			while(busquedas--){

				for(int col: indices){
					for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]--;

					for(char base: bases){
						calidadBase[base] = 0;
						if(base == sol[col] && rng()%(2*m) == 0) continue;
						
						for(int i=0; i<n; i++){
							int dif = 0;
							if(base != dataset[i][col]) dif = 1;
							if(hamming[i] + dif >= (int)(th*m)) calidadBase[base]++;
						}
					}
					sol[col] = encontrarMejorBase(col);
					
					for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]++;
				}
									
				if(calidadBase[ sol[indices.back()] ] > mejorCal){
					mejorSol = sol;
					mejorCal = calidadBase[ sol[indices.back()] ];
				}
			}

			return pair<string, int>{mejorSol, mejorCal};
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
						if(hamming[i] + dif >= th_m) calidadBase[base]++;
					}
				}
				sol[col] = encontrarMejorBase(col);
				
				for(int i=0; i<n; i++) if(dataset[i][col] != sol[col]) hamming[i]++;
			}

			int cal = 0;
			for(int h: hamming) if(h >= th_m) cal++;

			return pair<string, int>{sol, cal};
		}

		int calidad(string sol){
			vector<int> hamming(n, m);
			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) hamming[j]--;

			int cal = 0;
			for(int h: hamming) if(h >= th_m) cal++;

			return cal;
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
			if(rng()%100 < pg*100) p = d->miniGrasp(nb);
			else p = d->solRandom();

			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void mutar(){
			for(int i=0; i<solucion.size(); i++){
				if(prob(pm)) solucion[i] = dataset->seguirFeromonas(i);
				//if(prob(0.005)) solucion[i] = dataset->baseDif(solucion[i]);
			}
			actualizado = false;
		}

		void busquedaLocal(){
			auto p = dataset->busquedaLocal(solucion);
			solucion = p.first;
			fitness = p.second;
			actualizado = true;
		}

		void actualizarFitness(){
			if(!actualizado){
				fitness = dataset->calidad(solucion);
				actualizado = true;
			}
		}

		void lazaro(Bacteria *b1, Bacteria *b2){
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
			for(int i=0; i<poblacion; i++) bacterias.push_back( Bacteria(dataset) );
		}

		void evaluarBacterias(){
			int b;
			mejorIter = -1;

			for(int i=0; i<poblacion; i++){
				bacterias[i].actualizarFitness();
				if(bacterias[i].fitness > mejorIter){
					mejorIter = bacterias[i].fitness;
					b = i;
				}
			}

			if(mejorIter > mejor){
				mejor = mejorIter;
				mejorSol = bacterias[b].solucion;

				tiempoMejor = time(NULL) - ti;
				if(!tuning) cout << "Nueva solucion: " << mejor << "  Tiempo: " << tiempoMejor << endl;
			}

			if(!tuning && update) cout << mejorIter << " " << antibiotico << endl;
		}

		void crearAntibiotico(){
			int antiBio = dataset->n;

			for(int i=0; i<torneos; i++){
				int r1 = rng()%poblacion;
				int r2 = (r1 + 1 + rng()%(poblacion-1))%poblacion;
				antiBio = min(antiBio, max(bacterias[r1].fitness, bacterias[r2].fitness));
			}

			antibiotico = antiBio;
		}

		void clasificacion(){
			receptoras.clear();
			donadoras.clear();
			
			for(int i=0; i<poblacion; i++){
				if(bacterias[i].fitness > antibiotico) donadoras.push_back(i);
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
			if(donadoras.size() > 1){
				for(int i: receptoras){
					int r1 = rng()%donadoras.size();
					int r2 = (donadoras[r1] + 1 + rng()%(donadoras.size()-1))%donadoras.size();
					bacterias[i].lazaro(&bacterias[donadoras[r1]], &bacterias[donadoras[r2]]);
				}
			}
		}

		void actualizarFeromonas(){

			if(!donadoras.empty()){
				int b = donadoras[rng()%donadoras.size()];
				agregarFeromonas(bacterias[b].solucion, bacterias[b].fitness);
			}else{
				int b = rng()%poblacion;
				agregarFeromonas(bacterias[b].solucion, bacterias[b].fitness);			
			}
		}

		void agregarFeromonas(string sol, int cal){
			for(int i=0; i<sol.size(); i++) dataset->feromonas[i][sol[i]] += cal;
		}

	public:
		Sim(string instancia){
			ti = time(NULL);
			dataset = new Dataset(instancia);
			antibiotico = 0;
			mejor = 0;
		}

		auto iniciar(){
			crearPoblacion();

			while(time(NULL) - ti <= tiempoMaximo && mejor != dataset->n){
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


int main2(int argc, char *argv[]){
	rng.seed(time(NULL));

	for(int i=0; i<argc; i++){
		// Parametros
		if( !strcmp(argv[i], "-i" ) ) instancia = argv[i+1];
		if( !strcmp(argv[i], "-th") ) threshold = atof(argv[i+1]);
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

	Sim s(instancia);
	s.iniciar();

	return 0;
}
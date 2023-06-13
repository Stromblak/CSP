#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <unordered_map>
#include <algorithm>
#include <random>
using namespace std;

string instancia = "instancias/10-1000.txt";
minstd_rand rng;
double p = 0.5;             // Probabilidad de seleccion aleatoria de string a comparar en la busqueda local
int N = 2000;				// Limite de busquedas sin mejoras para la heuristica
int ejecuciones = 20;		// Numero de veces que se ejecuta la heuristica


bool prob(double p){
	if(rng()%1000 < p*1000) return true;
	else return false;
}

class Dataset{
	private:
		vector<int> indices;
        vector<char> bases;
		vector<string> dataset;
		unordered_map<char, int> calidadBase;
		vector<unordered_map<char, int>> contador;
		vector<unordered_map<char, vector<int>>> posiciones;

		void procesarDatos(){
			bases 		= vector<char>{'A', 'C', 'G', 'T'};
			contador 	= vector<unordered_map<char, int>>(m);
			posiciones 	= vector<unordered_map<char, vector<int>>>(m);

			// posiciones y contador de cada base en cada columna
			for(int col=0; col<m; col++){
				for(int j=0; j<n; j++) posiciones[col][dataset[j][col]].push_back(j);
				for(char c: bases) contador[col][c] = posiciones[col][c].size();
			}
		}

	public:
		int n, m;

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
			string mejor;
			int cal = m+1;

			for(int i=0; i<n; i++){
				int calidadNueva = calidad(dataset[i]);
				if(calidadNueva < cal ){
					mejor = dataset[i];
					cal = calidadNueva;
				}
			}

			return improve_solution(mejor);
		}

		pair<string, int> improve_solution(string s){
			vector<short> d(n, m), d2(n);

			for(int i=0; i<m; i++) for(int j: posiciones[i][s[i]]) d[j]--;
			for(int i=0; i<n; i++) d2[i] = d[i];

			int dc = *max_element(d.begin(), d.end());
			int OldMax = dc;
			int NoImprov = 0;
			while (NoImprov <= N){
                int b;
				for(int i=0; i<n; i++){
					if(d[i] == dc){
						b = i;
						if(prob(p)) b = (b + 1 + rng()%(dataset.size()-1))%dataset.size();
						break;
					}
				}

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

		int calidad(string sol){
			vector<short> hamming(n, m);

			for(int i=0; i<m; i++) for(int j: posiciones[i][sol[i]]) hamming[j]--;
			return *max_element(hamming.begin(), hamming.end());
		}
};

int main(int argc, char *argv[]){
	rng.seed(time(NULL));
	for(int i=0; i<argc; i++){
		if( !strcmp(argv[i], "-i" ) ) instancia = argv[i+1];
	}
	
	Dataset d(instancia);

	int Dmin = d.m + 1;
	string s;
	#pragma omp parallel for
		for(int i=0; i<ejecuciones; i++){
			auto sol = d.heuristica();

			#pragma omp critical
				if(sol.second < Dmin){
					Dmin = sol.second;
					s = sol.first;
				}
		}

	cout << "Dmin: " << Dmin << endl;
	return 0;
}
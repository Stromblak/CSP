#include <iostream>
#include <vector>
#include <math.h>
#include <cooperative_groups.h>
#include <algorithm>
#include <random>
#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <chrono>

using namespace std;
using namespace cooperative_groups;
using namespace chrono;


__global__ void kernelHamming(char *s1, char *s2, int *dist){
	extern __shared__ int temp[];
	*dist = 0;
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

vector<string> dataset;
int n, m;
char *solDst;
vector<char*> datasetDst;
vector<int*> dist;
vector<cudaStream_t> stream;
int bloques, hebras, sharedBytes;


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

vector<int> hammingCUDA(string sol){
	cudaMemcpyAsync(solDst, &sol[0], m*sizeof(char), cudaMemcpyHostToDevice);
	for(int i=0; i<n; i++) *dist[i] = 0;

	for(int i=0; i<n; i++) kernelHamming<<<bloques, hebras, sharedBytes, stream[i]>>>(solDst, datasetDst[i], dist[i]);
	cudaDeviceSynchronize();

	vector<int> distHamming(n);
	for(int i=0; i<n; i++) distHamming[i] = *dist[i];

	return distHamming;
}

vector<int> hamming(string &sol){
	vector<int> distHamming(n);
	for(int i=0; i<n; i++) for(int j=0; j<m; j++) if( sol[j] != dataset[i][j] ) distHamming[i]++;

	return distHamming;
}


int main(){
	minstd_rand rng(time(NULL));

	n = 30;
	m = 5000;
	vector<char> bases = {'A', 'C', 'T', 'G'};

	string sol(m, ' ');
	dataset = vector<string>(n, string(m, ' '));
	for(char &c: sol) c = bases[rng()%4];
	for(string &s: dataset) for(char &c: s) c = bases[rng()%4];

	procesarCUDA();

	vector<int> cpu = hamming(sol);
	for(int i: cpu) cout << i << " ";
	cout << endl;

	vector<int> gpu = hammingCUDA(sol);
	for(int i: gpu) cout << i << " ";
	cout << endl;

	double it = 1000;
	auto start = high_resolution_clock::now();
	
	for(int i=0; i<it; i++) hamming(sol);

	auto finish = high_resolution_clock::now();
	auto f = duration_cast<milliseconds>(finish - start).count();

	cout << f/it << endl;
	
	start = high_resolution_clock::now();
	
	for(int i=0; i<it; i++) hammingCUDA(sol);

	finish = high_resolution_clock::now();
	f = duration_cast<milliseconds>(finish - start).count();

	cout << f/it << endl;

	return 0;
}





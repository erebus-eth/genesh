#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <chrono>
#include <iostream>

#include "gen_compare_gpu.h"

#include "defs.h"

using namespace std;

/**
* Generate a random base
* @return a random base ('C', 'T', 'G' or 'A')
*/ 
Base randomBase() {
	Base bases[] = { 'C', 'T', 'G', 'A' };
	return bases[rand() % 4];
}

/**
* Generate a random sequence of genes
* @param size the number of bases to add
* 
* @return the sequence generated
*/ 
GeneSequence randomGeneSequence(int size) {
	GeneSequence genes;
	genes.reserve(size);

	for (int n = 0; n < size; ++n) {
		genes.push_back(randomBase());
	}

	return genes;
}

/**
* Show matching results
* @param genes1 sequence num. 1
* @param genes2 sequence num. 2
* @param matches sequence of tuples (p1, p2, l). See findMatches() function for explanation.
*/
void showResults(GeneSequence& genes1, GeneSequence& genes2, vector<int>& matches) {
	for (int c = 0; c < matches.size(); c += 3) {
		cout << "g1: " << matches[c] << " g2: " << matches[c + 1] << " l: " << matches[c + 2] << " ";
		for (int i = 0; i < matches[c + 2]; ++i) {
			cout << (char)genes1[matches[c] + i];
		}
		cout << " ";
		for (int i = 0; i < matches[c + 2]; ++i) {
			cout << (char)genes2[matches[c + 1] + i];
		}
		cout << endl;
	}
}

/**
* Find matching genes between two sequences of a minimun given length
* @param genes1 First sequence
* @param genes2 Second sequence
* @min minimum gene length detected
* 
* @return a list of tuples (p1, p2, l) where p1 is the position of the gene in the first sequence, 
* p2 the position in the second sequence and l the length of the gene (that has to be longer than min)
*/
vector<int> findMatches(GeneSequence& genes1, GeneSequence& genes2, int min) {
	vector<int> result;
	
	int maxcg1 = genes1.size() - min;
	int maxcg2 = genes2.size() - min;
	// Iterate over both sequences
	for (int cg1 = 0; cg1 < maxcg1; ++cg1) {
		for (int cg2 = 0; cg2 < maxcg2; ++cg2) {

			// Start matching bases until reaching the end of a sequence or different bases
			int icg1 = cg1;
			int icg2 = cg2;
			while (icg1 < genes1.size() && icg2 < genes2.size() && genes1[icg1] == genes2[icg2]) {
				++icg1;
				++icg2;
			}

			// Compute gen length and add to result if is longer than min
			int gl = icg1 - cg1;
			if (gl >= min) {
				result.push_back(cg1);
				result.push_back(cg2);
				result.push_back(gl);
			}
		}
	}

	return result;
}


int main() {
	// Uncomment to generate different sequences per execution
	// srand(time(0));

	// Create random gene sequences
	GeneSequence genes1 = randomGeneSequence(10000);
	GeneSequence genes2 = randomGeneSequence(10000);

	// Compute in CPU
	auto beginCPU = std::chrono::high_resolution_clock::now();
	
	vector<int> matches = findMatches(genes1, genes2, 9);

	auto endCPU = std::chrono::high_resolution_clock::now();

	showResults(genes1, genes2, matches);
	cout << "Total secuencias encontradas: " << matches.size() / 3 << endl;

	float elapsedCPU = std::chrono::duration_cast<std::chrono::milliseconds>(endCPU - beginCPU).count() / 1000.0f;
	cout << "Tiempo CPU: " << elapsedCPU << " s." << endl;

	// Compute in GPU
	auto beginGPU = std::chrono::high_resolution_clock::now();

	matches = findMatchesGPU(genes1, genes2, 9);

	auto endGPU = std::chrono::high_resolution_clock::now();

	showResults(genes1, genes2, matches);
	cout << "Total secuencias encontradas: " << matches.size() / 3 << endl;

	float elapsedGPU = std::chrono::duration_cast<std::chrono::milliseconds>(endGPU - beginGPU).count() / 1000.0f;
	cout << "Tiempo GPU: " << elapsedGPU << " s." << endl;

	// Show speedup
	cout << "Aceleracion GPU: " << elapsedCPU / elapsedGPU << "X" << endl;

	return 0;
}
#include <cuda_runtime.h>
#include <vector>

#include "gen_compare_gpu.h"

__global__ void findMatchesKernel(int* genes1, int* genes2, int* result, int min, int maxcg1, int maxcg2) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x; // thread id
    // gen conts
    int cg1 = tid / maxcg2;
    int cg2 = tid % maxcg2; // así conseguimos que cada thread tenga un par de contadores únicos

    // no hacemos solapamientos en caso de que N_thread > N_blocks
    if (cg1 >= maxcg1 || cg2 >= maxcg2) {
        return;
    }

    // básicamente es que, una vez se sale del bucle, tenemos una coincidencia en el rango [cgX, icgX]
    int icg1 = cg1;
    int icg2 = cg2;
    while (icg1 < maxcg1 && icg2 < maxcg2 && genes1[icg1] == genes2[icg2]) {
        ++icg1;
        ++icg2;
    }

    int gl = icg1 - cg1; // gen len
    if (gl >= min) {
        int index = atomicAdd(result, 3); // no WAW
        result[index + 1] = cg1;
        result[index + 2] = cg2;
        result[index + 3] = gl;
    }
}

std::vector<int> findMatchesGPU(GeneSequence &genes1, GeneSequence &genes2, int min) {
    int maxcg1 = genes1.size() - min; // si str de len 10 y min 4, entonces el índice máx es 10 - 4 = 6 para no hacer un OOB
    int maxcg2 = genes2.size() - min;

    // cast a * para la GPU
    int* d_genes1;
    int* d_genes2;
    cudaMalloc(&d_genes1, genes1.size() * sizeof(int));
    cudaMalloc(&d_genes2, genes2.size() * sizeof(int));
    cudaMemcpy(d_genes1, genes1.data(), genes1.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_genes2, genes2.data(), genes2.size() * sizeof(int), cudaMemcpyHostToDevice);

    // para el vector [[g1, g2, len], [,,,], ...]
    int* d_result;
    cudaMalloc(&d_result, 3 * maxcg1 * maxcg2 * sizeof(int));
    cudaMemset(d_result, 0, 3 * maxcg1 * maxcg2 * sizeof(int));

    // GTX 3070 TI
    int blockSize = 256;
    int gridSize = (maxcg1 * maxcg2 + blockSize - 1) / blockSize;
    findMatchesKernel<<<gridSize, blockSize>>>(d_genes1, d_genes2, d_result, min, maxcg1, maxcg2);

    // ver el número de resultados
    int resultSize;
    cudaMemcpy(&resultSize, d_result, sizeof(int), cudaMemcpyDeviceToHost);

    // un copy de toda la vida
    std::vector<int> result(resultSize);
    cudaMemcpy(result.data(), d_result + 1, resultSize * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(d_genes1);
    cudaFree(d_genes2);
    cudaFree(d_result);

    return result;
}

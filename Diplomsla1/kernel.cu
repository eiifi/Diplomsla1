

/*
1. so vsi proteini enako dolgi ?
2. mutacija
*/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h> 
#include <random>

#define maxLenProtein 10
//#define minLenProtein 6
#define numProtein 10
#define Cr 6

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

//nikamor je namjnen zgolj prvemu clenu
enum smer{naprej, nazaj, gor, dol, levo, desno, nikamor};
enum vrsta{H, P};

struct tocka {
    int x = 0;
    int y = 0;
    int z = 0;
};

struct protein
{
    int dolzinaProteina = maxLenProtein;
    int proteinSmer[maxLenProtein];
    bool proteinVrsta[maxLenProtein];
    struct tocka tocke [maxLenProtein];
    int hevristika;
};

__global__ void addKernel(int *c, struct protein* Arr)
{
    int i = threadIdx.x;
    c[i] = Arr[i].proteinSmer[1];
}

bool enako(struct protein ena, struct protein dva) {
    bool enak = true;
    for (int i = 0; i < maxLenProtein; i++) {
        if (ena.proteinSmer[i] != dva.proteinSmer[i]) {
            return false;
        }
    }
    return true;
}

void tvori_mrezo(struct protein prot) {
    for (int i = 0; i < maxLenProtein; i++) {
        struct tocka t;
        if (i != 0) {
            t = prot.tocke[i - 1];
        }
        
        if (prot.proteinSmer[i] == nikamor) {
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == naprej) {
            t.x++;
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == nazaj) {
            t.x--;
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == gor) {
            t.z++;
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == dol) {
            t.z--;
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == levo) {
            t.y++;
        }
        if (prot.proteinSmer[i] == desno) {
            t.y--;
            prot.tocke[i] = t;
        }
    }
}

void doloci_hevristiko(struct protein prot) {

}

int main()
{
    
    //struct protein *arr = new struct protein[numProtein];
    struct protein arr1[numProtein];
    srand(time(NULL));
    for (unsigned int i = 0; i < numProtein; i++) {
        struct protein generiranProtein;
        //generiranProtein.dolzinaProteina = rand() %(maxLenProtein - minLenProtein) + minLenProtein;

        for (unsigned int j = 0; j < generiranProtein.dolzinaProteina; j++) {
            generiranProtein.proteinSmer[j] = rand() % 6;
            generiranProtein.proteinVrsta[j] = rand() % 2;
        }
        arr1[i] = generiranProtein;
    }

    arr1[0].proteinSmer[0] = 6;
 
    for (int i = 0; i < numProtein; i++) {

        // KRIZANJE
        int r1 = rand() % numProtein;

        while (r1 == i)
        {
            r1 = rand() % numProtein;
        }

        struct protein generiranProtein1;
        struct protein generiranProtein2;

        int crosoverPoint = rand() % 10;

        for (int j = 0; j < maxLenProtein; j++)
        {
            if (crosoverPoint < j )
            {
                generiranProtein1.proteinSmer[j] = arr1[i].proteinSmer[j];
                generiranProtein1.proteinVrsta[j] = arr1[i].proteinVrsta[j];
                generiranProtein2.proteinSmer[j] = arr1[r1].proteinSmer[j];
                generiranProtein2.proteinVrsta[j] = arr1[r1].proteinVrsta[j];

            }
            else {
                generiranProtein1.proteinSmer[j] = arr1[r1].proteinSmer[j];
                generiranProtein1.proteinVrsta[j] = arr1[r1].proteinVrsta[j];
                generiranProtein2.proteinSmer[j] = arr1[i].proteinSmer[j];
                generiranProtein2.proteinVrsta[j] = arr1[i].proteinVrsta[j];
            }
        }
        if (enako(generiranProtein1, arr1[i]) || enako(generiranProtein1, arr1[r1]) || rand() % 100 < 1) {
            int tockaObrta1 = rand() % maxLenProtein;
            int tockaObrta2 = rand() % maxLenProtein;
            while (tockaObrta1 == tockaObrta2) {
                tockaObrta2 = rand() % maxLenProtein;
            }
            int temp = generiranProtein1.proteinSmer[tockaObrta1];
            generiranProtein1.proteinSmer[tockaObrta1] = generiranProtein1.proteinSmer[tockaObrta2];
            generiranProtein1.proteinSmer[tockaObrta2] = temp;
        }
        if (enako(generiranProtein2, arr1[i]) || enako(generiranProtein2, arr1[r1]) || rand() % 100 < 1) {
            int tockaObrta1 = rand() % maxLenProtein;
            int tockaObrta2 = rand() % maxLenProtein;
            while (tockaObrta1 == tockaObrta2) {
                tockaObrta2 = rand() % maxLenProtein;
            }
            int temp = generiranProtein1.proteinSmer[tockaObrta1];
            generiranProtein1.proteinSmer[tockaObrta1] = generiranProtein1.proteinSmer[tockaObrta2];
            generiranProtein1.proteinSmer[tockaObrta2] = temp;
        }

        tvori_mrezo(generiranProtein1);
        tvori_mrezo(generiranProtein2);
        tvori_mrezo(arr1[i]);

        doloci_hevristiko(generiranProtein1);
        doloci_hevristiko(generiranProtein2);
        doloci_hevristiko(arr1[i]);


    }













    int c[10] = {0};

    printf(" {%d,%d,%d,%d,%d,%d,%d}\n",
        arr1[1].proteinSmer[1], arr1[2].proteinSmer[1], arr1[3].proteinSmer[1], arr1[4].proteinSmer[1], arr1[5].proteinSmer[1], arr1[6].proteinSmer[1], arr1[7].proteinSmer[1]);
    printf(" {%d,%d,%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4], c[5], c[6]);
    // Add vectors in parallel.
    int* dev_b = 0;
    struct protein* arrr = 0;
    cudaSetDevice(0);
    cudaMalloc((void**)&dev_b, 10 * sizeof(int));
    cudaMalloc((void**)&arrr, 10 * sizeof(struct protein));
    cudaMemcpy(arrr, arr1, 10 * sizeof(struct protein), cudaMemcpyHostToDevice);
    addKernel <<<1, 10 >>> (dev_b, arrr);

    cudaMemcpy(c, dev_b, 10 * sizeof(int), cudaMemcpyDeviceToHost);
    printf(" {%d,%d,%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4], c[5], c[6]);

    cudaFree(c);
    cudaFree(arrr);
    cudaFree(dev_b);

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.

    
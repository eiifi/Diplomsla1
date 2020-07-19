

/*
1. so vsi proteini enako dolgi ?
2. mutacija
*/
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h> 
#include <random>
#include <algorithm>
// dolzina proteina 5 - 256
#define maxLenProtein 256
#define minLenProtein 5
#define numProtein 210
#define Cr 6

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

//nikamor je namjnen zgolj prvemu clenu
enum smer{nikamor, naprej, nazaj, gor, dol, levo, desno};
//enum smerobratno{nikamor, nazaj, naprej, dol, gor, desno, levo};
enum vrsta{H, P};

struct tocka {
    int x = 0;
    int y = 0;
    int z = 0;
};

struct protein
{
    int dolzinaProteina = maxLenProtein;
    int proteinSmer[maxLenProtein] = {0};
    bool proteinVrsta[maxLenProtein];
    struct tocka tocke [maxLenProtein];
    int hevristika = 0;
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

void tvori_mrezo(struct protein &prot) {
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
            prot.tocke[i] = t;
        }
        if (prot.proteinSmer[i] == desno) {
            t.y--;
            prot.tocke[i] = t;
        }
    }
}

bool cikelj(struct tocka a, struct tocka b) {
    if (a.x == b.x && a.y == b.y && a.z == b.z) {
        return true;
    }
    return false;
}

bool hevristicna_povezava(struct tocka a, struct tocka b) {

    //povezav za nazaj ni potrebno preverjati, saj tocke v katerih se isce pocvezava so vsaj 3 narazen
    if (a.x == b.x && a.y == b.y && (a.z == b.z + 1 || a.z == b.z - 1)) {
        return true;
    }
    if (a.x == b.x && a.z == b.z && (a.y == b.y + 1 || a.y == b.y - 1)) {
        return true;
    }
    if (a.z == b.z && a.y == b.y && (a.x == b.x + 1 || a.x == b.x - 1)) {
        return true;
    }

    return false;
}

void doloci_hevristiko(struct protein &prot) {
    int br = 0;
    prot.hevristika = 0;
    for (int i = 0; i < maxLenProtein-3; i++) {

        for (int j = i+3; j < maxLenProtein; j++) {
            if (cikelj(prot.tocke[i], prot.tocke[j])) {
                // v primeru cikla minus dolzina proteina
                prot.hevristika -= prot.dolzinaProteina;
            }

            if (prot.proteinVrsta[i] == H && prot.proteinVrsta[j] == H) {
                if (hevristicna_povezava(prot.tocke[i], prot.tocke[j])) {
                    prot.hevristika++;
                }
            }
        }
    }
}

int main()
{
    //struct protein *arr = new struct protein[numProtein];
    int nfes = 0;
    std::cin >> nfes;
    nfes *= 1000;
    int count = 0;
    struct protein arr1[numProtein];
    srand(5);
    for (int i = 0; i < numProtein; i++) {
        struct protein generiranProtein ;
        // generiramo dolzino proteina
        generiranProtein.dolzinaProteina = rand() %(maxLenProtein - minLenProtein) + minLenProtein;
        int prev = 20;
        int premik = 0;
        //tvorimo smeri proteina pri cemer preverjamo da netvorimo kontra strani od prejsne
        for ( int j = 0; j < generiranProtein.dolzinaProteina; j++) {

            generiranProtein.proteinSmer[j] = rand() % 6+1;
            if (prev % 2 == 0) {
                premik = -1;
            }
            else {
                premik = 1;
            }
            while (prev + premik == generiranProtein.proteinSmer[j]) {
               generiranProtein.proteinSmer[j] = rand() % 6+1;
            }
            prev = generiranProtein.proteinSmer[j];


            generiranProtein.proteinVrsta[j] = rand() % 2;
        }
        arr1[i] = generiranProtein;
    }

    while (count < nfes) {

        for (int i = 0; i < numProtein; i++) {

            // KRIZANJE
            int r1 = rand() % numProtein;

            while (r1 == i)
            {
                r1 = rand() % numProtein;
            }

            struct protein generiranProtein1;
            struct protein generiranProtein2; 
            
            struct protein najdenI = arr1[i];
            struct protein najdenR1 = arr1[r1];

            int crospoint = 0;
            if (najdenI.dolzinaProteina < najdenR1.dolzinaProteina) {
                crospoint = najdenI.dolzinaProteina;
            }
            else {
                crospoint = najdenR1.dolzinaProteina;
            }
            crospoint -= 3;
            int crosoverPoint = rand() % crospoint + 3;

            
            std::copy(najdenI.proteinSmer, najdenI.proteinSmer + crosoverPoint, generiranProtein1.proteinSmer);
            std::copy(najdenR1.proteinSmer + crosoverPoint, najdenR1.proteinSmer + najdenR1.dolzinaProteina, generiranProtein1.proteinSmer + crosoverPoint);

            std::copy(najdenI.proteinVrsta, najdenI.proteinVrsta + crosoverPoint, generiranProtein1.proteinVrsta);
            std::copy(najdenR1.proteinVrsta + crosoverPoint, najdenR1.proteinVrsta + najdenR1.dolzinaProteina, generiranProtein1.proteinVrsta + crosoverPoint);

            generiranProtein1.dolzinaProteina = najdenR1.dolzinaProteina;
            

            std::copy(najdenR1.proteinSmer, najdenR1.proteinSmer + crosoverPoint, generiranProtein2.proteinSmer);
            std::copy(najdenI.proteinSmer + crosoverPoint, najdenI.proteinSmer + najdenI.dolzinaProteina, generiranProtein2.proteinSmer + crosoverPoint);

            std::copy(najdenR1.proteinVrsta, najdenR1.proteinVrsta + crosoverPoint, generiranProtein2.proteinVrsta);
            std::copy(najdenI.proteinVrsta + crosoverPoint, najdenI.proteinVrsta + najdenI.dolzinaProteina, generiranProtein2.proteinVrsta + crosoverPoint);

            generiranProtein2.dolzinaProteina = najdenI.dolzinaProteina;
            /*
            for (int j = 0; j < maxLenProtein; j++)
            {
                if (crosoverPoint < j)
                {
                    // memory copy std::copy
                    generiranProtein1.proteinSmer[j] = arr1[i].proteinSmer[j];
                    generiranProtein1.proteinVrsta[j] = arr1[i].proteinVrsta[j];
                     generiranProtein2.dolzinaProteina = arr1[i].dolzinaProteina;generiranProtein2.dolzinaProteina = arr1[i].dolzinaProteina;
                    generiranProtein2.proteinSmer[j] = arr1[r1].proteinSmer[j];
                    generiranProtein2.proteinVrsta[j] = arr1[r1].proteinVrsta[j];

                }
                else {
                    generiranProtein1.proteinSmer[j] = arr1[r1].proteinSmer[j];
                    generiranProtein1.proteinVrsta[j] = arr1[r1].proteinVrsta[j];
                    generiranProtein2.dolzinaProteina = arr1[i].dolzinaProteina;generiranProtein2.dolzinaProteina = arr1[i].dolzinaProteina;
                    generiranProtein2.proteinSmer[j] = arr1[i].proteinSmer[j];
                    generiranProtein2.proteinVrsta[j] = arr1[i].proteinVrsta[j];
                }
            }
            */
            if (enako(generiranProtein1, najdenI) || enako(generiranProtein1, najdenR1) || rand() % 100 < 2) {

                int tockaObrta1 = rand() % generiranProtein1.dolzinaProteina;
                int tockaObrta2 = rand() % generiranProtein1.dolzinaProteina;
                int tockaObrta3 = rand() % generiranProtein1.dolzinaProteina;
                int premik = 0;

                generiranProtein1.proteinSmer[tockaObrta1] = rand() % 6+1;
                if (generiranProtein1.proteinSmer[tockaObrta1 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein1.proteinSmer[tockaObrta1 - 1] + premik == generiranProtein1.proteinSmer[tockaObrta1]) {
                    generiranProtein1.proteinSmer[tockaObrta1 - 1] = rand() % 6+1;
                }

                generiranProtein1.proteinSmer[tockaObrta2] = rand() % 6+1;
                if (generiranProtein1.proteinSmer[tockaObrta2 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein1.proteinSmer[tockaObrta2 - 1] + premik == generiranProtein1.proteinSmer[tockaObrta2]) {
                    generiranProtein1.proteinSmer[tockaObrta2 - 1] = rand() % 6+1;
                }

                generiranProtein1.proteinSmer[tockaObrta3] = rand() % 6+1;
                if (generiranProtein1.proteinSmer[tockaObrta3 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein1.proteinSmer[tockaObrta3 - 1] + premik == generiranProtein1.proteinSmer[tockaObrta3]) {
                    generiranProtein1.proteinSmer[tockaObrta3 - 1] = rand() % 6+1;
                }

                generiranProtein1.proteinVrsta[tockaObrta1] = rand() % 2;
                generiranProtein1.proteinVrsta[tockaObrta2] = rand() % 2;
                generiranProtein1.proteinVrsta[tockaObrta3] = rand() % 2;

            }
            if (enako(generiranProtein2, najdenI) || enako(generiranProtein2, najdenR1) || rand() % 100 < 1) {
                int tockaObrta1 = rand() % generiranProtein2.dolzinaProteina;
                int tockaObrta2 = rand() % generiranProtein2.dolzinaProteina;
                int tockaObrta3 = rand() % generiranProtein2.dolzinaProteina;
                int premik = 0;
                generiranProtein2.proteinSmer[tockaObrta1] = rand() % 6+1;
                if (generiranProtein2.proteinSmer[tockaObrta1 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein2.proteinSmer[tockaObrta1 - 1] + premik == generiranProtein2.proteinSmer[tockaObrta1]) {
                    generiranProtein2.proteinSmer[tockaObrta1 - 1] = rand() % 6+1;
                }

                generiranProtein2.proteinSmer[tockaObrta2] = rand() % 6+1;
                if (generiranProtein2.proteinSmer[tockaObrta2 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein2.proteinSmer[tockaObrta2 - 1] + premik == generiranProtein2.proteinSmer[tockaObrta2]) {
                    generiranProtein2.proteinSmer[tockaObrta2 - 1] = rand() % 6+1;
                }

                generiranProtein2.proteinSmer[tockaObrta3] = rand() % 6+1;
                if (generiranProtein2.proteinSmer[tockaObrta3 - 1] % 2 == 0) {
                    premik = -1;
                }
                else {
                    premik = 1;
                }
                while (generiranProtein2.proteinSmer[tockaObrta3 - 1] + premik == generiranProtein2.proteinSmer[tockaObrta3]) {
                    generiranProtein2.proteinSmer[tockaObrta3 - 1] = rand() % 6+1;
                }

                generiranProtein2.proteinVrsta[tockaObrta1] = rand() % 2;
                generiranProtein2.proteinVrsta[tockaObrta2] = rand() % 2;
                generiranProtein2.proteinVrsta[tockaObrta3] = rand() % 2;

            }

            tvori_mrezo(generiranProtein1);
            tvori_mrezo(generiranProtein2);
            tvori_mrezo(najdenI);

            doloci_hevristiko(generiranProtein1);
            doloci_hevristiko(generiranProtein2);
            doloci_hevristiko(najdenI);
            
            if ( generiranProtein1.hevristika > generiranProtein2.hevristika) {
                if (najdenI.hevristika < generiranProtein1.hevristika) {
                    arr1[i] = generiranProtein1;
                }
            }
            else {
                if (najdenI.hevristika < generiranProtein2.hevristika) {
                    arr1[i] = generiranProtein2;
                }
            }

            count++;
        }
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

    
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

#include "wolffsys.h"

void ex2b();
void ex2cdf();

int main(int argc, char* argv[]){


    //ex2b();
    ex2cdf();

return 0;
}

void ex2b(){
    
    srand(time(NULL));      // Seed dependent on time
    int len = 2;
    //double temp = 0.25;
    double temp[len] = {0.25, 0.5};
    //double maxtemp = 2.0;
    double J = 1.0;
    int cycles = 1000000;
    int L = 16;

    double *mag = new double [len];
    double *amag2 = new double [len];

    /*
    for (int i=0; i<len; i++)
        temp[i] = (double) i*maxtemp/(double)len;
    */

    for (int i=0; i<len; i++){
        System_one_D sys(L, temp[i], J);

        sys.run_MC_simulation(cycles);
        sys.get_correlation(); 
        sys.get_magnetization();

        mag[i] = sys.avg_mag;
        amag2[i] = sys.avg_mag_sqrd;

    }

    std::string filename = "1Dmaggy.txt";

    std::ofstream ofile;
    ofile.open(filename);

    for (int i=0; i<len; i++){
        ofile << mag[i] << "," << amag2[i] << std::endl;
    }

    ofile.close();

    delete[] mag;
    delete[] amag2;
}   // End ex2b


void ex2cdf(){

    srand(time(NULL));      // Seed dependent on time
    int len = 100;
    double *temp = new double [len];
    double maxtemp = 2.0;
    double J = 1.0;
    int cycles = 100000;
    //int L = 16;
    //int size[2] = {5, 7};
    int size[3] = {8, 16, 32};

    double *real_avg_mag = new double [len];
    double *avg_mag_sqrd = new double [len];
    double *avg_mag_fth = new double [len];

    omp_set_num_threads(4);
    for (int i=0; i<len; i++){
        temp[i] = (double) i * maxtemp/(double)len;
        real_avg_mag[i] = 0.0;
        avg_mag_sqrd[i] = 0.0;
        avg_mag_fth[i] = 0.0;
    }


    for (int s=0; s<3; s++){
        
        int L = size[s];

        for (int i=0; i<len; i++){ 
            real_avg_mag[i] = 0.0;
            avg_mag_sqrd[i] = 0.0;
            avg_mag_fth[i] = 0.0;
        }

        int sys_size[2] = {L, L};
        int nbins = 1;
        for (int n=0; n<nbins; n++){
            #pragma omp parallel for
            for (int i=0; i<len; i++){
            
                // Initialize
                System_two_D sys(sys_size, temp[i], J);

                // Equilibrate system
                sys.equilibrate(5000);
                // Perform simulation
                sys.run_MC_simulation(cycles);

                // Store values of interest
                real_avg_mag[i] += sys.real_avg_mag;
                avg_mag_sqrd[i] += sys.avg_mag_sqrd;
                avg_mag_fth[i] += sys.avg_mag_fth;
            }
        }

        printf("Finished main loop for L = %d\n", L);

        std::string filename = "r_maggy" + std::to_string(L) + ".txt";
        std::ofstream ofile;
        ofile.open(filename);

        for (int i=0; i<len; i++){
            ofile << real_avg_mag[i]/(double)nbins << "," << avg_mag_sqrd[i]/(double)nbins
               << "," << avg_mag_fth[i]/(double)nbins << std::endl;
        }
        ofile.close();
    }

    delete[] real_avg_mag;
    delete[] avg_mag_sqrd;
    delete[] avg_mag_fth;
}

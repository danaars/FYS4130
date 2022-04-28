#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

#include "wolffsys.h"

void ex2b();
void ex2cde();

int main(int argc, char* argv[]){


    ex2b();
    //ex2cde();

return 0;
}

void ex2b(){
    
    srand(time(NULL));      // Seed dependent on time
    int len = 100;
    //double temp = 0.25;
    double *temp = new double [len];
    double maxtemp = 2.0;
    double J = 1.0;
    int cycles = 100000;
    int L = 16;

    double *mag = new double [len];
    double *amag2 = new double [len];

    for (int i=0; i<len; i++)
        temp[i] = (double) i*maxtemp/(double)len;

    for (int i=0; i<len; i++){
        System_one_D sys(L, temp[i], J);

        sys.run_MC_simulation(cycles);
        //sys.get_correlation(); 
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

    delete[] temp;
    delete[] mag;
    delete[] amag2;
}   // End ex2b


void ex2cde(){

    srand(time(NULL));      // Seed dependent on time
    int len = 100;
    double *temp = new double [len];
    double maxtemp = 2.0;
    double J = 1.0;
    int cycles = 100000;
    int L = 16;
    //int size[2] = {5, 7};
    int size[2] = {L, L};

    for (int i=0; i<len; i++){
        temp[i] = (double) i * maxtemp/(double)len;
    }

    double *real_avg_mag = new double [len];
    double *avg_mag_sqrd = new double [len];
    double *avg_mag_fth = new double [len];

    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i=0; i<len; i++){
    
        System_two_D sys(size, temp[i], J);
        //System_two_D sys(size, 1.5, J);
        //std::cout << sys.p << std::endl;
        //sys.printstate();
        //sys.equilibrate(10000);
        sys.run_MC_simulation(cycles);

        /*
        std::cout << "m: " << sys.m << std::endl;
        std::cout << "m2: " << sys.m2 << std::endl;
        std::cout << "m4: " << sys.m4 << std::endl;

        printf("Real part of average magnetization: %lf\n", sys.real_avg_mag);
        printf("Average squared magnetization: %lf\n", sys.avg_mag_sqrd);
        */
    
        /*
        real_avg_mag[i] = (sys.m/(double) cycles).real();
        avg_mag_sqrd[i] = (std::conj(sys.m) * sys.m).real() / (double) cycles;
        avg_mag_fth[i] = sys.avg_mag_fth;
        */
        real_avg_mag[i] = sys.real_avg_mag;
        avg_mag_sqrd[i] = sys.avg_mag_sqrd;
        avg_mag_fth[i] = sys.avg_mag_fth;
    }
    printf("Finished main loop\n");

    std::string filename = "maggy.txt";
    std::ofstream ofile;
    ofile.open(filename);

    for (int i=0; i<len; i++){
        ofile << real_avg_mag[i] << "," << avg_mag_sqrd[i] << std::endl;
    }
    ofile.close();

}

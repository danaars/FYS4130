#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex>
//#include <cmath>
#include <fstream>
#include <string>

class System_one_D
{
    public:

    int len;
    int *state;
    int MC_cycles;
    double beta;
    double J;
    double p;

    std::complex<double> *mag;
    std::complex<double> *prod_mag;
    double *C;

    double avg_mag, avg_mag_sqrd;

    System_one_D(int size, double temp, double coupling){
        len = size; 
        state = new int [len];
        mag = new std::complex<double> [len];
        prod_mag = new std::complex<double> [len];
        
        // Chain random initialization
        for (int i=0; i<len; i++){
            state[i] = rand() % 3;
            mag[i] = 0.0;
        }   // End for loop

        // Calculate beta and probability
        beta = 1.0/temp;
        J = coupling;
        p = calculate_prob();

        avg_mag = 0.0;
        avg_mag_sqrd = 0.0;
    }   // End System_one_D

    void printstate(){
        for (int i=0; i<len; i++){
            printf("%5d", state[i]);
        }
        printf("\n"); 
    }   // End printstate

    double calculate_prob(){
        // Probability obtained by the detailed balance equations
        p = 1.0 - std::exp(-beta * J);
        return p;
    }

    void flip(int idx, int flipto, bool firstflip = false){ 
        int prev = state[idx];
        double r;

        if (firstflip == true){
            r = (double) rand()/RAND_MAX;
            if (state[idx] == 0){
                state[idx] = (r > 0.5) ? 1:2;
            }
            else if (state[idx] == 1){
                state[idx] = (r > 0.5) ? 0:2;
            }
            else if (state[idx] == 2){
                state[idx] = (r > 0.5) ? 0:1;
            }

            flipto = state[idx];
        }
        else {
            state[idx] = flipto;
        }

        int next_idx;
        // Periodic boundary condition
        // Right
        if (prev == state[(len + idx + 1) % len]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_idx = (len + idx + 1) % len;
                //printf("Going to index: %d\n", next_idx);
                flip(next_idx, flipto);
            }
        }
        // Left
        if (prev == state[(len + idx - 1) % len]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_idx = (len + idx - 1) % len;
                //printf("Going to index: %d\n", next_idx);
                flip(next_idx, flipto);
            }
        } 
    }   // End flip

    void store_magnetization(){
        // calculate and store the magnetization
        const std::complex<double> imag(0.0, 2*M_PI/3);  // i*2*pi/3
        std::complex<double> arg, m;
        std::complex<double>  m0 = std::exp(imag * (double) state[0]);

        for (int i=0; i<len; i++){

            arg = imag * (double) state[i]; 
            m = std::exp(arg);

            mag[i] += m;
            prod_mag[i] += std::conj(m0) * m;
        } 
    } 

    /*
    void average_magnetization(){
        // Calculate the average magnetization for the sim
        
        std::complex<double> cycles = std::complex<double> (MC_cycles, 0);
        for (int i=0; i<len; i++){
            mag[i] = mag[i] / cycles;
            //real_mag[i] = (double) real(mag[i]);
            //abs_mag_sqrd[i] = abs_mag_sqrd[i] / MC_cycles;
        }
    }   // End average_magnetization
    */

    void run_MC_simulation(int cycles){
    
        MC_cycles = cycles;
        int idx;
        int flipto;
        bool firstflip = true;
        for (int c=0; c<MC_cycles; c++){
            idx = rand() % len;     // Choose random index
            flip(idx, flipto, firstflip);// Perform one cluster flip with the Wolff algo
            store_magnetization();  // Store values for statistics
            //printstate();
        }
    }   // End run_MC_simulation

    void get_correlation(){
        // Correlation between site 0 and site r
        C = new double [len];

        double mc = (double) MC_cycles;

        for (int i=0; i<len; i++){
            C[i] = (double) (prod_mag[i]/mc - std::conj(mag[0])/mc * mag[i]/mc).real();
        }

        // Write C(r) to file
        std::string filename = "corr1D_" + std::to_string(1/beta) + ".txt";

        std::ofstream ofile;
        ofile.open(filename);

        for (int i=0; i<len; i++){
            ofile << C[i] << std::endl;
        }

        ofile.close();
        
    }   // End get_correlation

    void get_magnetization(){

        double mc = (double) MC_cycles;
        double m = 0.0;
        double m2 = 0.0;

        //std::string filename = "1Dmaggy.txt";

        //std::ofstream ofile;
        //ofile.open(filename); 

        for (int i=0; i<len; i++){
            m += (mag[i]/mc).real();
            m2 += (std::conj(mag[i])*mag[i]/mc).real();
        }

        avg_mag = m;
        avg_mag_sqrd = m2;
    
    }   // End get_magnetization
    ~System_one_D(){
        delete[] state;
        delete[] mag;
        delete[] prod_mag;
        delete[] C;
    }   // End ~System_one_D
};



class System_two_D
{
    public:

    //int dims;
    int *size;
    int rows, cols, MC_cycles;
    int **state;
    double beta, J, p;
    double real_avg_mag, avg_mag_sqrd, avg_mag_fth;
    std::complex<double> m;
    double m2, m4;

    //std::complex<double> **mag;
    std::complex<double> *LUTmag;
    //System_two_D(int dims, int *size){

    System_two_D(int *size, double temp, double coupling){

        rows = *(size);
        cols = *(size + 1);

        // Grid memory allocation
        state = new int*[rows];
        //mag = new std::complex<double>*[rows];

        for (int i=0; i<rows; i++){
            state[i] = new int[cols];
            //mag[i] = new std::complex<double>[rows];
        }

        // Grid random spin initialization
        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                state[i][j] = rand() % 3;       // Rand int between 0 and 3
                //mag[i][j] = 0.0;
            }
        } // End for

        // Calculate beta and probability
        beta = 1.0/temp;
        J = coupling;
        p = calculate_prob();

        m = std::complex<double>(0.0, 0.0);
        m2 = 0.0;
        m4 = 0.0;

        real_avg_mag = 0.0;
        avg_mag_sqrd = 0.0;
        avg_mag_fth = 0.0;

        const std::complex<double> imag(0.0, 2*M_PI/3);  // i*2*pi/3

        LUTmag = new std::complex<double> [3];
        for (int i=0; i<3; i++){
            LUTmag[i] = std::exp(imag * (double)i);
        }

    }   // End System class declarator

    void printstate(){

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){
                printf("%5d", state[i][j]);
            }
            printf("\n");
        }
    }   // End printstate

    double calculate_prob(){
        // Probability obtained by the detailed balance equations
        return 1.0 - std::exp(-beta * J);
    }   // End calculate_prob

    void flip(int x, int y, int flipto, bool firstflip = false){
    
        int prev = state[x][y];
        double r;

        if (firstflip == true){
            r = (double) rand()/RAND_MAX;
            // Change spin at position x, y
            if (state[x][y] == 0){
                //state[x][y] = (r > 0.5) ? 1:2;    
                state[x][y] = 1;    
                //printf("flip\n");
            }
            else if (state[x][y] == 1){
                //state[x][y] = (r > 0.5) ? 0:2;
                state[x][y] = 2;    
                //printf("flip\n");
            }
            else if (state[x][y] == 2){
                //state[x][y] = (r > 0.5) ? 0:1;
                state[x][y] = 0;    
                //printf("flip\n");
            }

            flipto = state[x][y];
        }
        else {
            state[x][y] = flipto;
        }

        //printf("Changed spin at (%d, %d) from %d -> %d\n", x, y, prev, state[x][y]);
        // Check neighbouring spins, move if favourable

        int next_x, next_y;
        // Periodic boundary conditions
        // Right
        if (prev == state[(rows + x + 1) % rows][y]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_x = (rows + x + 1) % rows;
                //printf("Going to (%d, %d)\n", next_x, y);
                flip(next_x, y, flipto);
                //printstate();
            }
        }
        // Down
        if (prev == state[x][(cols + y + 1) % cols]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_y = (cols + y + 1) % cols;
                //printf("Going to (%d, %d)\n", x, next_y);
                flip(x, next_y, flipto);
                //printstate();
            }
        }
        // Left
        if (prev == state[(rows + x - 1) % rows][y]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_x = (rows + x - 1) % rows;
                //printf("Going to (%d, %d)\n", next_x, y);
                flip(next_x, y, flipto);
                //printstate();
            }
        }
        // Up
        if (prev == state[x][(cols + y - 1) % cols]){
            r = (double) rand()/RAND_MAX;
            if (r < p){
                next_y = (cols + y - 1) % cols;
                //printf("Going to (%d, %d)\n", x, next_y);
                flip(x, next_y, flipto);
                //printstate();
            }
        }

    }   // End flip

    void store_magnetization(){

        std::complex<double> tmp_m(0.0, 0.0);
        double tmp_m2;

        for (int i=0; i<rows; i++){
            for (int j=0; j<cols; j++){

                tmp_m += (LUTmag[ state[i][j] ] / ((double) (rows * cols)));
                //tmp_m2 = (std::conj(tmp_m) * tmp_m).real();

                //m = std::exp(arg);

                //mag[i][j] += m;
            }
        }
        tmp_m2 = abs(tmp_m) * abs(tmp_m);

        m += tmp_m;
        m2 += tmp_m2;
        m4 += tmp_m2 * tmp_m2;
    }   // End store_magnetization

    void run_MC_simulation(int cycles){

        //std::cout << "In run mc:" << std::endl;
        //std::cout << m << std::endl;
        
        MC_cycles = cycles;
        int flipto, init_x, init_y;
        bool firstflip = true;
        
        for (int c=0; c<MC_cycles; c++){
            //std::cout << "In main loop:" << std::endl;
            //std::cout << m << std::endl;
            init_x = rand() % rows;
            init_y = rand() % cols;
            
            //printf("MC cycle nr: %d, initial guess (%d, %d)\n", c, init_x, init_y);
            //printstate();
            flip(init_x, init_y, flipto, firstflip);
            //printstate();
            store_magnetization();
            //std::cout << "iter: " << c+1<< ", m = " << m << std::endl;
            //printf("\n");
        }

        //printstate();
        get_magnetizations();
//        std::cout << "Inside run_MC:" << std::endl;
//        std::cout << "real average magnetization: " << real_avg_mag << std::endl;
//        std::cout << "absolute magnetization squared: " << avg_mag_sqrd << std::endl;
        //printf("\n");
    }   // End run_MC_simulation

    void equilibrate(int cycles){
        int flipto, init_x, init_y;
        bool firstflip = true;

        for (int i=0; i<cycles; i++){ 
            init_x = rand() % rows;
            init_y = rand() % cols;
            flip(init_x, init_y, flipto, firstflip);
        }

    }
    
    void get_magnetizations(){
        //double L_sqrd = (double) (rows * cols); 
        //std::complex<double> m = 0.0;
        //double m2 = 0.0;

        //for (int i=0; i<rows; i++){
            //for (int j=0; j<cols; j++){
                //m += mag[i][j]/(double) MC_cycles;
                //m2 += (std::conj(mag[i][j]) * mag[i][j]).real() / (double) MC_cycles;
            //}
        //}
        double s = (double) MC_cycles;
        //printf("double MC_cycles = %lf\n", s);
        
        //std::cout << m/s << std::endl;
        //std::cout << (m/s).real() << std::endl;
        real_avg_mag = (m/s).real();
        avg_mag_sqrd = (m2/s);
        avg_mag_fth = (m4/s);

    }   // End get_magnetization

    ~System_two_D(){
        // Deallocate heap memory
        for (int i=0; i<rows; i++){
            delete[] state[i];
            //delete[] mag[i];
        }
        delete[] state;
        delete[] LUTmag;
        //delete[] mag;
    }   // End ~System_two_D
};

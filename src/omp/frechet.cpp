#include <iostream>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <omp.h>

using namespace std;

double* lecture (const char* filename, int* n)
{
    ifstream f(filename);
    int size;
    f>>size;
    (*n) = size;
    double* t = new double[2*size];
    string line;
    int ptr = 0;
    f>>line;
    while (!f.eof()) {
        int pos = line.find(";");
        t[ptr] = stod(line.substr(0, pos));
        ptr++;
        t[ptr] = stod(line.substr(pos + 1));
        ptr++;
        f>>line;
    }
    f.close();
    return t;
}

double dist_e(double* v1, double* v2) {
    return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]));
}

double* matrice_distance(double* traj1, int n1, double* traj2, int n2) {
    double* distances = new double[n1*n2];
    for (int i=0; i<n1; i++){
        for (int j=0; j<n2; j++) distances[i*n2+j] = dist_e(traj1+2*i, traj2+2*j);
    }
    return distances;
}

double* matrice_distanceOpti(double* traj1, int n1, double* traj2, int n2) {

    double* distances = new double[n1*n2];
    int k1, k2, k;
    
    double dMax;

    k = min( n1, n2);
 
    #pragma omp parallel for num_threads(4)
    for (int i=0; i < n1*n2; i++) {
        distances[i] = -1;
    }
   
   
    // La diagonale
    dMax = 0.0f;
    
    #pragma omp parallel for reduction(max:dMax) num_threads(4)
    for (int i=0; i<k; i++){
        int ii = 2 * i;
        int diag = i*n2 + i;
        distances[diag] = dist_e(traj1 + ii, traj2 + ii);
        dMax = fmax( distances[diag], dMax);
    }

    if( n1 < n2 ) {
        k1 = 2*(n1-1);
        k2 = k1*n2;

	// derniere ligne
        #pragma omp parallel for reduction(max:dMax) num_threads(4)
        for (int j=n1; j<n2; j++){
            int diag = k2 + j;
            distances[diag] = dist_e(traj1 + k1, traj2 + 2 * j);
            dMax = fmax( distances[diag], dMax);
        }
    } else {
        k2 = 2*(n2 - 1);
	
	// derniere colonne
        #pragma omp parallel for reduction(max:dMax) num_threads(4)
        for (int j=n2; j < n1; j++){
            int diag = j * n2 + n2 - 1;
            distances[diag] = dist_e(traj1 + 2 * j, traj2 + k2);
            dMax = fmax( distances[diag], dMax);
        }
    }

    // Les lignes
    k2 = k - 1;

     #pragma omp parallel for num_threads(4)
    for (int i = 0; i < n1; i++) {
    	int diag = i * n2;
    	int i2 = 2 * i;
        for (int j = i+1; j < n2; j++) {
            int ligne = diag + j;
            distances[ligne] = dist_e(traj1 + i2, traj2 + 2 * j);
            if (distances[ligne] > dMax) {
            	distances[ligne] = -1;
            	break;
            }
        }
    }
    
    // Les colonnes
    k2 = k - 1;

     #pragma omp parallel for num_threads(4)
    for (int j = 0; j < n2; j++) {
        for (int i = j+1; i < n1; i++) {
            int diag = i * n2 + j;
            distances[diag] = dist_e(traj1 + 2 * i, traj2 + 2 * j);
            if (distances[diag] > dMax) {
            	distances[diag] = 0;
            	break;
            }
        }
    }
    
    cout << "dMax = " << dMax << endl;

        return distances;
    }


    double min(double a, double b) {
        return ((a<b)?a:b);
    }
    double max(double a, double b) {
        return ((a<b)?b:a);
    }

    double* matrice_frechet(double* distances, int n1, int n2)
    {
        double* frechet = new double[n1*n2];
        int i = 0, j = 0, jmin = 0;


        for (int i=0; i < n1*n2; i++) {
            frechet[i] = -1;
        }

        double mini;
        for(; i < n1; i++){
            int i2 = i * n2;
            for(j = 0;  j < n2; ) {
                if( distances[i2 + j] == -1 ) j++;
                else break;
            }
            jmin = j;
            int iM_jM, iM_j, i_jM, i_j;
            for( ; j < n2 && distances[i2 + j] != -1; j++){
                int mini;
                if( i > 0){
                    i_j = i2 + j;
                    iM_j = i_j - n2;
                    if( j >  0){
                        iM_jM = iM_j - 1;
                        i_jM = i_j - 1;
                        if( distances[iM_jM] != -1 ){
                            mini = distances[iM_jM];
                            if( distances[iM_j] != -1 ) mini = fmin(mini, distances[iM_j]);
                            if( distances[i_jM] != -1 ) mini = fmin(mini, distances[i_jM]);
                            frechet[i_j] = fmax(distances[i_j], mini);
                        } else if( distances[iM_j] != -1 ){
                            mini = distances[iM_j];
                            if( distances[i_jM] != -1 ) mini = fmin(mini, distances[i_jM]);
                            frechet[i_j] = fmax(distances[i_j], mini);
                        }else if( distances[i_jM] != -1 ) frechet[i_j] = fmax(distances[i_j], distances[i_jM]);
                        else frechet[i_j] = distances[i_j];

                    }else if( distances[iM_j] != -1 ) frechet[i_j] = fmax(distances[i_j], distances[iM_j]);
                    else frechet[i_j] = distances[i_j];
                }else if( j >  0){
                    i_jM = j - 1;
                    i_j = j;
                    if( distances[i_jM] != -1 ) frechet[i_j] = fmax(distances[i_j], distances[i_jM]);
                    else frechet[i_j] = distances[i_j];
                }else frechet[0] = distances[0];
            }
        }
        return frechet;

    }

    int main(int argc, char*argv[]) {
        int n_traj1, n_traj2;
        double* traj1 = lecture(argv[1],&n_traj1);
        double* traj2 = lecture(argv[2],&n_traj2);

        for (int i=0; i<n_traj1; i++)
            cout << traj1[2*i] << " " << traj1[2*i+1] << endl;
        for (int i=0; i<n_traj2; i++)
            cout << traj2[2*i] << " " << traj2[2*i+1] << endl;

        double* mat_dist = matrice_distance(traj1,n_traj1,traj2,n_traj2);

        cout << "la matrice de distance" << endl;

        for (int i=0; i<n_traj1; i++) {
            for (int j = 0; j < n_traj2; j++)
                cout << mat_dist[i * n_traj2 + j] << " ";
            cout << endl;
        }

        double* mat_dist2 = matrice_distanceOpti(traj1,n_traj1,traj2,n_traj2);

        cout << "la matrice de distance opti" << endl;
        for (int i=0; i<n_traj1; i++) {
            for (int j = 0; j < n_traj2; j++)
                cout << mat_dist2[i * n_traj2 + j] << " ";
            cout << endl;
        }

        double* frechet = matrice_frechet(mat_dist,n_traj1,n_traj2);

        cout << "la matrice de fréchet" << endl;

        for (int i=0; i<n_traj1; i++) {
            for (int j = 0; j < n_traj2; j++)
                cout << frechet[i * n_traj2 + j] << " ";
            cout << endl;
        }
        
        double* frechet2 = matrice_frechet(mat_dist2,n_traj1,n_traj2);

        cout << "la matrice de fréchet Opti" << endl;

        for (int i=0; i<n_traj1; i++) {
            for (int j = 0; j < n_traj2; j++)
                cout << frechet2[i * n_traj2 + j] << " ";
            cout << endl;
        }

        return 0;
    }

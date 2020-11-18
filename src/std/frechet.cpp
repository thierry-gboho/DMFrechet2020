#include <iostream>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

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


double min(double a, double b) {
    return ((a<b)?a:b);
}
double max(double a, double b) {
    return ((a<b)?b:a);
}

double* matrice_frechet(double* distances, int n1, int n2)
{
    double* frechet = new double[n1*n2];

    frechet[0] = distances[0];

    for (int i=1; i<n1; i++)
        frechet[i*n2] = max(distances[i*n2],frechet[(i-1)*n2]);

    for (int j=1; j<n2; j++)
        frechet[j] = max(distances[j], frechet[(j-1)]);

    for (int i=1; i<n1; i++)
        for (int j=1; j<n2; j++)
            frechet[i*n2+j] = max(distances[i*n2+j],min(frechet[(i-1)*n2+j],min(frechet[i*n2+(j-1)],frechet[(i-1)*n2+(j-1)])));

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


    double* frechet = matrice_frechet(mat_dist,n_traj1,n_traj2);

    cout << "la matrice de fréchet" << endl;

    for (int i=0; i<n_traj1; i++) {
        for (int j = 0; j < n_traj2; j++)
            cout << frechet[i * n_traj2 + j] << " ";
        cout << endl;
    }

    return 0;
}

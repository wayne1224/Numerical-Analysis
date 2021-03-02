// different boundary or source

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

#define EPSILON 0.00001
#define N 21
#define WIDTH 1.0
#define SOURCE 0.1

void print_matrix(vector<vector<double>>);
void finite(vector<vector<double>>);
bool norm_inf(vector<vector<double>> , vector<vector<double>>);

int main(){
    vector<vector<double>> T(N);

    for(int i = 0 ; i < N ; i++){
        T[i].resize(N);

        for(int j = 0 ; j < N ; j++){
            T[i][j] = 0;
        }
    }

    finite(T);

    system("pause");
    return 0;
}

void print_matrix(vector<vector<double>> a){
    for(int i = 0 ; i < N ; i++){
        for(int j = 0 ; j < N ; j++){
            printf("%8.5f " , a[i][j]);
        }
        cout << ";" << endl;
    }
}

void finite(vector<vector<double>> T){ 
    double s = SOURCE;
    double h = WIDTH / (N - 1);
    double d = 0;
    double w = 1;  
    int times = 0;

    while(true){
        vector<vector<double>> Told = T;

        for(int i = 1 ; i < N - 1 ; i++){ // (1 , 1) ~ (N - 2 , N - 2)
            for(int j = 1 ; j < N - 1 ; j++){                
                if(i == (N - 1) / 2 && j == (N - 1) / 2){ // s(x , y) = 5
                    T[i][j] = T[i][j] + (w / 4) * (-s * h * h + (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) - (4 * T[i][j]));
                }
                else{ // s(x , y) = 0
                    T[i][j] = T[i][j] + (w / 4) * (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1] - (4 * T[i][j]));
                }
            }      
        }

        for(int i = 0 ; i < N ; i++){
            T[i][0] = 30;     // left boundary
            T[i][N - 1] = 30; // right boundary
            T[0][i] = 20;     // bottom boundary
            // T[N - 1][i] = 30; // top boundary
        }

        for(int i = 1 ; i < N - 1 ; i++){ // top boundary
            T[N - 1][i] = T[N - 1][i] + w * (T[N - 2][i] + h * d - T[N - 1][i]);
        }  
      
        times = times + 1; 

        if(norm_inf(Told , T)){
            break;
        }
    }
    cout << endl;
    print_matrix(T);
    cout << endl;
    cout << times << endl;
}

bool norm_inf(vector<vector<double>> Told , vector<vector<double>> T){ // max row sum of entries
    double row_sum[N - 2];
    
    for(int i = 1 ; i < N - 1 ; i++){ // (1 , 1) ~ (N - 1 , N - 1)
        double tmp = 0;

        for(int j = 1 ; j < N - 1 ; j++){ 
            tmp = tmp + (T[i][j] - Told[i][j]);
        }

        row_sum[i] = fabs(tmp);
    }

    double max = 0; 

    for(int i = 0 ; i < N - 2 ; i++){
        if(max < row_sum[i]){
            max = row_sum[i];
        }
    }

    return (max <= EPSILON);
}
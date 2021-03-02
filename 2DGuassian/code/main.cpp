#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

#define PI 3.14159265358979323846
#define ANS 0.160429671237858

// f(x , y) = (sin(2 * pi * x) / (2 * pi * x)) * (sin(3 * pi * y) / (3 * pi * y))

vector<vector<double>> divide(double[] , int);
double fxy(double , double);
double Gaussian(int , vector<double>);
void requirementB(double[]);
void requirementC(double[]);
void requirementD(double[]);
void print_point(double[]);

double Root[7][7] = {{0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0}, // N = 1
                    {0.5773502691896257 , -0.5773502691896257 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0}, // N = 2
                    {0.0, 0.7745966692414834 , -0.7745966692414834 , 0.0 , 0.0 , 0.0 , 0.0}, // N = 3
                    {0.3399810435848563 , -0.3399810435848563 , 0.8611363115940526 , -0.8611363115940526 , 0.0 , 0.0 , 0.0}, // N = 4
                    {0.0 , -0.5384693101056831 , 0.5384693101056831 , -0.9061798459386640 , 0.9061798459386640 , 0.0 , 0.0}, // N = 5
                    {0.6612093864662645 , -0.6612093864662645 , -0.2386191860831969 , 0.2386191860831969 , -0.9324695142031521 , 0.9324695142031521 , 0.0}, // N = 6
                    {0.0 , 0.4058451513773972 , -0.4058451513773972 , -0.7415311855993945 , 0.7415311855993945 , -0.9491079123427585 , 0.9491079123427585}  // N = 7
                    };

double weight[7][7] = { {2.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0}, 
                        {1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0}, 
                        {0.888888888888888888888, 0.555555555555555556, 0.555555555555555556, 0.0 , 0.0 , 0.0 , 0.0},
                        {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538 , 0.0 , 0.0 , 0.0},
                        {0.5688888888888889 , 0.4786286704993665 , 0.4786286704993665 , 0.2369268850561891 , 0.2369268850561891 , 0.0 , 0.0},
                        {0.3607615730481386 , 0.3607615730481386 , 0.4679139345726910 , 0.4679139345726910 , 0.1713244923791704 , 0.1713244923791704 , 0.0},
                        {0.4179591836734694 , 0.3818300505051189 , 0.3818300505051189 , 0.2797053914892766 , 0.2797053914892766 , 0.1294849661688697 , 0.1294849661688697},
                        };

int main(){ 
    double domain[4] = {-1 , -1 , 1 , 1}; // xmin ymin xmax ymax
    
    requirementB(domain);
    requirementC(domain);    
    requirementD(domain);
    
    //print_point(domain); 

    system("pause");
    return 0;
}

vector<vector<double>> divide(double domain[] , int interval){ 
    vector<vector<double>> domains(interval);

    for(int i = 0 ; i < interval ; i++){ 
        domains[i].resize(4);
    }

    int k = 0;
    double xh = (domain[2] - domain[0]) / sqrt(interval); // x方向每個間隔的距離
    double yh = (domain[3] - domain[1]) / sqrt(interval); // y方向每個間隔的距離

    // domain (xmin ymin xmax ymax)
    // (xmin , ymin) = (domain[0] , domain[1]) 最左下的點

    for(int i = 0 ; i < sqrt(interval) ; i++){
        for(int j = 0 ; j < sqrt(interval) ; j++){      
            // 間隔裡的 xmin ymin xmax ymax
            domains[k][0] = domain[0] + xh * j;        
            domains[k][1] = domain[1] + yh * i;        
            domains[k][2] = domain[0] + xh * (j + 1);  
            domains[k][3] = domain[1] + yh * (i + 1);  
 
            k = k + 1;
        }
    }

    // cout << " xmin  ymin  xmax  ymax" << endl;

    // for(int i = 0 ; i < interval ; i++){
    //     for(int j = 0 ; j < 4 ; j++){
    //         cout << std::setw(5) << std::fixed << std::setprecision(2) << domains[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    return domains;
}

double fxy(double x , double y){
    double a = 2 * PI * x;
    double b = 3 * PI * y;

    if(a == 0){
        a = 1;
    }
    else{
        a = sin(a) / (a);
    }

    if(b == 0){
        b = 1;
    }
    else{
        b = sin(b) / (b);
    }

    return a * b;
}

double Gaussian(int N , vector<double> domain){ // N = sample points的數量 , domain = 積分範圍 {xmin ymin xmax ymax}
    double sum = 0;

    for(int i = 0 ; i < N ; i++){  // sum算出來是一個interval的積分
        for(int j = 0 ; j < N ; j++){
            double x = (((domain[2] - domain[0]) / 2) * Root[N - 1][i]) + ((domain[0] + domain[2]) / 2);
            double y = (((domain[3] - domain[1]) / 2) * Root[N - 1][j]) + ((domain[1] + domain[3]) / 2);
            sum = sum + (weight[N - 1][i] * weight[N - 1][j] * fxy(x , y));
            // cout << fxy(x , y) << endl;
            // cout << (weight[N - 1][i] * weight[N - 1][j] * fxy(x , y)) << endl;
        }
    }
    
    sum = sum * ((domain[2] - domain[0]) / 2) * ((domain[3] - domain[1]) / 2);

    return sum;
}

void requirementB(double domain[4]){
    double result[4][3] = {0};

    for(int i = 1 ; i <= 4 ; i++){
        int interval = pow(i , 2);
        vector<vector<double>> domains = divide(domain , interval); // 每個間格的domain

        for(int j = 2 ; j <= 4 ; j++){   // N = 2 ~ 4
            double sum = 0;

            for(int k = 0 ; k < domains.size() ; k++){
                sum = sum + Gaussian(j , domains[k]); 
            }
    
            result[i - 1][j - 1] = sum; 
        }
    }

    for(int i = 0 ; i < 4 ; i++){
        cout << "interval = " << std::setw(2) << std::fixed << std::setprecision(0) << pow(i + 1 , 2) << "       error" << endl;

        for(int j = 0 ; j < 3 ; j++){
            cout << std::setw(10) << std::fixed << std::setprecision(16) << result[i][j] << "  " << abs(result[i][j] - ANS) << endl;
        }
        
        cout << endl;
    } 
}

void requirementC(double domain[4]){
    double result[7][3] = {0};

    for(int i = 1 ; i <= 7 ; i++){  // interval = 1 ~ 49
        int interval = i * i;
        vector<vector<double>> domains = divide(domain , interval); // 每個間格的domain

        for(int j = 2 ; j <= 4 ; j++){   // N = 2 ~ 4
            double sum = 0;

            for(int k = 0 ; k < domains.size() ; k++){
                sum = sum + Gaussian(j , domains[k]); 
            }
    
            result[i - 1][j - 1] = sum; 
        }
    }

    for(int i = 0 ; i < 7 ; i++){
        cout << "interval = " << std::setw(2) << std::fixed << std::setprecision(0) << pow(i + 1 , 2) << "       error" << endl;

        for(int j = 0 ; j < 3 ; j++){
            cout << std::setw(10) << std::fixed << std::setprecision(16) << result[i][j] << "  " << abs(result[i][j] - ANS) << endl;
        }

        cout << endl;
    } 
}

void requirementD(double domain[4]){
    double result[7][7] = {0};

    for(int i = 1 ; i <= 7 ; i++){  // interval = 1 ~ 49
        int interval = i * i;
        vector<vector<double>> domains = divide(domain , interval); // 每個間格的domain

        for(int j = 1 ; j <= 7 ; j++){   // N = 1 ~ 7
            double sum = 0;

            for(int k = 0 ; k < domains.size() ; k++){
                sum = sum + Gaussian(j , domains[k]); 
            }
    
            result[i - 1][j - 1] = sum; 
        }
    }

    for(int i = 0 ; i < 7 ; i++){
        cout << "interval = " << std::setw(2) << std::fixed << std::setprecision(0) << pow(i + 1 , 2) << "       error" << endl;

        for(int j = 0 ; j < 7 ; j++){
            cout << std::setw(10) << std::fixed << std::setprecision(16) << result[i][j] << "  " << abs(result[i][j] - ANS) << endl;
        }

        cout << endl;
    } 
}

void print_point(double domain[4]){
    double result[7][7] = {0};

    for(int i = 1 ; i <= 7 ; i++){  // interval = 1 ~ 49
        int interval = i * i;
        vector<vector<double>> domains = divide(domain , interval); // 每個間格的domain

        for(int j = 1 ; j <= 7 ; j++){   // N = 1 ~ 7
            double sum = 0;

            for(int k = 0 ; k < domains.size() ; k++){
                sum = sum + Gaussian(j , domains[k]); 
            }
    
            result[i - 1][j - 1] = sum; 
        }
    }

    for(int i = 0 ; i < 7 ; i++){
        // cout << i + 1 << endl;
        for(int j = 0 ; j < 7 ; j++){
            // cout << j + 1 << endl;
            cout << std::fixed << std::setprecision(16) << "  " << abs(result[j][i] - ANS) / ANS << endl;
        }

        cout << endl;
    } 
} 

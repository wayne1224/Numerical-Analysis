#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

void methodA(vector<double> , vector<double>);
void methodB();

int main(){
    vector<double> x = {9.5 , 5 , 3.5 , 1.5 , 3 , 5.5 , 10 , 7 , 6.5 , 9.5};
    vector<double> y = {24 , 22 , 20 , 12 , 7 , 4 , 2 , 6 , 18 , 24};
    
    cout << " x = ";

    for(int i = 0 ; i < 10 ; i++){
        cout << std::setw(8) << x[i] << " ";
    }
    cout << endl;

    cout << " y = ";

    for(int i = 0 ; i < 10 ; i++){
        cout << std::setw(8) << y[i] << " ";
    }
    cout << endl;

    methodA(x , y);
    methodB();

    system("pause");
    return 0;
}

void methodA(vector<double> x , vector<double> y){ // Chord - length
    vector<double> t(10 , 0);
    vector<double> l(10 , 0);

    for(int i = 1 ; i < 10 ; i++){
        l[i - 1] = sqrt(pow(x[i - 1] - x[i] , 2) + pow(y[i - 1] - y[i] , 2));
        t[i] = t[i - 1] + l[i - 1];
    }

    cout << "t1 = ";

    for(auto i : t){ 
        cout << std::setw(8) << std::fixed << std::setprecision(5) << i << " ";
    }
    cout << endl;
}

void methodB(){ // Uniform
    cout << "t2 = ";
    
    for(double i = 0 ; i < 10 ; i++){
        cout << std::setw(8) << std::fixed << std::setprecision(5) << i << " ";
    }
}

/*

 x =      9.5        5      3.5      1.5        3      5.5       10        7      6.5      9.5
 y =       24       22       20       12        7        4        2        6       18       24
t1 =  0.00000  4.92443  7.42443 15.67064 20.89079 24.79592 29.72035 34.72035 46.73076 53.43896
t2 =  0.00000  1.00000  2.00000  3.00000  4.00000  5.00000  6.00000  7.00000  8.00000  9.00000

*/

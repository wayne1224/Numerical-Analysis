#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::pair;

vector<double> methodA(vector<double> , vector<double>);
vector<double> methodB();
vector<vector<double>> divided_difference(vector<double> , vector<double>);
double Horner(double , vector<double> , vector<vector<double>>); // 給參數t會有 x or y 的位置

int main(){
    vector<double> x = {9.5 , 5 , 3.5 , 1.5 , 3 , 5.5 , 10 , 7 , 6.5 , 9.5};
    vector<double> y = {24 , 22 , 20 , 12 , 7 , 4 , 2 , 6 , 18 , 24}; 
    vector<double> t1;
    vector<double> t2;
    vector<vector<double>> df1x;
    vector<vector<double>> df1y;
    vector<vector<double>> df2x;
    vector<vector<double>> df2y;

    t1 = methodA(x , y);
    t2 = methodB();
    df1x = divided_difference(t1 , x);
    df1y = divided_difference(t1 , y);
    df2x = divided_difference(t2 , x);
    df2y = divided_difference(t2 , y);

    vector<double> pointX;
    vector<double> pointY;
    vector<double> pointT;

    cout << "t1" << endl;

    for(int i = 0 ; i < 54 ; i++){
        pointX.push_back(Horner(1 * i , t1 , df1x));
        pointY.push_back(Horner(1 * i , t1 , df1y));
        pointT.push_back(1 * i); 
    }

    for(int i = 0 ; i < pointX.size() ; i++){
        cout << pointT[i] << " " << pointX[i] << endl;
    }

    cout << endl;

    for(int i = 0 ; i < pointY.size() ; i++){
        cout << pointT[i] << " " << pointY[i] << endl;
    }


    pointX.clear();
    pointY.clear();
    pointT.clear();

    cout << endl << "t2" << endl;

    for(int i = 0 ; i < 40 ; i++){
        pointX.push_back(Horner(0.225 * i , t2 , df2x));
        pointY.push_back(Horner(0.225 * i , t2 , df2y));
        pointT.push_back(0.225 * i); 
    }

    for(int i = 0 ; i < pointX.size() ; i++){
        cout << pointT[i] << " " << pointX[i] << endl;
    }

    cout << endl;

    for(int i = 0 ; i < pointY.size() ; i++){
        cout << pointT[i] << " " << pointY[i] << endl;
    }

    system("pause");
    return 0;
}

vector<double> methodA(vector<double> x , vector<double> y){ // Chord - length
    vector<double> t(10 , 0);
    vector<double> l(10 , 0);

    for(int i = 1 ; i < 10 ; i++){
        l[i - 1] = sqrt(pow(x[i - 1] - x[i] , 2) + pow(y[i - 1] - y[i] , 2));
        t[i] = t[i - 1] + l[i - 1];
    }

    return t;
}

vector<double> methodB(){ // Uniform
    vector<double> t(10 , 0);

    for(int i = 0 ; i < 10 ; i++){
        t[i] = i;
    }

    return t;
}

vector<vector<double>> divided_difference(vector<double> x , vector<double> y){ // 2D 
    vector<vector<double>> df(10);

    for(int i = 0 ; i < 10 ; i++){
        df[i].resize(10);
    }
    
    for(int i = 0 ; i < 10 ; i++){
        df[0][i] = y[i];
    }

    for(int i = 1 ; i < 10 ; i++){
        for(int j = 0 ; j < 10 - i ; j++){
            df[i][j] = (df[i - 1][j + 1] - df[i - 1][j]) / (x[i + j] - x[j]);
        }
    }

    return df;
}

double Horner(double t , vector<double> tt , vector<vector<double>> coef){ // t = 參數 , tt = methodA or methodB 的參數陣列 , coef = 係數 
    double sum = coef[9][0];
    
    for(int i = 8 ; i >= 0 ; i--){
        sum = sum * (t - tt[i]) + coef[i][0];
    }
 
    return sum;
}

/*
用多項式產生的 x , y ， 與原本 x , y 做比較

set xrange[0 : 60]
set yrange[-20 : 80]
plot 't1_TtoX.txt' with line , 't1andx.txt' with line
plot 't1_TtoY.txt' with line , 't1andy.txt' with line

set xrange[0 : 10]
set yrange[0 : 30]
plot 't2_TtoX.txt' with line , 't2andx.txt' with line
plot 't2_TtoY.txt' with line , 't2andy.txt' with line

*/
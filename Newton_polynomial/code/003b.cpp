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
    vector<double> y = {50 , 22 , 20 , 12 , 7 , 4 , 2 , 6 , 18 , 50}; 
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

    cout << "t1 points" << endl;

    for(int i = 0 ; i < 40 ; i++){  // t = 0 - 201.9
        double t = 2.5 * i;
        double xx = Horner(t , t1 , df1x);
        double yy = Horner(t , t1 , df1y);

        // cout << std::setw(3) << t << "  ";
        cout << xx << "  " << yy << endl;      
    }

    cout << Horner(t1[9] , t1 , df1x) << "  " << Horner(t1[9] , t1 , df1y) << endl << endl;

    ///////////////////////////////////////////////
    
    cout << "t2 points" << endl;

    for(int i = 0 ; i < 40 ; i++){ // t = 0 - 9
        double t = 0.225 * i;
        double xx = Horner(t , t2 , df2x);
        double yy = Horner(t , t2 , df2y);

        // cout << std::setw(3) << t << "  ";
        cout << xx << "  " << yy << endl;     
    }

    cout << Horner(t2[9] , t2 , df2x) << "  " << Horner(t2[9] , t2 , df2y) << endl << endl;

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

    cout << "t1 = ";

    for(auto i : t){ 
        cout << std::setw(8) << std::fixed << std::setprecision(5) << i << " ";
    }
    cout << endl;

    return t;
}

vector<double> methodB(){ // Uniform
    vector<double> t(10 , 0);

    for(int i = 0 ; i < 10 ; i++){
        t[i] = i;
    }

    cout << "t2 = ";

    for(double i = 0 ; i < 10 ; i++){
        cout << std::setw(8) << std::fixed << std::setprecision(5) << i << " ";
    }
    cout << endl;

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

set xrange[-2 : 12]
set yrange[-20 : 70]
set xtics 1.0
set ytics 10.0
plot 'moon2.txt' with linespoints pointtype 7

t1 =  0.00000 28.35930 30.85930 39.10551 44.32567 48.23079 53.15522 58.15522 70.16563 102.30595
t2 =  0.00000  1.00000  2.00000  3.00000  4.00000  5.00000  6.00000  7.00000  8.00000  9.00000

vector<double> pointX;
    vector<double> pointY;
    vector<double> pointT;

    cout << "t1" << endl;

    for(int i = 0 ; i < 40 ; i++){
        pointX.push_back(Horner(5 * i , t1 , df1x));
        pointY.push_back(Horner(5 * i , t1 , df1y));
        pointT.push_back(5 * i); 
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



*/
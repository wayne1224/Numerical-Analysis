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
void print_point(vector<double> , vector<double> , vector<vector<double>> , vector<vector<double>> , vector<vector<double>> , vector<vector<double>>);

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
    print_point(t1 , t2 , df1x , df1y , df2x , df2y);

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

void print_point(vector<double> t1 , vector<double> t2 , vector<vector<double>> df1x , vector<vector<double>> df1y , vector<vector<double>> df2x , vector<vector<double>> df2y){
    cout << "t1 divided_difference" << endl;

    cout << "t1 to x =  " ;

    for(auto i : df1x){
        cout << i[0] << " ";
    }
    cout << endl;

    cout << "t1 to y =  " ;

    for(auto i : df1y){
        cout << i[0] << " ";
    }
    cout << endl << endl;

    ///////////////////////////////////////////////

    cout << "t2 divided_difference" << endl;

    cout << "t2 to x =  " ;

    for(auto i : df2x){
        cout << i[0] << " ";
    }
    cout << endl;

    cout << "t2 to y =  " ;

    for(auto i : df2y){
        cout << i[0] << " ";
    }
    cout << endl << endl;

    ///////////////////////////////////////////////

    cout << "t1 points" << endl;

    for(int i = 0 ; i < 54 ; i++){ // t = 0 - 53.4
        double t = 1 * i;
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
}

/*

 x =      9.5        5      3.5      1.5        3      5.5       10        7      6.5      9.5
 y =       24       22       20       12        7        4        2        6       18       24
t1 =  0.00000  4.92443  7.42443 15.67064 20.89079 24.79592 29.72035 34.72035 46.73076 53.43896
t2 =  0.00000  1.00000  2.00000  3.00000  4.00000  5.00000  6.00000  7.00000  8.00000  9.00000

t1 divided_difference
9.5 -0.913812 0.0422674 -0.000574527 4.57429e-005 -2.69795e-006 8.86261e-008 -2.53391e-008 1.88749e-009 -9.35578e-011
24 -0.406138 -0.0530494 0.00237493 -6.34734e-005 2.75169e-006 -8.09994e-008 9.67657e-009 -7.13176e-010 3.71339e-011

t2 divided_difference
9.5 -4.5 1.5 -0.583333 0.3125 -0.116667 0.0333333 -0.00952381 0.00287698 -0.000789517
24 -2 0 -1 0.625 -0.208333 0.0486111 -0.0077381 0.000694444 -1.65344e-005


set xrange[-5 : 20]
set yrange[0 : 25]
set xtics 1.0
set ytics 1.0
plot 'moon.txt' with linespoints pointtype 7


t1 gnuplot set
set xrange[-20 : 65]
set yrange[0 : 25]
set xtics 5.0
set ytics 1.0
plot 't1.txt' with linespoints pointtype 14 , 'moon.txt' with linespoints pointtype 7


t2 gnuplot set
set xrange[-5 : 20]
set yrange[0 : 25]
set xtics 5.0
set ytics 1.0
plot 't2.txt' with linespoints pointtype 14 , 'moon.txt' with linespoints pointtype 7

*/
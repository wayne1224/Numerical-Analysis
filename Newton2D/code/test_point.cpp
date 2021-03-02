#include <iostream>
#include <math.h>
#include <iomanip>

#define EPSILON 0.000001

using std::cout;
using std::cin;
using std::endl;
using std::fixed;
using std::setprecision;

// Function : f(x , y) = (x - 3)^2 / 25 + (y - 4)^2 / 25 - 1 = 0
//            g(x , y) = -x^2 + 6x + y - 10 = 0

void newton2D(double , double);
double f(double , double);
double g(double , double);
double fx(double , double);
double fy(double , double);
double gx(double , double);
double gy(double , double);
void print(int , double , double , double);

int main(){
    double x = 8 , y = 9;

    while(y >= -1){ // x的範圍 -2 ~ 8 y的範圍 -1 ~ 9 
        while(x >= -2){
            newton2D(x , y);
            x = x - 0.25;
        }

        cout << "  " << fixed << setprecision(2) << y << endl;

        x = 8;
        y = y - 0.25;    
    }

    

    newton2D(x , y);

    system("pause");
    return 0;
}

void newton2D(double x , double y){
    double error , xnew , ynew , xold , yold , fold , gold;
    double fdx , fdy , gdx , gdy , delta;
    int i = 1;
    
    xold = x; // X0
    yold = y; // Y0
    
    while(1){
        fold = f(xold, yold);
        gold = g(xold, yold);
        fdx = fx(xold, yold);
        fdy = fy(xold, yold);
        gdx = gx(xold, yold);
        gdy = gy(xold, yold);

        delta = fdy * gdx - fdx * gdy; 
        xnew = xold + (gdy * fold - fdy * gold) / delta;  // Xn+1 = Xn + h  
        ynew = yold + (-gdx * fold + fdx * gold) / delta; // Yn+1 = Yn + k
        error = pow(pow(xnew - xold , 2) + pow(ynew - yold , 2) , 0.5);
        

        if(error < EPSILON){
            cout << " o ";
            break;
        }

        xold = xnew; 
        yold = ynew;
        i = i + 1;

        if(i > 200){
            cout << " x ";
            break;
        }
    }

} 

double f(double x , double y){   // Function : f(x , y) = (x - 3)^2 / 25 + (y - 4)^2 / 25 - 1 = 0 
    return (pow(x - 3.0 , 2) / 25.0) + (pow(y - 4.0 , 2) / 25.0) - 1.0;
}

double g(double x , double y){   // Function : g(x , y) = -x^2 + 6x + y - 10 = 0
    return -pow(x , 2) + (6.0 * x) + y - 10.0;
}

double fx(double x , double y){
    return (2.0 * x - 6.0) / 25.0;
}

double fy(double x , double y){
    return (2.0 * y - 8.0) / 25.0;
}

double gx(double x , double y){
    return (-2.0 * x) + 6.0;
}

double gy(double x , double y){
    return 1.0;
}


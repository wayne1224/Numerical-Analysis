#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <time.h>  
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

#define EPSILON 0.000001
#define PI 3.141592653

class matrix{
    public:
        vector<vector<double>> data;
        int row_size;
        int column_size;
        
        matrix(int , int);
        void print_matrix();
        void identity();

        matrix operator-(const matrix &);
        matrix operator*(const matrix &);
};

void generate_data(matrix &);
void Jacobi(matrix , matrix);
double max_off_diagonal(matrix , int &, int &);
double sum_off_diagonal(matrix);
void A_update(matrix & , int , int , double , double);
void P_update(matrix & , int , int , double , double);
void wFile(vector<double> , vector<matrix> , int);

int main(){
    int N = 10;
    
    class matrix A(N , N);
    class matrix P(N , N);

    generate_data(A);
    P.identity();

    // cout << "matrix A :" << endl;
    // A.print_matrix();

    Jacobi(A , P);

    system("pause");
    return 0;
}

matrix::matrix(int row , int column){
    data.resize(row);

    for(int i = 0 ; i < row ; i++){
        data[i].resize(column);
    }

    for(int i = 0 ; i < row ; i++){
        for(int j = 0 ; j < column ; j++){ 
            data[i][j] = 0;
        }
    }

    row_size = row;
    column_size = column;
}

void matrix::print_matrix(){
    for(int i = 0 ; i < row_size ; i++){
        for(int j = 0 ; j < column_size ; j++){ 
            printf("%9.5f " , data[i][j]);
        }
        cout << endl;
    }

    cout << endl;
}

void matrix::identity(){
    for(int i = 0 ; i < row_size ; i++){
        data[i][i] = 1;
    }
}

matrix matrix::operator-(const matrix &other){
    class matrix result(row_size , column_size);

    for(int i = 0 ; i < row_size ; i++){
        for (int j = 0 ; j < column_size ; j++){
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }

    return result;
}

matrix matrix::operator*(const matrix &other){
    class matrix result(row_size , other.column_size);

    for (int i = 0 ; i < row_size ; i++) {
        for (int j = 0 ; j < other.column_size ; j++) {
            for(int k = 0 ; k < column_size ; k++) {
                result.data[i][j] = result.data[i][j] + (data[i][k] * other.data[k][j]); 
            }
        }
    }

    return result;
}

void generate_data(matrix &A){
    srand((unsigned)time(NULL));

    for(int i = 0 ; i < A.row_size ; i++){
        for(int j = i ; j < A.column_size ; j++){
            A.data[i][j] = A.data[j][i] = (double)(rand() % 20);
        }
    }
}

void Jacobi(matrix A , matrix P){
    int p , q;
    int times = 0;
    class matrix AA = A;
    vector<double> offD;
    vector<matrix> As;

    double Apq = max_off_diagonal(A , p , q);

    while(fabs(Apq) > EPSILON){
        double theta;
        double a = 2.0 * Apq;
		double b = A.data[p][p] - A.data[q][q];

        theta = atan(a / b) / 2;

        double c = cos(theta);
		double s = sin(theta);

        A_update(A , p , q , c , s);
        P_update(P , p , q , c , s);

        As.push_back(A);
        offD.push_back(sum_off_diagonal(A));

        Apq = max_off_diagonal(A , p , q);
        times = times + 1;
    }

    // cout << "times = " << times << endl;

    // for(int i = 0 ; i < offD.size() ; i++){
    //     printf("%11.5f " , offD[i]);
    //     if((i + 1) % 10 == 0 or i == offD.size() - 1){
    //         cout << endl;
    //     }
    // }

    for(int i = 0 ; i < offD.size() ; i++){
        printf("%3d %15.9f \n" , i , offD[i]);
    }

    //while(1);

    cout << endl;
}

double max_off_diagonal(matrix A , int &p , int &q){
    double max = A.data[0][1];
	p = 0;
	q = 1;

	for(int i = 0 ; i < A.row_size ; i++){
        for(int j = i + 1 ; j < A.column_size ; j++){
            if(fabs(A.data[i][j]) > fabs(max)){
				p = i;
				q = j;
				max = A.data[i][j];
			}
        }		
    }
		
	return max;
}

double sum_off_diagonal(matrix A){
    double sum = 0;

    for(int i = 0 ; i < A.row_size ; i++){
        for(int j = i + 1 ; j < A.column_size ; j++){
            sum = sum + fabs(A.data[i][j]);
        }		
    }

    return sum;
}

void A_update(matrix &A , int p , int q , double c , double s){
    double Bp[A.row_size] = {0};
    double Bq[A.row_size] = {0};

    for(int i = 0 ; i < A.row_size ; i++){
        if(i != p && i != q){
            Bp[i] = (c * A.data[i][p]) + (s * A.data[i][q]);
            Bq[i] = (-s * A.data[i][p]) + (c * A.data[i][q]);
        }
    }

    Bp[p] = (c * c * A.data[p][p]) + (2 * s * c * A.data[p][q]) + (s * s * A.data[q][q]);
    Bq[q] = (s * s * A.data[p][p]) - (2 * s * c * A.data[p][q]) + (c * c * A.data[q][q]);

    for(int i = 0 ; i < A.row_size ; i++){
        A.data[p][i] = A.data[i][p] = Bp[i];
        A.data[q][i] = A.data[i][q] = Bq[i];
    }

    // A.data[p][q] = A.data[q][p] = 0;
}

void P_update(matrix &P , int p , int q , double c , double s){
    double Bp[P.row_size] = {0};
    double Bq[P.row_size] = {0};

    for(int i = 0 ; i < P.row_size ; i++){
	   Bp[i] = c * P.data[i][p] + s * P.data[i][q];
	   Bq[i] = (-s) * P.data[i][p] + c * P.data[i][q];
   }

   for(int i = 0 ; i < P.row_size ; i++){
	   P.data[i][p] = Bp[i];
	   P.data[i][q] = Bq[i];
   }
}

// set xrange[0 : 140]
// set yrange[-10 : 400]
// plot '004_1.txt' with linespoints pointtype 7

// set xrange[0 : 600]
// set yrange[-100 : 2000]
// plot '004_2.txt' with linespoints pointtype 7
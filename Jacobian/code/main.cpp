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
void A_update(matrix & , int , int , double , double);
void P_update(matrix & , int , int , double , double);
double norm(matrix);
matrix check_orthogonal(matrix);

int main(){
    int N = 4;
    
    class matrix A(N , N);
    class matrix P(N , N);

    generate_data(A);
    P.identity();
    
    cout << "matrix A :" << endl;
    A.print_matrix();

    Jacobi(A , P);

    // while(1){

    // }
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
    class matrix AA = A; // original A

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

        cout << "times = " << times << endl;
        cout << "eigenvalue , A after iteration" << endl;
        A.print_matrix();
        cout << "eigenvector" << endl;
        P.print_matrix();

        class matrix Av = AA * P;
        class matrix Lv = P * A;

        cout << "Av" << endl;
        Av.print_matrix();
        cout << "Lv" << endl;
        Lv.print_matrix();
        cout << "norm = " << norm(Av - Lv) << endl << endl;

        cout << "Check Orthogonal" << endl;
        check_orthogonal(P).print_matrix();

        Apq = max_off_diagonal(A , p , q);
        times = times + 1;
    }
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
    Bp[q] = Bq[p] = ((c * c - s * s) * A.data[p][q]) + (s * c * (A.data[q][q] - A.data[p][p]));

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

double norm(matrix A){ // infinite norm
    double row_sum[A.row_size];

    for(int i = 0 ; i < A.row_size ; i++){
        double tmp = 0;

        for(int j = 0 ; j < A.row_size ; j++){
            tmp = tmp + fabs(A.data[i][j]);
        }

        row_sum[i] = tmp;
    }

    // for(int i = 0 ; i < A.row_size ; i++){
    //     cout << row_sum[i] << " ";
    // }
    // cout << endl;

    double max = 0; 

    for(int i = 0 ; i < A.row_size ; i++){
        if(max < row_sum[i]){
            max = row_sum[i];
        }
    }

    return max;
}

matrix check_orthogonal(matrix v){
    int n = v.row_size;

    class matrix result(n , n);

    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < n ; j++){
            double tmp = 0;

            for(int k = 0 ; k < n ; k++){
                tmp = tmp + (v.data[i][k] * v.data[j][k]);
            }

            if(tmp < 0.00000000000000001){ // tmp < 10 ^ -17
                tmp = 0;
            }

            result.data[i][j] = sqrt(tmp);    
        }
    }
    return result;
}


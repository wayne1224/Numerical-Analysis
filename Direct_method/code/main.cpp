#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::vector;

class matrix{
    public:
        vector<vector<double>> data;
        int row_size;
        int column_size;
        
        matrix(int , int);
        void print_matrix();
        void identity();
        
        matrix operator+(const matrix &); // 只有同列行數的才能加減 不然會錯
        matrix operator-(const matrix &);
        matrix operator*(const matrix &);
};

void requirement();

void generate_data(matrix & , matrix & , matrix & , int , int);
double Horner(double , int);

matrix gaussian_elimination(matrix , matrix);
matrix QR_decomposition(matrix , matrix);
matrix backward_substitution(matrix , matrix);

void compare_error(matrix , matrix , matrix);
double norm_2(matrix);
double norm_inf(matrix);

int main(){
    cout << std::fixed << std::setprecision(12);

    requirement();

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
            cout << std::setw(8) << data[i][j] << " ";
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

matrix matrix::operator+(const matrix &other){
    class matrix result(row_size , column_size);

    for(int i = 0 ; i < row_size ; i++){
        for (int j = 0 ; j < column_size ; j++){
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }

    return result;
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

void requirement(){
    for(int i = 7 ; i <= 12 ; i++){  // i = degrees , i + 4 = sample points      
        class matrix A(i + 4 , i + 1);
        class matrix AT(i + 1 , i + 4);
        class matrix Y(i + 4 , 1);

        generate_data(A , AT , Y , i , i + 4);

        class matrix B = AT * A;
        class matrix D = AT * Y;

        class matrix result1 = gaussian_elimination(B , D);
        class matrix result2 = QR_decomposition(B , D);
        class matrix result3 = QR_decomposition(A , Y);

        // cout << "sample points = " << i + 4 << " , degrees = " << i << endl << endl;

        // cout << "result1" << endl;
        // result1.print_matrix();

        // cout << "result2" << endl;
        // result2.print_matrix();

        // cout << "result3" << endl;
        // result3.print_matrix();

        compare_error(result1 , result2 , result3);
    }
}

void generate_data(matrix &A , matrix &AT , matrix &Y , int degrees , int points){
    vector<double> x;
    vector<double> y;

    for(double i = 0 ; i < points ; i++){     
        x.push_back(2 + i * 0.2);
        y.push_back(Horner((2 + i * 0.2) , degrees - 1));
    }
    
    // cout << " x     y" << endl;

    // for(int i = 0 ; i < x.size() ; i++){
    //     cout << std::fixed << std::setprecision(2) << x[i] << " " << y[i] << endl;
    // }

    // cout << endl;

    for(int i = 0 ; i < A.row_size ; i++){
        for(int j = 0 ; j < A.column_size ; j++){
            A.data[i][j] = pow(x[i] , j);
        }
    }

    for(int i = 0 ; i < AT.row_size ; i++){
        for(int j = 0 ; j < AT.column_size ; j++){
            AT.data[i][j] = pow(x[j] , i);
        }
    }

    for(int i = 0 ; i < Y.row_size ; i++){
        Y.data[i][0] = y[i];
    }
}

double Horner(double t , int n){
    double sum = 1;

    for(int i = n ; i >= 0 ; i--){
        sum = sum * (t - 0) + 1;
    }

    return sum;
}

matrix gaussian_elimination(matrix A , matrix B){ // Ax = B
    int n = A.column_size;

    // Partial pivoting and Forward elimination

    for(int i = 0 ; i < n - 1 ; i++){      
        
        // Partial pivoting

        double max = fabs(A.data[i][i]);
        int p = i;

        for(int j = i ; j < n ; j++){
            if(fabs(A.data[j][i]) > max){
                p = j;
                max = fabs(A.data[j][i]);
            }
        }

        if(p != i){
            double tmp;

            for(int j = i ; j < n ; j++){
                tmp = A.data[p][j];
				A.data[p][j] = A.data[i][j];
				A.data[i][j] = tmp;
            }

            tmp = B.data[p][0];
			B.data[p][0] = B.data[i][0];
			B.data[i][0] = tmp;

        }


        // Forward elimination

        for(int j = i + 1 ; j < n ; j++){
            if(A.data[j][i] == 0.0){
                continue;
            }

            double r = A.data[j][i] / A.data[i][i];
            
            for(int k = i ; k < n ; k++){
                A.data[j][k] = A.data[j][k] - (r * A.data[i][k]);              
            }

            B.data[j][0] = B.data[j][0] - (r * B.data[i][0]);
        }
    }

   

    class matrix x = backward_substitution(A , B);

    return x;
}

matrix QR_decomposition(matrix A , matrix B){ 
    int n = 0;

    if(A.row_size > A.column_size){
        n = A.column_size;     // over-constraint system 消 A.column_size 個 columns
    }
    else{
        n = A.column_size - 1; // 消 A.column_size - 1 個 columns
    }

    // if A is over-constraint system ex: A = 11 x 8 , H = 11 x 11 

    for(int i = 0 ; i < n ; i++){ 
        class matrix H(A.row_size , A.row_size);
        class matrix I(A.row_size , A.row_size);
        class matrix e(A.row_size , 1);
        class matrix t(A.row_size , 1);
        class matrix v(A.row_size , 1);
        class matrix vt(1 , A.row_size);

        double tt = 0;

        for(int j = i ; j < A.row_size ; j++){
            t.data[j][0] = A.data[j][i];
            tt = tt + (t.data[j][0] * t.data[j][0]);
        }
        
        tt = sqrt(tt);
        e.data[i][0] = tt;

        if(t.data[i][0] >= 0){ // v = t + (|t| * e)
            v = t + e;
        }
        else{ // v = t - (|t| * e)
            v = t - e;
        }

        for(int j = 0 ; j < v.row_size ; j++){
            vt.data[0][j] = v.data[j][0];
        }

        I.identity();

        class matrix vvT = v * vt; // A.row_size x A.row_size
        class matrix vTv = vt * v; // 1 x 1

        for(int j = 0 ; j < vvT.row_size ; j++){
            for(int k = 0 ; k < vvT.column_size ; k++){
                vvT.data[j][k] = 2.0 * (vvT.data[j][k] / vTv.data[0][0]); // 2 * (v * vt / vt * v)
            }
        }
    
        H = I - vvT;
        A = H * A;
        B = H * B;
    }

    class matrix x = backward_substitution(A , B);

    return x;
}

matrix backward_substitution(matrix A , matrix B){
    class matrix x(A.column_size , B.column_size);
    int n = A.column_size;

    for(int i = n - 1 ; i >= 0 ; i--){
        x.data[i][0] = B.data[i][0] / A.data[i][i];

        for(int j = i - 1 ; j >= 0 ; j--){
            B.data[j][0] = B.data[j][0] - (A.data[j][i] * x.data[i][0]);
        }
    }

    return x;
}

void compare_error(matrix result1 , matrix result2 , matrix result3){
    class matrix ans(result1.row_size , 1);

    for(int i = 0 ; i < ans.row_size ; i++){
        ans.data[i][0] = 1;
    }

    class matrix error1 = ans - result1;
    class matrix error2 = ans - result2;
    class matrix error3 = ans - result3;

    // cout << "2 norm" << endl;
    // cout << "error1 = " << norm_2(error1) << endl;
    // cout << "error2 = " << norm_2(error2) << endl;
    // cout << "error3 = " << norm_2(error3) << endl << endl;

    // cout << "00 norm" << endl;
    // cout << "error1 = " << norm_inf(error1) << endl;
    // cout << "error2 = " << norm_inf(error2) << endl;
    // cout << "error3 = " << norm_inf(error3) << endl << endl;

    // cout << norm_2(error1) << endl;
    // cout << norm_2(error2) << endl;
    // cout << norm_2(error3) << endl;

    // cout << norm_inf(error1) << endl;
    // cout << norm_inf(error2) << endl;
    // cout << norm_inf(error3) << endl << endl;

}

double norm_2(matrix A){
    double result = 0;

    for(int i = 0 ; i < A.row_size ; i++){
        result = result + pow(A.data[i][0] , 2);
    }

    return sqrt(result);
}

double norm_inf(matrix A){
    double result = 0;

    for(int i = 0 ; i < A.row_size ; i++){
        if(result < fabs(A.data[i][0])){
            result = fabs(A.data[i][0]);
        }
    }

    return result;
}


/*

// cout << "A(" << j << " , " << i << ") " << A.data[j][i] << endl;
// cout << "A(" << i << " , " << i << ") " << A.data[i][i] << endl;
// cout << "r = " << r << endl;

// cout << "A(" << j << " , " << k << ") " << A.data[j][k] << endl;
// cout << "A(" << i << " , " << k << ") " << A.data[i][k] << endl;
// cout << "r * A.data[i][k] = " << r * A.data[i][k] << endl;

*/
#ifndef FUNC
#define FUNC

using namespace std;

#define cosd(x) cos((x)*M_PI/180)
#define sind(x) sin((x)*M_PI/180)
#define asind(x) asin((x))*180/M_PI

vector<vector<complex<double>>> get_D(double,double,double,double,double,double);

vector<vector<complex<double>>> get_M(double,double,double,double,double);

vector<vector<complex<double>>> matrix_product(vector<vector<complex<double>>> A,vector<vector<complex<double>>> B);

void show_matrix( vector<vector<complex<double>>>);

vector<double> theta_cal(vector<double> , double);



//---------------------------------------------------------------------------------------//

void show_matrix( vector<vector<complex<double>>> matrix){

    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            cout<<matrix[i][j]<<"|";
        }
        cout<<endl;
    }cout<<endl;
}
vector<double> theta_cal(vector<double> n , double incident_angle){
    vector<double> theta(n.size());
    theta[0]=incident_angle;
    for (int k=1; k<n.size();k++){
        theta[k]=asind((n[k-1]/n[k])*sind(theta[k-1]));
    }
    return theta;

}
vector<vector<complex<double>>> get_M(double n,double a,double l,double landa, double theta){

        complex<double> bi ((2*M_PI*a*l)/landa , (2*M_PI*n*l)/landa);
        bi=bi* cosd(theta);
        vector<vector<complex<double>>>  M_matrix {{0,0},{0,0}};


        M_matrix[0][0]=exp(bi);
        M_matrix[1][1]=exp(-bi);

        return M_matrix;

}
vector<vector<complex<double>>> get_D(double n1,double a1,double n2,double a2,double theta1, double theta2){
    vector<vector<complex<double>>>  D_matrix {{0,0},{0,0}};
    complex<double> n1_complex (n1,a1);
    complex<double> n2_complex (n2,a2);
    n1_complex=n1_complex*cosd(theta1);
    n2_complex=n2_complex*cosd(theta2);
    complex<double> r12,r21,t12,t21;

    r12=(n1_complex-n2_complex)/(n1_complex+n2_complex);
    r21=(n2_complex-n1_complex)/(n1_complex+n2_complex);
    t12=(2.0*n1_complex)/(n1_complex+n2_complex);
    t21=(2.0*n2_complex)/(n1_complex+n2_complex);

    D_matrix[0][0]=1.0/(t12);
    D_matrix[0][1]=(-r21)/(t12);
    D_matrix[1][0]=(r12)/(t12);
    D_matrix[1][1]=(t12*t21-r12*r21)/(t12);
    return D_matrix;

}
vector<vector<complex<double>>> matrix_product(vector<vector<complex<double>>> A,vector<vector<complex<double>>> B){
    vector<vector<complex<double>>>  C {{0,0},{0,0}};

    	    for(int i = 0; i < 2; ++i)
                for(int j = 0; j < 2; ++j)
                    for(int k = 0; k < 2; ++k)
                        C[i][j] += A[i][k] * B[k][j];
	        return C;
}

#endif // FUNC

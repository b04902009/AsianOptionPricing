#include <iostream>
#include<fstream>
#include<math.h>
using namespace std;

double S,X,H,t,v,r;
int n,k;
double delta_t ,r_bar ,u ,d ,p;

double AvgMax(int j, int i){
    double A = (1 - pow(u,j-i+1)) / (1-u);
    cout<<"A_MAX:"<<A<<endl;

    double B = pow(u,j-i)*d*((1-pow(d,i))/(1-d));
    cout<<"B_MAX:"<<B<<endl;

    return (S * A + S * B)/(j+1);
}
double AvgMin(double j, double i){

    double A = (1 - pow(d,i+1)) / (1-d);
    cout<<"A_Min:"<<A<<endl;

    double B = pow(d,i)*u*((1-pow(u,j-i))/(1-u));
    cout<< pow(d,i)<<endl;
    cout<<u<<endl;
    cout<<((1-pow(u,j-i))/(1-u))<<endl;
    cout<<"B_Min:"<<B<<endl;

    return (S * A + S * B)/(j+1);
}
double InterpoStates(int j, int i, int m){
     return ((k-m)/(double)k)*AvgMin(j,i) + (m/(double)k)*AvgMax(j,i);
}


int main(){
    ifstream fin("test.txt");
    fin>>S>>X>>H>>t>>v>>r>>n>>k;

    delta_t = t / n;
    r_bar = r*delta_t;
    u = 4;
    d = 1/u;
    p = (exp(r_bar)-d)/(u-d);

    int i=1;
    int j=1;

    cout<< AvgMax(1,1)<<endl;
    cout<< AvgMin(1,1)<<endl;

    return 0;
}

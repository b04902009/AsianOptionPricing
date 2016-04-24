#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;

// Input for Asian American Barrier Pricing
double S,X,H,t,v,r;
int n,k;

// BOPM and Black-Scholes Model parameter
double delta_t ,r_bar ,u ,d ,p;


double ***Create3DArray(int a, int b, int c) {
    double ***A = new double**[a];
    for (int i = 0 ; i < a ; i++) {
        A[i] = new double*[b];
        for (int j = 0 ; j < b ; j++) {
            A[i][j] = new double[c];
        }
    }
    return A;
}


double AvgMax(int j, int i){
    double A = (1 - pow(u,j-i+1)) / (1-u);
    double B = pow(u,j-i)*d*((1-pow(d,i))/(1-d));
    return (S * A + S * B)/(j+1);
}
double AvgMin(int j, int i){
    double A = (1 - pow(d,i+1)) / (1-d);
    double B = pow(d,i)*u*((1-pow(u,j-i))/(1-u));
    return (S * A + S * B)/(j+1);
}
double InterpoStates(int j, int i, int m){
     return ((k-m)/(double)k)*AvgMin(j,i) + (m/(double)k)*AvgMax(j,i);
}


double RunAvgAu(double cur_node, int j, int i){
    return ((j+1)*cur_node + S*pow(u,j+1-i)*pow(d,i)) / (j+2);
}
double RunAvgAd(double cur_node, int j, int i){
    return ((j+1)*cur_node + S*pow(u,j-i)*pow(d,i+1)) / (j+2);
}


int FindFootL(double A, int j, int i){
    double delta_max_min = ( AvgMax(j,i)-AvgMin(j,i) )/k;
    return (int) floor( (A - AvgMin(j,i))/ delta_max_min );
}
double FindInterpoXu(double Au, double*** asianMTree, int j, int i, int l){
    //if(asianMTree[j+1][i][l] - asianMTree[j+1][i][l+1] < 0.001)
    //    return 0;
    //else
    return (Au - asianMTree[j+1][i][l+1])/(asianMTree[j+1][i][l] - asianMTree[j+1][i][l+1]);
}
double FindInterpoXd(double Ad, double*** asianMTree, int j, int i, int l){
    //if(asianMTree[j+1][i+1][l] - asianMTree[j+1][i+1][l+1] < 0.001)
    //    return 0;
    //else
    return (Ad - asianMTree[j+1][i+1][l+1])/(asianMTree[j+1][i+1][l] - asianMTree[j+1][i+1][l+1]);
}
double InterpoCu(double*** CTree, double x, int j, int i, int l){
     return x*CTree[j+1][i][l] + (1-x)*CTree[j+1][i][l+1];
}
double InterpoCd(double*** CTree, double x, int j, int i, int l){
     return x*CTree[j+1][i+1][l] + (1-x)*CTree[j+1][i+1][l+1];
}


double AsianOptions(){
    delta_t = t / n;
    r_bar = r*delta_t;
    u = exp(v*sqrt(delta_t));
    d = 1/u;
    p = (exp(r_bar)-d)/(u-d);
    int node_num = n+1;

    double ***asianMTree = Create3DArray(node_num, node_num, k+1);
    double ***CTree = Create3DArray(node_num, node_num, k+1);

    for(int j=0; j<node_num; j++){
        for(int i=0; i<node_num; i++){
            //if(j==1 && i==1){
                //cout<<(1 - pow(u,j-i+1)) / (1-u)<<endl;
                //cout<<pow(u,j-i)*d*((1-pow(d,i))/(1-d))<<endl;
                //cout<< AvgMin(1,1)<<endl;
                //cout<< AvgMax(1,1)<<endl;
                //cout<<asianMTree[1][1][0]<<endl;
                //cout<<"$$$$$$$$"<<endl;
            //}
            for(int m=0; m<=k; m++){
                asianMTree[j][i][m] = InterpoStates(j, i, m);
            }
            //if(j==node_num-2)
            //cout<<"*****"<<endl;
        }
    }


    double Au, Ad, interpo_xu=0, interpo_xd=0, Cu, Cd, asianC, americanC;
    int foot_lu=0, foot_ld=0;
    for(int j=node_num-1; j>=0; j--){
        for(int i=0; i<=j; i++){
            //americanC = S*pow(u,j-i)*pow(d,i) - X;
            for(int m=0; m<=k; m++){
                Au = RunAvgAu(asianMTree[j][i][m], j, i);
                Ad = RunAvgAd(asianMTree[j][i][m], j, i);
                if(j==node_num-1){
                    Cu = fmax(Au - X, 0);
                    Cd = fmax(Ad - X, 0);
                }
                else{
                    //if(Au >= AvgMax(j,i)){
                    if(Au > AvgMax(j,i) || (Au-AvgMax(j,i))<0.0001 || (AvgMax(j,i)-Au)<0.0001){
                        foot_lu = k;
                        Cu = CTree[j+1][i][k];
                    }
                    else if(Au < AvgMin(j,i)){
                        foot_lu = 0;
                        Cu = CTree[j+1][i][0];
                    }
                    else{
                        foot_lu = FindFootL(Au, j, i);
                        interpo_xu = FindInterpoXu(Au, asianMTree, j, i, foot_lu);
                        Cu = InterpoCu(CTree, interpo_xu, j, i, foot_lu);
                    }

                    //if(Ad > AvgMax(j,i)){
                    if(Ad > AvgMax(j,i) || (Ad-AvgMax(j,i))<0.0001 || (AvgMax(j,i)-Ad)<0.0001){
                        foot_ld = k;
                        Cd = CTree[j+1][i+1][k];
                    }
                    else if(Ad < AvgMin(j,i)){
                        foot_ld = 0;
                        Cd = CTree[j+1][i+1][0];
                    }
                    else{
                        foot_ld = FindFootL(Ad, j, i);
                        interpo_xd = FindInterpoXd(Ad, asianMTree, j, i, foot_ld);
                        Cd = InterpoCd(CTree, interpo_xd, j, i, foot_ld);
                    }
                }

                asianC = (p*Cu + (1-p)*Cd)/exp(r_bar);

                //if(asianMTree[j][i][m] >= H){
                //    CTree[j][i][m] = 0;
                //}
                //else{
                CTree[j][i][m] = asianC;
                //}
                //if(j==node_num-2 && m==0){
                    //cout<<"###j==node_num-2###"<<endl;
                    //cout<<"Au:"<<Au<<", l:"<<foot_lu<<", x:"<<interpo_xu<<", Cu:"<<Cu<<endl;
                    //cout<<"Ad:"<<Ad<<", l:"<<foot_ld<<", x:"<<interpo_xd<<", Cu:"<<Cd<<endl;
                    //cout<<"asianC:"<<asianC<<endl;
                    //cout<<"######"<<endl;
                    //int temp2;
                    //cin>>temp2;
                //}
            }
            cout<<CTree[j][i][0]<<" ";
        }
        //int temp;
        //cin>>temp;
        cout<<endl;
        //cout<<"**********"<<endl;
    }

    return CTree[0][0][0];
}

int main(){
    ifstream fin("test.txt");
    fin>>S>>X>>H>>t>>v>>r>>n>>k;

    cout<<AsianOptions()<<endl;

    return 0;
}


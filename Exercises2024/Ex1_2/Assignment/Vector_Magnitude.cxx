// created by Hong-Gang on 13/11/24
#include <iostream>
#include <cmath>
using namespace std;

float My_Mag_Function(float,float);

int main(){
    float x = 2.3;
    float y = 4.5;
    cout<<"The magnitude of the vector is "<<sqrt(x*x+y*y)<<endl;
    cout<<"The result from my own function is "<<My_Mag_Function(x,y)<<endl;
    return 0;
}

float My_Mag_Function(float i, float j){
    return pow(i*i+j*j,0.5);
}
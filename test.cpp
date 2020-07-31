#include <iostream>
#include <vector>
#include <unistd.h>

using namespace std;

// [[Rcpp::export]]
void cpp_printer(){
    vector<string> vals(10,".. ");
    string a = "";
    unsigned int ct = 0;
    for(auto& v : vals){
        sleep(1);
        a += v;
        cout << a << ++ct << "\r" << flush;
    }
}
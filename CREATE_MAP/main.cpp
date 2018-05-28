#include <iostream>
#include <ctime>
#include <random>
#include <fstream>
#include <string>

const int Route_num = 1;
const int Map_x = 100;
const int Map_y = 100;
const int Map_z = 100;

void set_map(){
    srand((unsigned)time(NULL));
    std::string name("map.txt");
    std::ofstream map(name);

    for (int i = 0; i != Route_num; ++i){
        for (int ii = 0; ii != Map_x+1; ++ii){
            double x_tmp,y_tmp,z_tmp;
            x_tmp = ii;
            y_tmp = sin(x_tmp/10+10)*Map_y/2 + 50;
            //z_tmp = cos(x_tmp/8) * Map_z/2;

            map << x_tmp << '\t' << y_tmp << std::endl;//<< '\t' << z_tmp << std::endl;


        }
    }
}


int main() {
    set_map();
    return 0;
}
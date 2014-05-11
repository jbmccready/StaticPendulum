//#include "map_tools.h"
#include "pendulum_map.h"
#include "pendulum_system.h"
#include "Integrators/rk4.h"
#include "Integrators/ck45.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>
#include <thread>
#include <algorithm>
#include <QDomDocument>

typedef std::array< double , 4 > state_type;
int main()
{
    typedef ck45 integrator_type;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    pendulum_system mysystem;
    integrator_type myintegrator;
    pendulum_map<integrator_type> mymap;
//    mysystem.clear_attractors();
//    mysystem.add_attractor(0.5, 0.5);
//    mysystem.add_attractor(-3.0, 3.0);
//    mysystem.add_attractor(3.0, -3.0);
//    mysystem.add_attractor(-3.0, -3.0);
//    mysystem.add_attractor(1.0, 1.0);
//    mysystem.add_attractor(-1.0, 1.0);
//    mysystem.add_attractor(1.0, -1.0);
//    mysystem.add_attractor(-1.0, -1.0);
//    mysystem.set_all_attractor_strengths(3.0);
//    mysystem.g = 0.5;
//    mymap.set_converge_tol(0.5,0.1,10);
//    mymap.set_map(-10.0, 10.0, -10.0, 10.0, .01);
//    mymap.add_attractor_color(102, 51, 0);
//    mymap.add_attractor_color(120, 150, 250);
//    mymap.add_attractor_color(50, 200, 250);
//    mymap.add_attractor_color(120, 30, 0);
//    mymap.add_attractor_color(250, 150, 250);
//    mysystem.b = 0.5;
//    mymap.set_map(-10.0, 10.0, -10.0, 10.0, 0.025);
    QString count;
    QDomDocument mydoc("MapBatch");
    QDomElement root = mydoc.createElement("Maps");
    mydoc.appendChild(root);
    QDomElement map_element;
    for (int i = 0; i <= 0; i++) {
        if (i < 10) {
            count = "00" + QString::number(i);
        } else if (i < 100) {
            count = "0" + QString::number(i);
        } else {
            count = QString::number(i);
        }
//        mysystem.b = 0.1 + i*0.0008;
        map_element = mydoc.createElement("map" + count);
        root.appendChild(map_element);
        mymap.save_integrated_map(mysystem, myintegrator, count, map_element);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "\nTotal elapsed time: " << elapsed_seconds.count() << "s.\n";
    return 0;
}

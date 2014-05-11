#ifndef PENDULUM_SYSTEM_H
#define PENDULUM_SYSTEM_H


#include <cmath>
#include <array>
#include <vector>

typedef std::array< double , 4 > state_type;

struct attractor {
    double x; // x position
    double y; // y position
    double k; // Coulomb-like interaction k/r^2 force constant (N*m^2)
};

class pendulum_system {
public:
    double d = 0.05; // Distance from the end of the pendulum to the base plate (m)
    double m = 1.0; // Mass of the bob on the end of the pendulum (kg)
    double g = 9.8; // Acceleration due to gravity (m/s^2)
    double b = 0.2; // Linear drag coefficient (kg/s)
    double L = 10.0; // Length of the pendulum (m)
    std::vector<attractor> attractor_list; // Attractor position list

    // Default constructed magnet positions and strengths
    pendulum_system() {
        attractor_list.push_back(attractor{-0.5, sqrt(3.0)/2.0, 1.0});
        attractor_list.push_back(attractor{-0.5, -sqrt(3.0)/2.0, 1.0});
        attractor_list.push_back(attractor{1.0, 0.0, 1.0});
    }

    // Functor call for the system
    void operator()(const state_type &x, state_type &dxdt, const double /* t */) const; // no time dependence

    // methods for adding, removing and adjusting attractors
    void add_attractor(double x_position, double y_position, double attraction_strength);
    void set_attractor(int index, double x_position, double y_position, double attraction_strength);
    void set_all_attractor_strengths(double attraction_strength);
    void clear_attractors();
};

inline void pendulum_system::operator() (const state_type &x, state_type &dxdt, const double /* t */) const
{
    const double x_squared = x[0]*x[0];
    const double y_squared = x[1]*x[1];
    const double L_squared = L*L;
    const double norm_squared = x_squared + y_squared;

    const double g_value = -m*g/L * sqrt(1.0 - norm_squared/L_squared);

    double f_m_x = 0.0;
    double f_m_y = 0.0;

    const double a_value = (d+L-sqrt(L_squared-norm_squared));
    const double a_value_squared = a_value*a_value;

    double ax_value;
    double ay_value;
    double ax_value_squared;
    double ay_value_squared;
    double a_denom;
    for (auto attractor : attractor_list) {
        ax_value = x[0]-attractor.x;
        ay_value = x[1]-attractor.y;
        ax_value_squared = ax_value*ax_value;
        ay_value_squared = ay_value*ay_value;
        a_denom = -attractor.k/pow(ax_value_squared+ay_value_squared+a_value_squared,1.5);

        f_m_x += ax_value*a_denom;
        f_m_y += ay_value*a_denom;
    }

    dxdt[0] = x[2];
    dxdt[1] = x[3];
    dxdt[2] = (x[0]*g_value - b*x[2] + f_m_x) / m;
    dxdt[3] = (x[1]*g_value - b*x[3] + f_m_y) / m;
}

void pendulum_system::add_attractor(double x_position, double y_position, double attraction_strength = 1.0)
{
    attractor_list.push_back(attractor{x_position, y_position, attraction_strength});
}

void pendulum_system::set_attractor(int index, double x_position, double y_position, double attraction_strength)
{
    attractor_list[index].x = x_position;
    attractor_list[index].y = y_position;
    attractor_list[index].k = attraction_strength;
}

void pendulum_system::set_all_attractor_strengths(double attraction_strength)
{
    for (auto &attractor : attractor_list) {
        attractor.k = attraction_strength;
    }
}

void pendulum_system::clear_attractors()
{
    attractor_list.clear();
}

#endif // PENDULUM_SYSTEM_H

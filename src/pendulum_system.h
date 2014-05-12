#ifndef PENDULUM_SYSTEM_H
#define PENDULUM_SYSTEM_H


#include <cmath>
#include <array>
#include <vector>

/*! \mainpage
 *
 * \section intro_sec Introduction
 *
 * This is a modular code base for numerically integrating systems of differential equations across large sets of initial conditions.
 * Documentation is under construction.
 */

typedef std::array< double , 4 > state_type;

//! Pendulum function object that returns the derivative of the current state.

/*! The pendulum system is described by the following system of differential equations:
 *
 * \f$\vec{F}_{g}=-mg\frac{\sqrt{1-\frac{x^2+y^2}{L^2}}}{L}(x\hat{i}+y\hat{j})\f$ (force due to gravity)
 *
 * \f$\vec{F}_{a_n(x_n,y_n)}=\frac{-k[(x-x_n)\hat{i}+(y-y_n)\hat{j}]}{\left[(x-x_n)^2+(y-y_n)^2+\left(d+L-\sqrt{L^2-\left(x^2+y^2\right)}\right)\right]^{3/2}}\f$ (force due to attractor)
 *
 * \f$\vec{F}_{d}=-b(v_{x}\hat{i}+v_{y}\hat{j})\f$ (force due to dampening)
 *
 * \f$a_x=\frac{d^2x}{dt^2}=\left(F_{gx} + \sum\limits_{n=1}^{n=N}\vec{F}_{a_{nx}} + F_{dx}\right) \times \frac{1}{m}\f$ (acceleration in x direction)
 *
 * \f$a_y=\frac{d^2y}{dt^2}=\left(F_{gy} + \sum\limits_{n=1}^{n=N}\vec{F}_{a_{ny}} + F_{dy}\right) \times \frac{1}{m}\f$ (acceleration in y direction)
 *
 * Where m is the mass of the pendulum head, L is the length of the pendulum, g is the acceleration due to gravity, b is the dampening coefficient,
 * and d is the distance between the end of the pendulum at rest and the base plate of attractors.
 *
 * The system parameters are modified as public member variables while the attractors are modified by calling the respective member functions.
 *
 * The function is called using an overloaded () operator and returns through a reference parameter the derivative of the state passed in.
 */
class pendulum_system {
public:

    //! Container for an attractor, stores the position as an x-y coordinate, and an attractive force coefficient.
    struct attractor {
        double x; /*!< x coordinate position. */
        double y; /*!< y coordinate position. */
        double k; /*!< Attractive force coefficient where \f$F_{attractor} = \frac{-k}{x^2+y^2}\f$ */
    };

    double d = 0.05; /*!< Distance between the pendulum head at rest and the base plate. */
    double m = 1.0; /*!< Mass of the head of the pendulum. */
    double g = 9.8; /*!< Acceleration due to gravity. */
    double b = 0.2; /*!< Linear drag coefficient. */
    double L = 10.0; /*!< Length of the pendulum. */
    std::vector<attractor> attractor_list; /*!< List of attractors for the system. */

    //! Default constructor sets k = 1.0 for three attractors positioned at: \f$(-0.5, \sqrt{3}/2)\f$, \f$(-0.5, -\sqrt{3}/2)\f$, and \f$(1, 0)\f$.
    pendulum_system() {
        attractor_list.push_back(attractor{-0.5, sqrt(3.0)/2.0, 1.0});
        attractor_list.push_back(attractor{-0.5, -sqrt(3.0)/2.0, 1.0});
        attractor_list.push_back(attractor{1.0, 0.0, 1.0});
    }

    //! Function call that returns the derivative of the current state.
    void operator()(const state_type &x /*!< Current state input. */,
                    state_type &dxdt /*!< Derivative of the state, value modified by reference. */,
                    const double t /*!< Note: no time dependence. Parameter here to fit signature for integration.*/) const;

    //! Add an attractor at position (x_position, y_position) with attractive force coefficient attraction_strength.
    void add_attractor(double x_position, double y_position, double attraction_strength);

    //! Set new position and attraction strength for already existing attractor at an index.
    void set_attractor(int index, double x_position, double y_position, double attraction_strength);

    //! Set all the attractor strengths to the same value.
    void set_all_attractor_strengths(double attraction_strength);

    //! Clear all attractors.
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

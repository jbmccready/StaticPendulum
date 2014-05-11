#ifndef CK45_H
#define CK45_H
#include <array>
#include <algorithm>
// Cash and Karp Runge Kutta order 4/5 adapative step integrator
class ck45
{
public:
    ck45() {}

    // do_step attempts to perform one integration step, if successful returns 1 and updates the state, time and step size.
    // if unsuccessful returns 0 and adjusts the step size.
    template<typename system, typename state_type>
    int do_step (const system &dxdt, state_type &x, double &t, double &h) const;

    void set_tolerance(double relative_tolerance, double absolute_tolerance);
    void set_max_step_size(double max_step_size);
private:
    double m_rel_tol = 1e-6; // relative error tolerance for the adaptive integrator
    double m_abs_tol = 1e-6; // absolute error tolerance for the adaptive integrator
    double m_max_step_size = 0.1; // maximum stepsize

    // coefficients for method
    static constexpr double c[6] = {0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0};
    static constexpr double b_5th[6] = {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0};
    static constexpr double b_4th[6] = {2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0};
    static constexpr double b_diff[6] = {b_5th[0]-b_4th[0], b_5th[1]-b_4th[1], b_5th[2]-b_4th[2], b_5th[3]-b_4th[3], b_5th[4]-b_4th[4], b_5th[5]-b_4th[5]};
    static constexpr double a[6][5] = {
        {},
        {1.0/5.0},
        {3.0/40.0, 9.0/40.0},
        {3.0/10.0, -9.0/10.0, 6.0/5.0},
        {-11.0/54.0, 5.0/2.0, -70.0/27.0, 35.0/27.0},
        {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0}};
};

template<typename system, typename state_type>
int ck45::do_step (const system &dxdt, state_type &x, double &t, double &h) const
{
    const unsigned int state_size = x.size();
    std::array<state_type, 6> k;
    state_type temp_state; // used to store state for next k value and later used for 4th order solution
    dxdt(x, k[0], t);
    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = x[i]+h*a[1][0]*k[0][i];
    }
    dxdt(temp_state, k[1], t+c[1]*h);

    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = x[i]+h*(a[2][0]*k[0][i]+a[2][1]*k[1][i]);
    }
    dxdt(temp_state, k[2], t+c[2]*h);

    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = x[i]+h*(a[3][0]*k[0][i]+a[3][1]*k[1][i]+a[3][2]*k[2][i]);
    }
    dxdt(temp_state, k[3], t+c[3]*h);

    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = x[i]+h*(a[4][0]*k[0][i]+a[4][1]*k[1][i]+a[4][2]*k[2][i]+a[4][3]*k[3][i]);
    }
    dxdt(temp_state, k[4], t+c[4]*h);

    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = x[i]+h*(a[5][0]*k[0][i]+a[5][1]*k[1][i]+a[5][2]*k[2][i]+a[5][3]*k[3][i]+a[5][4]*k[4][i]);
    }
    dxdt(temp_state, k[5], t+c[5]*h);


    state_type order_5_solution;
    for (unsigned int i = 0; i < state_size; i++) {
        order_5_solution[i] = h*(b_5th[0]*k[0][i]+b_5th[1]*k[1][i]+b_5th[2]*k[2][i]+b_5th[3]*k[3][i]+b_5th[4]*k[4][i]+b_5th[5]*k[5][i]);
    }

    // difference between order 4 and 5, used for error check, reusing temp_state variable
    for (unsigned int i = 0; i < state_size; i++) {
        temp_state[i] = h*(b_diff[0]*k[0][i]+b_diff[1]*k[1][i]+b_diff[2]*k[2][i]+b_diff[3]*k[3][i]+b_diff[4]*k[4][i]+b_diff[5]*k[5][i]);
    }

    // boost odeint syle error step sizing method
    state_type error_val_list;
    for (unsigned int i = 0; i < state_size; i++) {
        error_val_list[i] = std::abs(temp_state[i]/(m_abs_tol + m_rel_tol * (x[i] + order_5_solution[i])));
    }


    double max_error_val = *(std::max_element(error_val_list.begin(), error_val_list.end()));

    if (max_error_val > 1.0) {
        // reject step and decrease step size
        h = h*std::max(0.9*std::pow(max_error_val, -0.25), 0.2);
        return 0;
    }
    else
    {
        if (max_error_val < 0.5) {
            // use step and increase step size
            t = t+h;
            for (unsigned int i = 0; i < state_size; i++) {
                x[i] = x[i] + order_5_solution[i];
            }
            h = std::min(h*std::min(0.9*std::pow(max_error_val, -0.20), 5.0), m_max_step_size);
            return 1;
        }
        else
        {
            // use step and keep same step size
            t = t+h;
            for (unsigned int i = 0; i < state_size; i++) {
                x[i] = x[i] + order_5_solution[i];
            }
            return 1;
        }
    }
}

void ck45::set_tolerance(double relative_tolerance, double absolute_tolerance)
{
    m_rel_tol = relative_tolerance;
    m_abs_tol = absolute_tolerance;
}

void ck45::set_max_step_size(double max_step_size)
{
    m_max_step_size = max_step_size;
}

#endif // CK45_H

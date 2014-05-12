#ifndef RK4_H
#define RK4_H
#include <array>

/*!
 * \brief Classic runge-kutta 4 method.
 *
 * Simple non-adaptive runge-kutta method, see Wikipedia: http://en.wikipedia.org/wiki/Runge-Kutta_methods#The_Runge.E2.80.93Kutta_method
 */
class rk4
{
public:
    rk4() {}

    //! Performs one step and returns 1.
    template<typename system, typename state_type>
    int do_step (const system &dxdt, state_type &x, double &t, const double h)
    {
        const unsigned int state_size = x.size();
        std::array<state_type, 4> k;
        state_type temp_state;
        dxdt(x, k[0], t);
        for (unsigned int i = 0; i < state_size; i++) {
            temp_state[i] = x[i]+0.5*h*k[0][i];
        }
        dxdt(temp_state, k[1], t+0.5*h);
        for (unsigned int i = 0; i < state_size; i++) {
            temp_state[i] = x[i]+0.5*h*k[1][i];
        }
        dxdt(temp_state, k[2], t+0.5*h);
        for (unsigned int i = 0; i < state_size; i++) {
            temp_state[i] = x[i]+h*k[2][i];
        }
        dxdt(temp_state, k[3], t+h);
        for (unsigned int i = 0; i < state_size; i++) {
            x[i] = x[i]+1.0/6.0*h*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i]);

        }
        t = t+h;
        return 1;
    }
};

#endif // RK4_H


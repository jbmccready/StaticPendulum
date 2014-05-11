#ifndef PENDULUM_MAP_H
#define PENDULUM_MAP_H
#include "pendulum_system.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <array>
#include <thread>
#include <chrono>
#include <QImage>
#include <QDomDocument>
#include <QString>
#include <QVector>

typedef std::array< double , 4 > state_type;

struct point_type
{
    state_type start_state;
    int converge_position = 255; // 255 reserved for points that do not converge
    double converge_time = 0.0;
    unsigned int step_count = 0;
};

typedef std::vector< std::vector<point_type> > map_type;
typedef std::vector< std::vector<point_type> >::iterator map_iter;


template <typename integrator_type>
class pendulum_map
{
public:
    pendulum_map();

    // parallel integrate and save colored map of convergence and grayscale map of integration time to convergence as png, file name is of the form position_map*filename*.png and time_map*filename*.png (without the asterisks)
    void save_integrated_map(pendulum_system &the_system, integrator_type &the_integrator, QString filename, QDomElement xml_element) const;

    // parallel integrate the map
    void parallel_integrate_map(const pendulum_system &the_system, const integrator_type &the_integrator, map_type &the_map) const;

    // integrate a single point
    void integrate_point(const integrator_type &the_integrator, const pendulum_system &the_system, point_type &the_point) const;

    // create a map of type map_type for the current x-y ranges and res
    map_type create_map_container() const;

    // integrate points and stop after converging to an attractor or the middle, stores relevent point integration information inside the points in the process
    void integrate_map(const integrator_type &the_integrator, const pendulum_system &the_system, map_iter first, map_iter last) const;

    // integration point with fixed integration time
    void fixed_integrate_point(const integrator_type &the_integrator, const pendulum_system &the_system, point_type &the_point) const;

    // integrate map with fixed integration time
    void fixed_integrate_map(const integrator_type &the_integrator, const pendulum_system &the_system, map_iter first, map_iter last) const;

    //parallel integrate the map with fixed integration time
    void fixed_parallel_integrate_map(const pendulum_system &the_system, const integrator_type &the_integrator, map_type &the_map) const;

    // basic property modifiers for the map, integration and colors for attractors NOTE: must have a color assigned for each attractor or the image will not form correctly
    void set_map(double x_start_position, double x_end_position, double y_start_position, double y_end_position, double resolution);
    void set_converge_tol(double position_tolerance, double mid_position_tolerance, double time_tolerance);
    void set_thread_count(unsigned int nthreads);
    void set_step_size(double step_size);
    void set_end_time(double end_time);
    void set_attractor_color(int index, int r, int g, int b);
    void add_attractor_color(int r, int g, int b);
    void clear_attractor_colors();
    void set_no_converge_color(int r, int g, int b);
    void set_mid_converge_color(int r, int g, int b);
private:
    double m_res = 0.05; // resolution of the map
    double m_xstart = -10.0; // start and end points for the map
    double m_ystart = -10.0;
    double m_xend = 10.0;
    double m_yend = 10.0;
    double m_tstart = 0.0; // starting integration time
    double m_tend = 20.0; // end integration time (only used for fixed integrations)
    double m_dt = 0.001; // starting integration step size
    unsigned int m_min_group = 1; // minimum group size for threading
    unsigned int m_nthreads = 32; // number of threads
    double m_pos_tol = 0.5; // position tolerance for checking magnet convergence
    double m_mid_tol = 0.1; // position tolerance for checking mid/gravity convergence
    double m_time_tol = 5.0; // time tolerance for checking convergence
    QVector<QRgb> attractor_colors; // index of colors to be assigned to the attractors
    QRgb no_converge_color = qRgb(255, 255, 255); // color for points that are outside bounds or do not converge to the middle or attractors
    QRgb mid_converge_color = qRgb(0, 0, 0); // color for points that converge to the middle
};

template <typename integrator_type>
pendulum_map<integrator_type>::pendulum_map()
{
    attractor_colors.push_back(qRgb(255, 140, 0));
    attractor_colors.push_back(qRgb(30, 144, 255));
    attractor_colors.push_back(qRgb(178, 34, 34));
}

template <typename integrator_type>
void pendulum_map<integrator_type>::parallel_integrate_map(const pendulum_system &the_system, const integrator_type &the_integrator, map_type &the_map) const
{
    // multithreaded integration of the map
    const unsigned int group = std::max(std::intptr_t(m_min_group), std::intptr_t((the_map.end()-the_map.begin())/m_nthreads));
    std::vector<std::thread> threads;
    threads.reserve(m_nthreads);
    auto it = the_map.begin();
    auto last = the_map.end();
    // function pointer to resolve compiler parsing UPDATE: Lambda cleaner and no pointer cleanup needed (compiled on gcc 4.8.1)
//    void (pendulum_map<integrator_type>::*fxn_ptr)(const integrator_type&, const pendulum_system&, map_iter, map_iter) const = &pendulum_map<integrator_type>::integrate_map;
    for (; it < last-group; it += group) {
//        threads.push_back(std::thread(fxn_ptr, this, std::ref(the_integrator), std::ref(the_system), it, std::min(it+group, last)));
        threads.push_back(std::thread([=,&the_integrator,&the_system]() {integrate_map(the_integrator, the_system, it, std::min(it+group, last));}));
    }
    integrate_map(the_integrator, the_system, it, the_map.end()); // left over chunks while we wait for other threads
    std::for_each(threads.begin(), threads.end(), [](std::thread& x){x.join();}); // wait for all threads to finish integrations
}

template <typename integrator_type>
void pendulum_map<integrator_type>::save_integrated_map (pendulum_system &the_system, integrator_type &the_integrator, QString filename, QDomElement xml_element) const
{
    // timer for computation time
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // dimensions of the map
    const int xdim = std::abs(std::round((m_xend-m_xstart)/m_res))+1;
    const int ydim = std::abs(std::round((m_yend-m_ystart)/m_res))+1;

    // blocks of memory for images
    uchar* position_solution_map = (uchar*)malloc(xdim * ydim);
    uchar* time_solution_map = (uchar*)malloc(xdim * ydim);

    //create map container and integrate the map
    map_type integration_map = create_map_container();
    parallel_integrate_map(the_system, the_integrator, integration_map);

    // collect general information about the map
    unsigned int buffer_index = 0;
    unsigned int total_count = 0;
    unsigned int mid_converge_count = 0;
    unsigned int outside_bounds_count = 0;
    double total_integration_time = 0.0;
    unsigned long long total_steps = 0;
    double max_time = 0;
    for (int j = ydim-1; j >= 0; j--) { // starting at upper left of map (ymax, xmin) to fill the memory with the correct orientation for QImage
        for (int i = 0; i < xdim; i++) {
            position_solution_map[buffer_index] = integration_map[i][j].converge_position;
            time_solution_map[buffer_index] = int(std::round(integration_map[i][j].converge_time));
            if (max_time < integration_map[i][j].converge_time) {
                max_time = integration_map[i][j].converge_time;
            }
            if (integration_map[i][j].converge_position == 255) {
                outside_bounds_count++;
            } else {
                total_count++;
                total_integration_time += integration_map[i][j].converge_time;
                total_steps += integration_map[i][j].step_count;
                if (integration_map[i][j].converge_position == 254) {
                    mid_converge_count++;
                }
            }
            buffer_index++;
        }
    }

    double avg_integration_time = total_integration_time/double(total_count);
    double avg_step_count = double(total_steps)/double(total_count);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start; // elapsed time for the process

    xml_element.setAttribute("points_integrated", total_count);
    xml_element.setAttribute("mid_converge_count", mid_converge_count);
    xml_element.setAttribute("points_outside_bounds", outside_bounds_count);
    xml_element.setAttribute("computation_time", elapsed_seconds.count());
    xml_element.setAttribute("avg_integration_time", avg_integration_time);
    xml_element.setAttribute("avg_number_of_steps", avg_step_count);
    xml_element.setAttribute("max_integration_time", max_time);
    std::cout << "\nTotal number of points: " << total_count << "\n";
    std::cout << "Points outside bounds: " << outside_bounds_count << "\n";
    std::cout << "Mid converge count: " << mid_converge_count << '\n';
    std::cout << "Average integration time: " << avg_integration_time << '\n';
    std::cout << "Average number of steps: " << avg_step_count << '\n';
    std::cout << "Max integration time: " << max_time << '\n';
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    QImage position_map_image(position_solution_map, xdim, ydim, xdim, QImage::Format_Indexed8);
    position_map_image.setColorTable(attractor_colors);
    position_map_image.setColor(254, mid_converge_color);
    position_map_image.setColor(255, no_converge_color);
    position_map_image.save("position_map" + filename + ".png");

    QImage time_map_image(time_solution_map, xdim, ydim, xdim, QImage::Format_Indexed8);
    for (unsigned int i = 0; i <= std::round(max_time); i++) {
        int scale_factor = std::floor(255.0/std::round(max_time));
        time_map_image.setColor(i, qRgb(255-i*scale_factor, 255-i*scale_factor, 255-i*scale_factor));
    }
    time_map_image.save("time_map" + filename + ".png");
}

template <typename integrator_type>
map_type pendulum_map<integrator_type>::create_map_container() const
{
    const int xdim = std::round((m_xend-m_xstart)/m_res)+1;
    const int ydim = std::round((m_yend-m_ystart)/m_res)+1;
    int xdim_factor = std::round(m_xstart/m_res); // create int multipliers to fill array to avoid floating math rounding error
    int ydim_factor = std::round(m_ystart/m_res);

    map_type integration_map(xdim, std::vector<point_type>(ydim));
    for (int i = 0; i < xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            integration_map[i][j].start_state[0] = double(xdim_factor)*m_res;
            integration_map[i][j].start_state[1] = double(ydim_factor)*m_res;
            integration_map[i][j].start_state[2] = 0.0;
            integration_map[i][j].start_state[3] = 0.0;
            ydim_factor++;
        }

        ydim_factor = std::round(m_ystart/m_res); // reset y value for next column
        xdim_factor++;
    }

    return integration_map;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::integrate_map(const integrator_type &the_integrator, const pendulum_system &the_system, map_iter first, map_iter last) const
{
    for (auto it = first; it != last; it++) {
        for(auto it_inside = it->begin(); it_inside != it->end(); it_inside++) {
            integrate_point(the_integrator, the_system, *it_inside);
        }
    }
}

template <typename integrator_type>
inline void pendulum_map<integrator_type>::fixed_integrate_point(const integrator_type &the_integrator, const pendulum_system &the_system, point_type &the_point) const
{
    double t = m_tstart;
    double h = m_dt;
    the_point.converge_time = m_tend-m_dt;
    unsigned int trial_count = 0;

    // check pendulum length boundary
    if (std::sqrt(std::pow(the_point.start_state[0], 2.0) + std::pow(the_point.start_state[1], 2.0)) < (the_system.L - 1e-10)) {
        const double absX = std::abs(the_point.start_state[0]);
        const double absY = std::abs(the_point.start_state[1]);
        // avoid undefined point at (0,0)
        if(absX > 1e-10 || absY > 1e-10) {
            // integration good to go, create local state for integration to keep start state
            state_type current_state = the_point.start_state;
            while (t < m_tend && trial_count < 1000000 ) {
                the_point.step_count += the_integrator.do_step(the_system, current_state, t, h);
                trial_count++;
            }

            for (unsigned int i = 0; i < the_system.attractor_list.size(); i++) {
                if ((the_system.attractor_list[i].x-m_pos_tol<current_state[0]) && (current_state[0]<the_system.attractor_list[i].x+m_pos_tol) && (the_system.attractor_list[i].y-m_pos_tol<current_state[1]) && (current_state[1]<the_system.attractor_list[i].y+m_pos_tol)) {
                    the_point.converge_position = i;
                    return;
                }
            }

            if ((0.0 - m_mid_tol < current_state[0]) && (current_state[0] < 0.0 + m_mid_tol) && (0.0 - m_mid_tol < current_state[1]) && (current_state[1] < 0.0 + m_mid_tol)) {
                the_point.converge_position = 254;
            }
        }
    }
    // ELSE: Position at (0,0) or outside of bounds or not found above attractor, undefined behavior for our pendulum system, leave converge position unset
}

template <typename integrator_type>
void pendulum_map<integrator_type>::fixed_integrate_map(const integrator_type &the_integrator, const pendulum_system &the_system, map_iter first, map_iter last) const
{
    for (auto it = first; it != last; it++) {
        for(auto it_inside = it->begin(); it_inside != it->end(); it_inside++) {
            fixed_integrate_point(the_integrator, the_system, *it_inside);
        }
    }
}

template <typename integrator_type>
void pendulum_map<integrator_type>::fixed_parallel_integrate_map(const pendulum_system &the_system, const integrator_type &the_integrator, map_type &the_map) const
{
    const unsigned int group = std::max(std::intptr_t(m_min_group), std::intptr_t((the_map.end()-the_map.begin())/m_nthreads));
    std::vector<std::thread> threads;
    threads.reserve(m_nthreads);
    auto it = the_map.begin();
    auto last = the_map.end();
    for (; it < last-group; it += group) {
        threads.push_back(std::thread([=,&the_integrator,&the_system]() {fixed_integrate_map(the_integrator, the_system, it, std::min(it+group, last));}));
    }
    fixed_integrate_map(the_integrator, the_system, it, the_map.end());
    std::for_each(threads.begin(), threads.end(), [](std::thread& x){x.join();});
}

template <typename integrator_type>
inline void pendulum_map<integrator_type>::integrate_point(const integrator_type &the_integrator, const pendulum_system &the_system, point_type &the_point) const
{
    double t = m_tstart;
    double h = m_dt;
    unsigned int trial_count = 0;

    // check pendulum length boundary
    if (std::sqrt(std::pow(the_point.start_state[0], 2.0) + std::pow(the_point.start_state[1], 2.0)) < (the_system.L - 1e-10)) {
        const double absX = std::abs(the_point.start_state[0]);
        const double absY = std::abs(the_point.start_state[1]);
        // avoid undefined point at (0,0)
        if(absX > 1e-10 || absY > 1e-10) {
            // integration good to go, create local state for integration to keep start state and hopefully optimize memory rather than calling a member variable every time
            state_type current_state = the_point.start_state;
            bool converged = false;
            std::vector<double> magnet_time;
            bool near_magnet = false;
            unsigned int current_magnet = 0;
            while (!converged && t < 1000 && trial_count < 1000000 ) {
                the_point.step_count += the_integrator.do_step(the_system, current_state, t, h);
                trial_count++;
                near_magnet = false;

                // check if pendulum head near an attractor and if it's been near for long enough time to consider converged
                for (unsigned int i = 0; i < the_system.attractor_list.size(); i++) {
                    if ((the_system.attractor_list[i].x-m_pos_tol<current_state[0]) && (current_state[0]<the_system.attractor_list[i].x+m_pos_tol) && (the_system.attractor_list[i].y-m_pos_tol<current_state[1]) && (current_state[1]<the_system.attractor_list[i].y+m_pos_tol)) {
                        near_magnet = true;
                        if (current_magnet == i) {
                            magnet_time.push_back(t);
                            if (magnet_time.back() - magnet_time.front() >= m_time_tol) {
                                magnet_time.clear();
                                the_point.converge_time = t;
                                the_point.converge_position = i;
                                converged = true;
                            }
                        } else {
                            magnet_time.clear();
                            current_magnet = i;
                            magnet_time.push_back(t);
                        }
                        break;
                    }
                }
                if (!near_magnet) {
                    // check if pendulum head near middle
                    if ((0.0 - m_mid_tol < current_state[0]) && (current_state[0] < 0.0 + m_mid_tol) && (0.0 - m_mid_tol < current_state[1]) && (current_state[1] < 0.0 + m_mid_tol)) {
                        near_magnet = true;
                        if (current_magnet == 254) {
                            magnet_time.push_back(t);
                            if (magnet_time.back() - magnet_time.front() >= m_time_tol) {
                                magnet_time.clear();
                                the_point.converge_time = t;
                                the_point.converge_position = 254;
                                converged = true;
                            }
                        } else {
                            magnet_time.clear();
                            current_magnet = 254;
                            magnet_time.push_back(t);
                        }
                    } else {
                        // not near middle or any attractor
                        magnet_time.clear();
                    }
                }
            }
        }
    }
    // ELSE: Position at (0,0) or outside of bounds, undefined behavior for our pendulum system, leave converge position unset
}



template <typename integrator_type>
void pendulum_map<integrator_type>::set_map(double x_start_position, double x_end_position, double y_start_position, double y_end_position, double resolution)
{
    m_xstart = x_start_position;
    m_xend = x_end_position;
    m_ystart = y_start_position;
    m_yend = y_end_position;
    m_res = resolution;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_converge_tol(double position_tolerance, double mid_position_tolerance, double time_tolerance)
{
    m_pos_tol = position_tolerance;
    m_mid_tol = mid_position_tolerance;
    m_time_tol = time_tolerance;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_thread_count(unsigned int nthreads)
{
    m_nthreads = nthreads;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_step_size(double step_size)
{
    m_dt = step_size;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_end_time(double end_time)
{
    m_tend = end_time;
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_attractor_color(int index, int r, int g, int b)
{
    attractor_colors[index] = qRgb(r, g, b);
}

template <typename integrator_type>
void pendulum_map<integrator_type>::add_attractor_color(int r, int g, int b)
{
    attractor_colors.push_back(qRgb(r, g, b));
}

template <typename integrator_type>
void pendulum_map<integrator_type>::clear_attractor_colors()
{
    attractor_colors.clear();
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_no_converge_color(int r, int g, int b)
{
    no_converge_color = qRgb(r, g, b);
}

template <typename integrator_type>
void pendulum_map<integrator_type>::set_mid_converge_color(int r, int g, int b)
{
    mid_converge_color = qRgb(r, g, b);
}

#endif // PENDULUM_MAP_H

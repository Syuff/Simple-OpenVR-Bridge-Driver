#pragma once

#include <stddef.h>

#include <array>
#include <cstdio>

class Observer_class {

    public:
        enum {
            N_ORDER  = 4, // order of model
            N_DOF  = 7, // degrees of freedom
            N_SENSORS = 3, // number of sensors
            POSITION = 0, // first order
            VELOCTY = 1, // second order
            ACCELERATION = 2, // third order
        };
        
        //output variable containing POSITION,VELOCTY, and ACCELERATION data
        std::array<std::array<double, N_DOF>, N_ORDER> output_states;      
        
        //reset states
        void reset_states();
        
        //add POSITION,VELOCTY, or ACCELERATION sensor
        void add_sensor(std::array<double,N_DOF> (&sensor),const int sensor_type);
        
        // update states from previous states and sensor data
        void update_states(const double Ts);

        //update output_states with time
        void output(const double Ts);
        
    private:

        template <class T, size_t ROW, size_t COL>
        using Matrix = std::array<std::array<T, COL>, ROW>;

        Matrix<double,N_ORDER,N_DOF> states = {{0}};
        Matrix<double,N_SENSORS,N_DOF> last_sensor = {0};  
        std::array<double,N_ORDER> R = {0.03,0.03,0.3, 0.03};
        std::array<int,N_SENSORS> sensors = {0};
        Matrix<double,N_ORDER,N_ORDER> P;
        Matrix<double,N_ORDER,N_ORDER> Q = 
	    {{
            {0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00},
            {0.000E+00, 1.058E-04, 3.203E-03, 2.405E-05},
            {0.000E+00, 3.203E-03, 9.620E-02, 4.810E-05},
            {0.000E+00, 2.405E-05, 4.810E-05, 9.620E-05}
        }};

        inline static const Matrix<double,N_ORDER,N_ORDER> A =
        {{
            {0.000E+00, 1.000E+00, 0.000E+00, 0.000E+00},
            {0.000E+00, 0.000E+00, 1.000E+00, 1.000E+00},
            {0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00},
            {0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00}
        }};
        
        inline static const Matrix<double,N_ORDER,N_ORDER> C =
        {{
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        }};
        
                

        
        Matrix<double,N_ORDER,N_ORDER> Matrix_add(const Matrix<double,N_ORDER,N_ORDER> (&input_1), const Matrix<double,N_ORDER,N_ORDER> (&input_2));
        Matrix<double,N_ORDER,N_ORDER> Matrix_mult(const Matrix<double,N_ORDER,N_ORDER> (&input_1), const Matrix<double,N_ORDER,N_ORDER> (&input_2));
        Matrix<double,N_ORDER,N_ORDER> Matrix_transpose(const Matrix<double,N_ORDER,N_ORDER> (&input));
        Matrix<double,N_ORDER,N_ORDER> inverse_S(const Matrix<double,N_ORDER,N_ORDER> (&input));
        int eye(const int i,const int j);
	
};



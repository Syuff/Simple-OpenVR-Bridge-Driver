#include "Observer.hpp"


void Observer_class::reset_states(){
    states = {{0}};
    last_sensor = {{0}};
    output_states = {{0}};
    P =     
    {{
        {1,     0,    0,    0},
        {0,     0,    0,    0},
        {0,     0,    0,    0},
        {0,     0,    0,    0}
    }};
    sensors = {0};

};

void Observer_class::add_sensor(std::array<double,N_DOF> (&sensor),const int sensor_type){
    std::copy(std::begin(sensor), std::end(sensor), last_sensor[sensor_type].begin());
    sensors[sensor_type] = 1;
};

void Observer_class::output(const double Ts){ 
    //std::array<double,N_DOF> temp_values = {0};
    Matrix<double,N_ORDER,N_DOF> next_states = {{0}};
    for (int k = 0; k < N_DOF; ++k) {
        for (int j = 0; j < N_ORDER; ++j) {
            for (int i = 0; i < N_ORDER; ++i) {
                next_states[j][k] += ((eye(j,i) + Ts * A[j][i]) * states[i][k]);
            }
        }
    }
    output_states = next_states;
};

void Observer_class::update_states(const double Ts){
    Matrix<double,N_ORDER,N_DOF> next_states = {{0}};
    Matrix<double,N_ORDER,N_ORDER> P_bar;
    Matrix<double,N_ORDER,N_ORDER> P_next = {{0}};
    Matrix<double,N_ORDER,N_ORDER> S = {{0}};
    Matrix<double,N_ORDER,N_ORDER> S_inv = {{0}};
    Matrix<double,N_ORDER,N_ORDER> K = {{0}};
    Matrix<double,N_ORDER,N_ORDER> A_d = {{0}};
    Matrix<double,N_ORDER,N_ORDER> C_i_T, A_d_T;
    Matrix<double,N_ORDER,N_ORDER> C_i = {{0}};
    Matrix<double,N_SENSORS,N_DOF> Sensor_diff;
    Matrix<double,N_ORDER,N_ORDER> Q_d;
    
    for (int j = 0; j < N_ORDER; ++j) {
        for (int i = 0; i < N_ORDER; ++i) {
            A_d[j][i] = eye(j,i) + Ts*A[j][i];
        }
    }
    
    A_d_T = Matrix_transpose(A_d);
    Q_d = Matrix_mult(Matrix_mult(A_d,Q),A_d_T);
    P_bar = Matrix_add(Matrix_mult(Matrix_mult(A_d,P),A_d_T),Q_d);
    
    for (int j = 0; j < N_ORDER; ++j) {
        if (sensors[j]){
                C_i[j][j] = 1;
        }
    }
    C_i_T = Matrix_transpose(C_i);
    S = Matrix_mult(Matrix_mult(C_i,P_bar),C_i_T);
    for (int j = 0; j < N_ORDER; ++j) {
        
        S[j][j] += R[j];
    }
    
    S_inv = inverse_S(S);

    K = Matrix_mult(Matrix_mult(P_bar,C_i_T),S_inv);
    
    for (int j = 0; j < N_ORDER; ++j){
        P_next[j][j] = 1;
    }
    for (int j = 0; j < N_ORDER; ++j) {
        for (int i = 0; i < N_ORDER; ++i) {
            for (int k = 0; k < N_ORDER; ++k){
                P_next[j][i] += - K[j][k] * C_i[k][i];
            }
        }
    }

    P = Matrix_mult(P_next,P_bar);
    
    Sensor_diff = last_sensor;
    last_sensor = {{0}};
    sensors = {0};
    
    for (int k = 0; k < N_DOF; ++k) {
        for (int j = 0; j < N_ORDER; ++j) {
            for (int i = 0; i < N_ORDER; ++i) {
                next_states[j][k] += A_d[j][i]*states[i][k];
            }
        }
    };

    for (int k = 0; k < N_DOF; ++k) {
        for (int j = 0; j < N_ORDER; ++j) {
            for (int i = 0; i < N_ORDER; ++i) {
                 Sensor_diff[j][k] += -C_i[j][i]*next_states[i][k];
            }
        }
    } 
    
    for (int k = 0; k < N_DOF; ++k) {
        for (int j = 0; j < N_ORDER; ++j) {
            for (int i = 0; i < N_ORDER; ++i) {
                next_states[j][k] += K[j][i]*Sensor_diff[i][k];
            }
        }
    };
    
    states = next_states;
};

int Observer_class::eye(const int i,const int j){
    if (i == j){
        return 1;
    }
    else{
        return 0;
    }
};

Observer_class::Matrix<double,Observer_class::N_ORDER,Observer_class::N_ORDER> Observer_class::Matrix_add(const Matrix<double,N_ORDER,N_ORDER> (&input_1), const Matrix<double,N_ORDER,N_ORDER> (&input_2)){
    Matrix<double,N_ORDER,N_ORDER> output;
    for (int j = 0; j < N_ORDER; ++j) {
        for (int i = 0; i < N_ORDER; ++i) {
            output[j][i] = input_1[j][i]+input_2[j][i];
        }
    }
    return output;
};

Observer_class::Matrix<double,Observer_class::N_ORDER,Observer_class::N_ORDER> Observer_class::Matrix_mult(const Matrix<double,N_ORDER,N_ORDER> (&input_1), const Matrix<double,N_ORDER,N_ORDER> (&input_2)){
    Matrix<double,N_ORDER,N_ORDER> output = {{0}};
    for (int j = 0; j < N_ORDER; ++j) {
        for (int i = 0; i < N_ORDER; ++i) {
            for (int k = 0; k < N_ORDER; ++k) { 
                output[j][i] += input_1[j][k]*input_2[k][i];
            }
        }
    }
    return output;
};

Observer_class::Matrix<double,Observer_class::N_ORDER,Observer_class::N_ORDER> Observer_class::Matrix_transpose(const Matrix<double,N_ORDER,N_ORDER> (&input)){
    Matrix<double,N_ORDER,N_ORDER> output = {{0}};
    for (int j = 0; j < N_ORDER; ++j) {
        for (int i = 0; i < N_ORDER; ++i) {
            output[j][i] = input[i][j];
        }
    }
    return output;
};
Observer_class::Matrix<double,Observer_class::N_ORDER,Observer_class::N_ORDER> Observer_class::inverse_S(const Matrix<double,N_ORDER,N_ORDER> (&input)){
    // Choose optimal inversion method depending on number of sensors.
    // diagonal inversion for num_sensors <= 1
    // otherwise solve determinant and lower triangle inversion
    Matrix<double,N_ORDER,N_ORDER> I = {{0}};
    int num_sensors = 0;
    double determinant = 0;
    
    for (int i = 0; i < N_SENSORS; ++i) {
        if (sensors[i]){
            ++num_sensors;
        }
    }
    if (num_sensors <= 1){
        I = {{0}};
        for (int j = 0; j < N_ORDER; ++j) {
            I[j][j] = 1/input[j][j];
        }
    }
    else{
        for(int i=0;i<3;i++){
            determinant = determinant + (input[0][i]*(input[1][(i+1)%3]*input[2][(i+2)%3] - input[1][(i+2)%3]*input[2][(i+1)%3]));
            }
        for(int i=0;i<3;i++){
            for(int j=0;j<=i;j++){
                I[i][j] = ((input[(i+1)%3][(j+1)%3] * input[(i+2)%3][(j+2)%3]) - (input[(i+1)%3][(j+2)%3]*input[(i+2)%3][(j+1)%3]))/ determinant;
                if (i!=j){
                    I[j][i] = I[i][j];
                }
            }
        }
    }
    return I;
};

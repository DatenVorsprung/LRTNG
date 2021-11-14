#!/bin/bash

# Brusselator
BENCHMARK=bruss
TIME_STEP=0.01
TIME_HORIZON=9
./LRTNG --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --time_step $TIME_STEP

# Dubins Car
BENCHMARK=dubins
TIME_STEP=0.00125
TIME_HORIZON=0.05
./LRTNG --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --time_step $TIME_STEP

# Cartpole with CTRNN DNN Controller
TIME_HORIZON=1
TIME_STEP=1e-5
BENCHMARK=cartpoleCTRNN
./LRTNG --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --time_step $TIME_STEP

# Cartpole with LTC DNN Controller
TIME_HORIZON=0.35
TIME_STEP=1e-6
BENCHMARK=cartpoleLTC_RK
./LRTNG --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --time_step $TIME_STEP

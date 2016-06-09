#!/bin/bash

dt=0.001
x=1
v=0
job=100
pt=1
c=1

printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.001.dat
dt=0.002
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.002.dat
dt=0.003
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.003.dat
dt=0.004
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.004.dat
dt=0.005
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.005.dat
dt=0.006
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.006.dat
dt=0.007
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.007.dat
dt=0.008
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.008.dat
dt=0.009
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.009.dat
dt=0.01
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.01.dat
dt=0.02
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.02.dat
dt=0.03
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.03.dat
dt=0.04
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.04.dat
dt=0.05
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.05.dat
dt=0.06
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.06.dat
dt=0.07
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.07.dat
dt=0.08
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.08.dat
dt=0.09
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.09.dat
dt=0.1
printf "$dt\n$x\n$v\n$job\n$pt" |./init_holf  
./holf > holf_dt=0.1.dat

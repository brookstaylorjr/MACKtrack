#!/bin/sh

mkdir -p $HOME/Documents/MATLAB/
touch startup.m $HOME/Documents/MATLAB/

stringOne=$PWD
stringTwo="addpath(genpath('"
stringThree="'))"

echo ${stringTwo}${stringOne}${stringThree} >> $HOME/Documents/MATLAB/startup.m


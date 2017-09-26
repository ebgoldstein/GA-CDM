clear all
close all
clc

%compile all the error/demon reg algorithms
cd demon_reg
compile_c_files
cd ..

%Given a:
%   -range of Beta condiitons: 0.001 10.0
%   -range of Hv condiitons: 0.001 1.0
%   -range of wind conditions: 0.1 0.4
%   -range of iterations: 100 1000

%Run a GA to optimize the parameters for CDM:

%%%%%%%%%%%%
%first move the grids
copyfile('SfMtogrid/April.dat','CDM/')
copyfile('SfMtogrid/AprilV.dat','CDM/')
copyfile('SfMtogrid/AprilVy.dat','CDM/')
copyfile('SfMtogrid/September.dat','CDM/')
copyfile('SfMtogrid/SeptemberV.dat','CDM/')
copyfile('SfMtogrid/November.dat','CDM/')
copyfile('SfMtogrid/NovemberV.dat','CDM/')

%%%%%%%%%%%%
%Then set up the bounds
bounds(1,:) = [0.001 10.0];
bounds(2,:) = [0.001 1.0];
bounds(3,:) = [0.2 0.4];
bounds(4,:) = [1 100];

LB=bounds(:,1);
UB=bounds(:,2);

POPSIZE=10;
MAXGEN=10;

options = optimoptions('ga');
options = optimoptions(options,'PopulationSize', POPSIZE);
options = optimoptions(options,'MaxGenerations', MAXGEN);
%options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @gaplotbestf });

%Now lets optimize
%This may take some time...
[x,fval,exitflag,output,population,score] = ga(@CDMfunction,4,[],[],[],[],LB,UB,[],[],options);

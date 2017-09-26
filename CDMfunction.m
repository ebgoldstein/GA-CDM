function [ val ] = CDMfunction(sol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% val - the fittness of this individual
% sol - the individual, returned to allow for Lamarckian evolution
% options - [current_generation]

cd CDM
%1. Build a parameter file (cdm_params.m)
%right now, param.par needs to be in this directory
%example for debugging:
%cdm_params( 1.0,0.001,0.2,1000)


%solution 4 is 1 to 10, but needs to be rounded to integer value and then
%multiplied by 100;
rsol=round(sol(4))*100;

%%%%%%%%%%%%

[SE]=cdm_params( sol(1),sol(2),sol(3),rsol);


%2. Run the model
%right now, you needD-M to be in the correct directory
system('./Dune param.par');
% system('./Dune params.par > output &'); %if you want to suppress output
% and store model info to output file. 'Try-Catch loop might be needed'

%3a. Pull the last topo file and assess error; GA finds a minima

LastSave=SE*100;
Model = importdata(sprintf('DATA/h.%05d.dat',LastSave));

%M=importdata(['DATA/h.01000.dat']);
D=importdata(['November.dat']);

%this is only over a region of the model
Lveg=55;
behindforedune=100;
numcellseval=((behindforedune-Lveg+1)*64);

%%%%%%%%Bosboom et al 2014 MSSES method%%%%%%%%

% I=importdata(['April.dat']);
% %df=M-Model;imagesc(df);colorbar
%
% delzp=Model(Lveg:behindforedune,:)-I(Lveg:behindforedune,:);
% delzo=D(Lveg:behindforedune,:)-I(Lveg:behindforedune,:);
%
% num=sum(sum((delzp-delzo).^2))/numcellseval;
% den=sum(sum((delzo).^2))/numcellseval;
%
% MSESS=1-(num/den);
% val=MSESS;

%%%%%%%%%%MSE%%%%%%%%%%%%%%%%
% diff=Model(Lveg:behindforedune,:)-D(Lveg:behindforedune,:);
% MSE=(sum(sum(diff.^2)))/((behindforedune-Lveg+1)*64);
%
% %%%%%%%%%RMSE%%%%%%%%%%%%%%%%
% RMSE=sqrt(MSE);
%
% % might want to normalize it by mean or range:
% % https://en.wikipedia.org/wiki/Root-mean-square_deviation
%
% val=RMSE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For the lines between the '%%%%%%%%%%%%%%%%DEMON%%%%%%%%%%%%%%%%' lines below, 
%this copyright applies (because the code is a modification:
%Copyright (c) 2009, Dirk-Jan Kroon 
%All rights reserved.

%Redistribution and use in source and binary forms, with or without 
%modification, are permitted provided that the following conditions are 
%met:

%* Redistributions of source code must retain the above copyright 
%notice, this list of conditions and the following disclaimer. 
%* Redistributions in binary form must reproduce the above copyright 
%notice, this list of conditions and the following disclaimer in 
%the documentation and/or other materials provided with the distribution

%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%DEMON%%%%%%%%%%%%%%%%
cd ..
cd demon_reg/functions_nonrigid

% Read two images
I1=Model(Lveg:behindforedune,:);  
I2=D(Lveg:behindforedune,:);

% Set static and moving image
S=I2; M=I1;
% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);

% The transformation fields
Tx=zeros(size(M)); Ty=zeros(size(M));

[Sy,Sx] = gradient(S);
for itt=1:200
    % Difference image between moving and static image
    Idiff=M-S;
    
    % Default demon force, (Thirion 1998)
    Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
    Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);
    
    % When divided by zero
    Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
    
    % Smooth the transformation field
    Uxs=3*imfilter(Ux,Hsmooth);
    Uys=3*imfilter(Uy,Hsmooth);
    
    % Add the new transformation field to the total transformation field.
    Tx=Tx+Uxs;
    Ty=Ty+Uys;
    M=movepixels(I1,Tx,Ty);
end
cd ..
p0=I1;     %predicted prior to warp
p1=M;   %deformed prediction
obs=I2; %observation
D=sqrt((Tx.^2)+(Ty.^2));    %displacement magnitude
Dmax=5;
littledelta=D/Dmax;
littledelta(littledelta>1)=1; %bosboom and reniers EQN5

SE0=(p0-obs).^2;    %locally squared error prior to warp
SE1=(p1-obs).^2;    %locally squared error after to warp

SEW=(SE1) + (littledelta.*(SE0-SE1)); %Bosboom adn Reniers EQN4

%get the number of cells in the array
[m,n] = size(SEW);
SEWsize=m*n;

RMSEW=sqrt((sum(sum(SEW)))/SEWsize);    %eqn 3 bosboom and reniers

val=RMSEW;
%3b. Alternatively, Pull the last veg file and assess error (MSE)
%VD=importdata(['SeptV.dat']);

cd ..
delete CDM/DATA/*.dat
%%%%%%%%%%%%%%%%DEMON%%%%%%%%%%%%%%%%'

end


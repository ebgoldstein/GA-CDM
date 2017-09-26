function [SE]=cdm_params( Beta,Hv,u, NT)
%CDM_PARAMS Write a CDM parameter file
%   Detailed explanation goes here
    %fid = fopen([project.directory, 'param.par'], 'w');
    %fid = fopen('/AGU2016simpleLatVeg/param.par','w');
    fid = fopen('param.par','w');

    
    fprintf(fid, '%s\n', ['NX = 128']);
    fprintf(fid, '%s\n', ['NY = 64']);
    fprintf(fid, '%s\n', ['dx = 1']);
    fprintf(fid, '%s\n', ['dt_max = 1210']); %time step (constant, max is not important) (sec)
    %fprintf(fid, '%s\n', ['Nt = 10000']); %number of iterations (such that total time is Nt*dt_max)
    fprintf(fid, '%s\n', ['Nt = ', num2str(NT)]);
    
    %in this case, total time is 12100000, or 1.21e+7, which is 140 days (time interval between SfM captures)
    
    fprintf(fid, '%s\n', ['Nt0 = 0']);
    
    %fprintf(fid, '%s\n', ['save.every = 1']);
    %fprintf(fid, '%s\n', ['save.every = 100']);
    
    %'save.every specific number of timesteps... need to make sure that;
    %number of saved files are 100, so Nt/save.every...
    SE=NT/100;
    fprintf(fid, '%s\n', ['save.every = ', num2str(SE)]);
    
    
    
    fprintf(fid, '%s\n', ['save.x-line = 0']);
    fprintf(fid, '%s\n', ['calc.x_periodic = 0']);
    fprintf(fid, '%s\n', ['calc.y_periodic = 1']);
    fprintf(fid, '%s\n', ['calc.shift_back = 0']);
    
    
    fprintf(fid, '%s\n', ['save.dir = ./DATA']);
    
    fprintf(fid, '%s\n', ['influx = const']);
    fprintf(fid, '%s\n', ['q_in = 0']);
    fprintf(fid, '%s\n', ['wind = const']);
    
    fprintf(fid, '%s\n', ['constwind.u = ', num2str(u)]);
   
    fprintf(fid, '%s\n', ['wind.fraction = 1']);
    fprintf(fid, '%s\n', ['veget.calc = 1']);
    fprintf(fid, '%s\n', ['veget.xmin = 55']);
    fprintf(fid, '%s\n', ['veget.sigma = 1.5']);
    fprintf(fid, '%s\n', ['veget.beta = 150']);
    fprintf(fid, '%s\n', ['veget.m = 0.16']);
    fprintf(fid, '%s\n', ['veget.season.t = 0']);
    
    fprintf(fid, '%s\n', ['veget.Vlateral.factor = ', num2str(Beta)]);
    fprintf(fid, '%s\n', ['veget.Tveg = ', num2str(Hv)]);
    
    fprintf(fid, '%s\n', ['beach.tau_t_L = 0.05']);
    fprintf(fid, '%s\n', ['shore.MHWL = 0']);
    fprintf(fid, '%s\n', ['veget.Init-Surf = init_h']);
     
    fprintf(fid, '%s\n', ['veget.init_h.file = AprilV.dat ']);
    
    
    fprintf(fid, '%s\n', ['veget.init_h.file_aux = AprilVy.dat']);
    
    fprintf(fid, '%s\n', ['Init-Surf = init_h']);
    fprintf(fid, '%s\n', ['init_h.file= April.dat']);
    
        %0 = save, 1 = don't save
    fprintf(fid, '%s\n', ['dontsave.veget = 0']);
    fprintf(fid, '%s\n', ['dontsave.u = 1']);
    fprintf(fid, '%s\n', ['dontsave.flux = 1']);
    fprintf(fid, '%s\n', ['dontsave.flux_s = 1']);
    fprintf(fid, '%s\n', ['dontsave.shear = 1']);
    fprintf(fid, '%s\n', ['dontsave.shear_pert = 1']);
    fprintf(fid, '%s\n', ['dontsave.stall = 1']);
    fprintf(fid, '%s\n', ['dontsave.rho = 1']);
    fprintf(fid, '%s\n', ['dontsave.h_deposit = 1']);
    fprintf(fid, '%s\n', ['dontsave.h_nonerod = 1']);
    fprintf(fid, '%s\n', ['dontsave.h_sep = 1']);
    fprintf(fid, '%s\n', ['dontsave.dhdt = 1']);
    fclose(fid); 
end


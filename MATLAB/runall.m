function runall()
    disp('Table 2, Column 2')
    [runtime1,~] = enzymekinetics_2D_IFETDRDP(1/100,99);
    fprintf('Time %f\n',runtime1);
    
    [runtime2,~] = enzymekinetics_2D_IFETDRDP(1/200,199);
    fprintf('Time %f\n',runtime2);
    
    [runtime3,~] = enzymekinetics_2D_IFETDRDP(1/400,399);
    fprintf('Time %f\n',runtime3);
    
    [runtime4,~] = enzymekinetics_2D_IFETDRDP(1/800,799);
    fprintf('Time %f\n',runtime4);
    
    disp('Table 4, Column 2')
    [runtime1,~] = BrusselatorLOD2D_IFETDRDP(1/100,101);
    fprintf('Time %f\n',runtime1);
    
    [runtime2,~] = BrusselatorLOD2D_IFETDRDP(1/200,201);
    fprintf('Time %f\n',runtime2);
    
    [runtime3,~] = BrusselatorLOD2D_IFETDRDP(1/400,401);
    fprintf('Time %f\n',runtime3);
    
    [runtime4,~] = BrusselatorLOD2D_IFETDRDP(1/800,801);
    fprintf('Time %f\n',runtime4);
    
    disp('Figure 5')
    BrusselatorNewmann3Da_IFETDRDP(0.001,11);
    
    disp('Figure 6 Left')
    BrusselatorNewmann3Db_IFETDRDPcheck(0.001,11);
    
    disp('Figure 6 Right')
    BrusselatorNewmann3Db_IFETDRDPcheck2(0.001,11);
    
    disp('Ginzburg 2D Figure 7 Left/Right and Timing')
    runtime = periodic2D(1/20,400);
    fprintf('Time %f\n',runtime);
    tic;ginzburg2D();runtime=toc;
    fprintf('ETDRK4 Time: %f\n',runtime)
    
    disp('Ginzburg 3D Timing')
    runtime = periodic3D(1/20,50);
    fprintf('Time %f\n',runtime);
    tic;ginzburg3D();runtime=toc;
    fprintf('ETDRK4 Time: %f\n',runtime)
    
    disp('Ginzburg 2D Figure 8 Left/Right and Timing')
    runtime = periodic2Dsmooth(1/20,200);
    fprintf('Time %f\n',runtime);
    disp('Table 5')
    runtime = periodic2Dsmooth(1/20,50);
    fprintf('Time %f\n',runtime);
    runtime = periodic2Dsmooth(1/20,100);
    fprintf('Time %f\n',runtime);
    runtime = periodic2Dsmooth(1/20,200);
    fprintf('Time %f\n',runtime);
    runtime = periodic2Dsmooth(1/20,400);
    fprintf('Time %f\n',runtime);

    disp('Ginzburg 3D Figure 9 Left/Right and Timing')
    runtime = ginz3d(1/20,200);
    fprintf('Time %f\n',runtime);
    disp('Table 6')
    runtime = ginz3d(1/20,50);
    fprintf('Time %f\n',runtime);
    runtime = ginz3d(1/20,100);
    fprintf('Time %f\n',runtime);
    runtime = ginz3d(1/20,200);
    fprintf('Time %f\n',runtime);

    disp('Schroedinger')
    [runtime,~]=Shrodinger1D_ETDRDP(1/100,1/2);
    fprintf('Time %f\n',runtime);
    [runtime,~]=Shrodinger_2DN_IFETDRDP(0.0125,79);
    fprintf('Time %f\n',runtime);
end

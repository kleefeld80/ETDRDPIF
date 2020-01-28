% Run in Matlab
convtest_Brusselator() % timings are printed
X=load('BrusselatorDATA.mat');
efficiency_plot_Brusselator2D(X.Time_mat,X.Error_mat)
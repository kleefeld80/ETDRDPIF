% Run in Matlab
conv_test_ekinetics() % timings are printed
X=load('EkineticsDATA0202.mat');
efficiency_plot_ekinetics(X.Time_mat,X.Error_mat)
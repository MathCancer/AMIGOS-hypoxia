M = read_microenvironment('output00000007_microenvironment0.mat');
Titles(1)={'Oxygen'};
Titles(2)={'ECM'};
Titles(3)={'VEGF'};
close all
figure()
plot_microenvironment(M,Titles)

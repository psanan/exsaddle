% An example of how to load and visualize data with MATLAB/Octave

% Update this according to your PETSc installation
addpath ~/petsc-maint/share/petsc/matlab

close all
clear

%% Read PETSc binary files (slow!)
sol=PetscBinaryRead('solution.petscbin');
%A = PetscBinaryRead('operator_0.petscbin');
BA=PetscBinaryRead('preconditioned_operator.petscbin');


%% Compute spectrum (slow!)
[V,D]=eig(full(BA));
e=diag(D);
max(max(abs(BA-BA')))
es = (sort(real(e))); 
max(abs(imag(e)))
epos = es.*(es>0); 
eneg=abs(es.*(es<0)); 

%% Plot Solution
figure('Position',[0 0 800 800])
mx = 20;
plotfield(sol,mx); % depends on knowing mxa nd mx=my! Will error unless mx/my are correct
shg

%% Plot Spectrum
figure('Position',[100 100 800 400])
subplot(1,2,1)
semilogy(epos,'r+');
hold on; 
semilogy(eneg,'bo'); 
title('All eigenvalues')
shg

subplot(1,2,2)
semilogy(epos(1:40),'r+'); hold on
semilogy(eneg(1:40),'bo'); 
title('A few eigenvalues')


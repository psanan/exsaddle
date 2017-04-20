function [] = plotfield(v,mx,my)
% Plot a Q2-Q1 2D solution 
% obviously this is very brittle and depends on assumptions about the
% ordering

if (nargin<3) my = mx; end
   
% Sizes
dim = 2;
nx_q2 = (2*mx+1);
ny_q2 = (2*my+1);
n_q2 = nx_q2*ny_q2;
nx_q1 = mx+1;
ny_q1 = my+1;
n_q1 = nx_q1*ny_q1;
n = length(v);
if (n ~= dim * n_q2 + n_q1)
   error('Size error: v incompatible with my and mx'); 
end

% Process data to plot
[xx,yy] = meshgrid(1:nx_q2,1:ny_q2);
v_uy = reshape(v(2:2:2*n_q2),ny_q2,nx_q2)';
v_ux = reshape(v(1:2:2*n_q2),ny_q2,nx_q2)';
v_umag = (v_ux.^2 + v_uy.^2).^(1/2);
[xx_p,yy_p] = meshgrid(1:nx_q1,1:ny_q1);
v_p = reshape(v(2*n_q2+1:end),ny_q1,nx_q1)';

%Plot
nplots_x = 2; nplots_y=2;

subplot(nplots_y,nplots_x,1)
%contour(xx,yy,v_umag); colorbar
surf(xx,yy,v_umag); colorbar
axis equal
title('|u|')

subplot(nplots_y,nplots_x,2)
quiver(xx,yy,v_ux,v_uy)
axis equal
title('u');

subplot(nplots_y,nplots_x,3)
contour(xx_p,yy_p,v_p); colorbar
axis equal
title('p')

subplot(nplots_y,nplots_x,4);
surf(xx_p,yy_p,v_p);colorbar
title('p')



end
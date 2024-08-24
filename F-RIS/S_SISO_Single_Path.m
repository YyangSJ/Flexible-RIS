clear all 
close all
f=10e9;
lambda=3e8/f;

alpha1=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1));
beta=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1));
gamma=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1));
theta_B=unifrnd(-1,1);
theta_U=unifrnd(-1,1);
phi_B=unifrnd(-1,1);
phi_U=unifrnd(-1,1);

alpha1=exp(1i*pi/4);
beta=exp(1i*pi/4);
gamma=exp(1i*pi/4);
theta_B=sqrt(2)/2;
theta_U=-sqrt(2)/2;
phi_B=sqrt(2)/2;
phi_U=-sqrt(2)/2;
 
zz=-lambda*2:lambda/150:lambda*2;
xx=-lambda*2:lambda/150:lambda*2;
vv=0:2*pi/150:2*pi;
h_d=gamma;
for i=1:length(zz)
    for j=1:length(xx)
        for k=1:length(vv)
            h_BR=alpha1*exp(-1i*2*pi/lambda*(theta_B*zz(i)+phi_B*xx(j)));
            h_RU=beta*exp(-1i*2*pi/lambda*(theta_U*zz(i)+phi_U*xx(j)));
            h_cas=h_RU'*h_BR; 
            y(i,j,k)=abs(h_d+h_cas*exp(1i*vv(k)))^2;
        end
    end
end
y=y/max(y(:)); 
n = 1000; % 设置渐变的步数
cmap1 = summer(n)/0.88;
cmap2 = parula(n);
alpha = linspace(0, 1, n)';
alpha_mat = repmat(alpha, 1, size(cmap1, 2)); % 将 alpha 扩展为矩阵

cmap_blend = (1 - alpha_mat) .* cmap1 + alpha_mat .* cmap2;colormap(cmap_blend);


 subplot(2, 3, [1 2 4 5]);
xslice = [zz(1)  zz(length(zz))];                               % define the cross sections to view
yslice = [xx(1) xx(length(zz))];          
  zslice = [0 2*pi];
SL=slice(zz, xx, vv, y, xslice, yslice, zslice)    % display the slices
set(SL, 'EdgeColor', 'none');
colormap(cmap_blend);
cb = colorbar;

 xticks([-2*lambda 0 2*lambda]);
 xlim([-2*lambda 2*lambda]);
xticklabels({'-2\lambda', '0', '2\lambda'});
xlabel('$z$-coordinate', 'Interpreter', 'latex')
 yticks([-2*lambda 0 2*lambda]);
 ylim([-2*lambda 2*lambda]);
yticklabels({'-2\lambda', '0', '2\lambda'});
ylabel('$x$-coordinate', 'Interpreter', 'latex')
zticks([0 pi 2*pi]);
zticklabels({'0', '\pi', '2\pi'});
zlim([0 2*pi]);
zlabel('Reflective phase $v$ [rad]', 'Interpreter', 'latex')
 set(gca, 'FontSize', 22, 'LineWidth', 2); 
subplot(2, 3, 3);

xslice = [];                               % define the cross sections to view
yslice = [];          
  zslice = [0  ];
SL=slice(zz, xx, vv, y, xslice, yslice, zslice)    % display the slices
title('$v=0$', 'Interpreter', 'latex')
set(SL, 'EdgeColor', 'none');
colormap(cmap_blend);
 set(gca, 'FontSize', 13, 'LineWidth', 1.2); 
 xticks([-2*lambda 0 2*lambda]);
xticklabels({'-2\lambda', '0', '2\lambda'});
 yticks([-2*lambda 0 2*lambda]);
yticklabels({'-2\lambda', '0', '2\lambda'});
zticks([0 pi 2*pi]);
zticklabels({'0', '\pi', '2\pi'});
zlim([0 2*pi]);
 set(gca, 'FontSize', 20, 'LineWidth', 1.8); 
subplot(2, 3,6);
xslice = [zz(ceil(length(zz)/2))];                               % define the cross sections to view
yslice = [xx(ceil(length(zz)/2))];          
  zslice = [  ];
SL=slice(zz, xx, vv, y, xslice, yslice, zslice)    % display the slices
title('$x=z=0$','Interpreter', 'latex')
set(SL, 'EdgeColor', 'none');
colormap(cmap_blend);
 set(gca, 'FontSize', 13, 'LineWidth', 1.2); 
 xticks([-2*lambda 0 2*lambda]);
xticklabels({'-2\lambda', '0', '2\lambda'});
 yticks([-2*lambda 0 2*lambda]);
yticklabels({'-2\lambda', '0', '2\lambda'});
zticks([0 pi 2*pi]);
zticklabels({'0', '\pi', '2\pi'});
 set(gca, 'FontSize', 20, 'LineWidth', 1.8); 
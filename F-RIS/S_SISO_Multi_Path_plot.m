clear all
close all

f=10e9;
lambda=3e8/f;
P=2;
L=2;
realization=1;
for reali=1:realization
     
    
    alpha1=[exp(1i*pi/6),exp(1i*pi/3)];
    beta=[exp(1i*pi/4),exp(1i*pi/2)];
    gamma=[exp(1i*pi/4)];
%     theta_B=[sqrt(3)/3,sqrt(3)/2];
%     theta_U=[-sqrt(3)/3,-sqrt(3)/2];
%  phi_B=[sqrt(3)/3,sqrt(3)/2];
%   phi_U=[-sqrt(3)/3,-sqrt(3)/2];
vartheta_B=[pi/4,-pi/3];
vartheta_U=[pi/2,pi/3*2];
varphi_B=[pi/3,pi/4];
varphi_U=[pi/6*5,pi/2];
      theta_B=cos(vartheta_B);
    theta_U=cos(vartheta_U);
 phi_B=cos(varphi_B).*sin(vartheta_B);
  phi_U=cos(varphi_U).*sin(vartheta_U); 
  
    zz=-lambda*2:lambda/150:lambda*2;
    xx=-lambda*2:lambda/150:lambda*2;
    vv=0:2*pi/150:2*pi;
    
    h_d=gamma;
    
    
    for i=1:length(zz)
        for j=1:length(xx)
            for k=1:length(vv)
                h_cas=0;
                for l=1:L
                for p=1:P
                    h_cas=h_cas+alpha1(l)*conj(beta(p))*exp(1i*2*pi/lambda*((theta_B(L)-theta_U(p))*zz(i)+(phi_B(L)-phi_U(p))*xx(j)));
                end
                end
                y(i,j,k)=abs(h_d+h_cas*exp(1i*vv(k)))^2;
            end
        end
    end
%     kk=zeros(P,1);
%     kk_len=-5:1:5;
%     for iter=1:5
%         for k_ind=1:P
%             y_position_temp=[];
%             for k_val=kk_len
%                 kk(k_ind)=k_val*2*pi;
%                 Psi=2*pi/lambda*[theta_cas, phi_cas];
%                 mu=angle(h_d)-angle(alpha*beta)+kk;
%                 bmp=inv(Psi'*Psi)*Psi'*mu;
%                 h_cas=0;
%                 for p=1:P
%                     h_cas=h_cas+alpha*beta(p)*exp(1i*2*pi/lambda*(theta_cas(p)*bmp(1)+phi_cas(p)*bmp(2)));
%                 end
%                 y_position_temp =[y_position_temp abs(h_d+h_cas)];
%             end
%             [k_max,indk_max]=max(y_position_temp);
%             kk(k_ind)=kk_len(indk_max);
%         end
%     end
%     y_position(reali)=k_max;
%     
%     position_only_max(reali)=max(max(y(:,:,34)));
%     phase_only_max(reali)=max(y(133,133,:));
%     position_phas_joint_max(reali)=max(vec(y));
end
% position_only_max_mean=mean(position_only_max)
% phase_only_max_mean=mean(phase_only_max)
% position_phas_joint_max_mean=mean(position_phas_joint_max)
% y_position_mean=mean(y_position)


 
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
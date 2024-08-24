function pw= MM_SNR(x1,x2,x3,x4,z1,z2,z3,z4,v1,v2,v3,v4,L,P,alpha,beta,gamma,lambda,theta_B,phi_B,theta_R,phi_R)
 
x=[x1,x2,x3,x4].';
z=[z1,z2,z3,z4].';
v=[v1,v2,v3,v4].';

% h_BR=zeros(4,1);h_RU=zeros(4,1);
% for l=1:L
%     h_BR=h_BR+alpha(l)*exp(1i*2*pi/lambda*(theta_B(l)*z+phi_B(l)*x));
% end
% for p=1:P
%     h_RU=h_RU+beta(p)*exp(1i*2*pi/lambda*(theta_R(p)*z+phi_R(p)*x));
% end
% pw= -abs(h_RU'*diag(h_BR)*exp(1i*v)+gamma)^2;

pw_cas=0;
for l=1:L
    for p=1:P
        for n=1:4
            pw_cas=pw_cas+alpha(l)*conj(beta(p))*exp(1i*v(n))*exp(1i*2*pi/lambda*((theta_B(l)-theta_R(p))*z(n)+(phi_B(l)-phi_R(p))*x(n)));
        end
    end
end
pw=-abs(pw_cas+gamma)^2;
end


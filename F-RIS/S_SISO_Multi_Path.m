clear all
close all

f=10e9;
lambda=3e8/f;
L=1;

realization=122222;
len1=50;
len2=100;
R=lambda*2;
for reali=1:realization
    reali
    

    for pp=1:1
        P=pp;
        
        gamma=1/sqrt(2)*(normrnd(0,3)+1i*normrnd(0,3));
        
        alpha1=1/sqrt(2)*(normrnd(0,2,L,1)+1i*normrnd(0,2,L,1));
        beta=1/sqrt(2)*(normrnd(0,4,P,1)+1i*normrnd(0,4,P,1));
        
        
        theta_B=unifrnd(-1,1,L,1);
        theta_U=unifrnd(-1,1,P,1);
        phi_B=unifrnd(-1,1,L,1);
        phi_U=unifrnd(-1,1,P,1);
        
        zz=-R:lambda/len1:R;
        xx=-R:lambda/len1:R;
        vv=0:2*pi/len2:2*pi;
        
        h_d=gamma;
        %         y=zeros(length(zz),length(xx),length(vv));
        %         for i=1:length(zz)
        %             for j=1:length(xx)
        %                 for k=1:length(vv)
        %                     h_cas=0;
        %                     for l=1:L
        %                         for p=1:P
        %                             h_cas=h_cas+sqrt(1/L/P)*alpha1(l)*conj(beta(p))*exp(1i*2*pi/lambda*((theta_B(l)-theta_U(p))*zz(i)+(phi_B(l)-phi_U(p))*xx(j)));
        %                         end
        %                     end
        %                     y(i,j,k)=abs(h_d+h_cas*exp(1i*vv(k)))^2;
        %                 end
        %             end
        %end
%                 position_only_max(reali,p)=max(max(y(:,:,1)));
%         phase_only_max(reali,p)=max(y(ceil(length(zz)/2),ceil(length(xx)/2),:));
%         position_phas_joint_max(reali,p)=max(vec(y));
        
        v=angle(gamma)-angle(sum(kron(alpha1,conj(beta))));
        h_cas=0;
        for l=1:L
            for p=1:P
                h_cas=h_cas+sqrt(1/L/P)*alpha1(l)*conj(beta(p))*exp(1i*2*pi/lambda*((theta_B(l)-theta_U(p))*0+(phi_B(l)-phi_U(p))*0));
            end
        end
        y(reali,pp)=abs((h_d)+(h_cas*exp(1i*v)))^2;
        
        
    end
end
position_only_max_mean=mean(position_only_max,1)
phase_only_max_mean=mean(phase_only_max,1)
position_phas_joint_max_mean=mean(position_phas_joint_max,1)
warning off
clear all;
close all;
 

N=16; % 4*4 RIS on a x-z plane
T=9; % number of element movement actions
Nx=12; Nz=12; % 12*12 virtual RIS formed by T element movement actions
NT=Nx*Nz;
f=10e9;
lambda=3e8/f;
d=lambda/2;
x=(-(Nx-1)/2:(Nx-1)/2).'*d;
z=(-(Nz-1)/2:(Nz-1)/2).'*d;

for nx=1:Nx
    for nz=1:Nz
        x_pos(nz,nx)=x(nx);
        z_pos(nz,nx)=z(nz);
    end
end
Gx=12;
Gz=12;
G=Gx*Gz;
 
A=kron(dftmtx(Gx),dftmtx(Gz));
A_temp=A;

cols_to_remove = all(A == 1); % remove the atom with theta=phi=0
A(:, cols_to_remove) = [];
P=4;
L=4;

pow=10^(20/10);
realization=1000;




for reali=1:realization
    reali
    QT_set=2:2:N;

    alpha=1/sqrt(2)*(normrnd(0,1,L,1)+1i*normrnd(0,1,L,1));
    beta=1/sqrt(2)*(normrnd(0,1,P,1)+1i*normrnd(0,1,P,1));
    gamma=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1));
 
    for qt=1:length(QT_set)
        QT=QT_set(qt);

        W=zeros(N*T,QT*T);
        for t=1:T
            W0=normrnd(0,1,N,QT)+1i*normrnd(0,1,N,QT);
            W0=W0./abs(W0);
            W((t-1)*N+1:t*N,(t-1)*QT+1:t*QT)=W0;
        end
        h_cas=zeros(NT,1);
        i=0; 
        for p=1:P
            for l=1:L
                i=i+1; 
                h_cas=h_cas+1/sqrt(L*P)*alpha(l)*conj(beta(p))*A(:,rand_angle(i));
            end
        end
        n=sqrt(1/pow)*1/sqrt(2)*(normrnd(0,1,QT*T,1)+normrnd(0,1,QT*T,1)*1i);

        y= W.'*h_cas+gamma*ones(QT*T,1)+n;
        II=400;
        PHI=[ones(QT*T,1), W.'*A];
        PHI_co=PHI'*PHI; 
        for n1=1:(1+size(A,2))
            for n2=1:(1+size(A,2))
                corr_matrix(n1,n2)=abs(PHI(:,n1)'*PHI(:,n2));
            end
        end

        K2=12; 
        [V, ~] = eigs((diag(sum(corr_matrix, 2))- corr_matrix+1e-6*eye((1+size(A,2)))), K2, 'sm'); 
        V_normalized = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))); 
        idx = kmeans((V_normalized), K2); 
        [cluster_idx,sort_idx]=sort(idx);
        for k=1:K2
            D(k)=sum(cluster_idx == k);
        end


        tic
        [ar1,hat_sup]=cs_somp(y,PHI,1+size(A,2),L*P);
        toc
        %mu_vsbl=V_SBL(y,PHI,II);
        tic
        mu_mfvsbl=MFV_SBL(y,PHI,PHI_co,II);
        toc 
        tic
        mu_fmfsbl=FMFSBL(PHI,y,II);
        toc 
        tic
        mu_fcmfsbl=CMFV_SBL(y,PHI,PHI_co,K2,II,D,sort_idx,idx);
        toc
        %mu_admm=ADMM(PHI,y,0.7,0.15,II);
        tic
        mu_FISTA=FISTA(PHI,y,zeros(G,1),1.05,1.01,0.15,II,1, 0);
        toc
        AT=[ones(Nx*Nz,1),A]; 

        nmse_cas_omp(reali,qt)=norm(h_cas-A*hat_sup(2:end))^2/norm(h_cas)^2;
        nmse_cas_mfvsbl(reali,qt)=norm(h_cas-A*mu_mfvsbl(2:end))^2/norm(h_cas)^2; 
        nmse_cas_fmfvsbl(reali,qt)=norm(h_cas-A*mu_fmfsbl(2:end))^2/norm(h_cas)^2; 
        nmse_cas_fcmfvsbl(reali,qt)=norm(h_cas-A*mu_fcmfsbl(2:end))^2/norm(h_cas)^2; 
        nmse_cas_fista(reali,qt)=norm(h_cas-A*mu_FISTA(2:end,II))^2/norm(h_cas)^2;
 
        nmse_direct_omp(reali,qt)=norm(gamma-hat_sup(1))^2/norm(gamma)^2; 
        nmse_direct_mfvsbl(reali,qt)=norm(gamma-mu_mfvsbl(1))^2/norm(gamma)^2; 
        nmse_direct_fmfvsbl(reali,qt)=norm(gamma-mu_fmfsbl(1))^2/norm(gamma)^2; 
        nmse_direct_fcmfvsbl(reali,qt)=norm(gamma-mu_fcmfsbl(1))^2/norm(gamma)^2;
        nmse_direct_fista(reali,qt)=norm(gamma-mu_FISTA(1,II))^2/norm(gamma)^2;

        %     [t2_index, t1_index] = ind2sub([Gx, Gz],ar1);
        %     phi_est=phi_dic(t1_index)
        %     phi_true=vec(phi_R)
        %     theta_est=phi_dic(t2_index)
        %     theta_true=vec(theta_R)
        %
        %     vars = [
        %         optimizableVariable('phi1',[phi_est(1)-2/Gx,phi_est(1)+2/Gx],'Type','real'),
        %         optimizableVariable('phi2',[phi_est(2)-2/Gx,phi_est(2)+2/Gx],'Type','real'),
        %         optimizableVariable('phi3',[phi_est(3)-2/Gx,phi_est(3)+2/Gx],'Type','real'),
        %         optimizableVariable('phi4',[phi_est(4)-2/Gx,phi_est(4)+2/Gx],'Type','real') ,
        %         optimizableVariable('theta1',[theta_est(1)-2/Gz,theta_est(1)+2/Gz],'Type','real'),
        %         optimizableVariable('theta2',[theta_est(2)-2/Gz,theta_est(2)+2/Gz],'Type','real'),
        %         optimizableVariable('theta3',[theta_est(3)-2/Gz,theta_est(3)+2/Gz],'Type','real'),
        %         optimizableVariable('theta4',[theta_est(4)-2/Gz,theta_est(4)+2/Gz],'Type','real')
        %         ];
        %     fun = @(x)BOMLE(x.phi1,x.phi2,x.phi3,x.phi4,x.theta1,x.theta2,x.theta3,x.theta4,z_pos,x_pos,y,lambda); %Create Function Handles
        %
        %     ang_ini = [phi_est(1),phi_est(2),phi_est(3),phi_est(4),theta_est(1),theta_est(2),theta_est(3),theta_est(4)];
        %     T_ang_ini = array2table(ang_ini);
        %     fun = @(x)BOMLE(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),z_pos,x_pos,y,lambda); %Create Function Handles
        %
        %     tic
        %     ang_est=fminsearch(fun,ang_ini)
        %     %ang_est=fminunc(fun,ang_ini)
        %
        %     toc
        %     %      results = bayesopt(fun,vars,'InitialX',T_ang_ini,'MaxObjectiveEvaluations',50,'PlotFcn',[])
        %     %
        %     %     ang_est= table2array(results.XAtMinObjective);
        %     phi_est2=[ang_est(1),ang_est(2),ang_est(3),ang_est(4)];
        %     theta_est2=[ang_est(5),ang_est(6),ang_est(7),ang_est(8)];
        %     A_est=[];
        %     for i=1:length(phi_est)
        %         A_est=[A_est exp(1i*2*pi/lambda*(theta_est2(i)*z_pos(:)+phi_est2(i)*x_pos(:)))];
        %     end
        %
        %     nmse5(reali)=norm(h-A_est*inv(A_est'*A_est)*A_est'*y)^2/norm(h)^2
    end
end

NMSE_OMP=mean(nmse_omp)
NMSE_MFVSBL=mean(nmse_mfvsbl,1)
NMSE_FMFVSBL=mean(nmse_fmfvsbl,1)
NMSE_FCMFVSBL=mean(nmse_fcmfvsbl,1)
NMSE_FISTA=mean(nmse_fista,1)

NMSE_CAS_OMP=mean(nmse_cas_omp,1)
%NMSE_CAS_VSBL=mean(nmse_cas_vsbl,1)
NMSE_CAS_MFVSBL=mean(nmse_cas_mfvsbl,1)
%NMSE_CAS_GMFVSBL=mean(nmse_cas_gmfvsbl,1)
NMSE_CAS_FMFVSBL=mean(nmse_cas_fmfvsbl,1)
%NMSE_CAS_CMFVSBL=mean(nmse_cas_cmfvsbl,1)
NMSE_CAS_FCMFVSBL=mean(nmse_cas_fcmfvsbl,1)
NMSE_CAS_FISTA=mean(nmse_cas_fista,1)

NMSE_DIRECT_OMP=mean(nmse_direct_omp,1)
%NMSE_DIRECT_VSBL=mean(nmse_direct_vsbl,1)
NMSE_DIRECT_MFVSBL=mean(nmse_direct_mfvsbl,1)
%NMSE_DIRECT_GMFVSBL=mean(nmse_direct_gmfvsbl,1)
NMSE_DIRECT_FMFVSBL=mean(nmse_direct_fmfvsbl,1)
%NMSE_DIRECT_CMFVSBL=mean(nmse_direct_cmfvsbl,1)
NMSE_DIRECT_FCMFVSBL=mean(nmse_direct_fcmfvsbl,1)
NMSE_DIRECT_FISTA=mean(nmse_direct_fista,1)
figure
semilogy(NMSE_CAS_FISTA)
hold on
semilogy(NMSE_CAS_OMP)
hold on
semilogy(NMSE_CAS_MFVSBL)
hold on
semilogy(NMSE_CAS_FCMFVSBL)
hold on
semilogy(NMSE_CAS_FMFVSBL)
figure
semilogy(NMSE_DIRECT_FISTA)
hold on
semilogy(NMSE_DIRECT_OMP)
hold on
semilogy(NMSE_DIRECT_MFVSBL)
hold on
semilogy(NMSE_DIRECT_FCMFVSBL)
hold on
semilogy(NMSE_DIRECT_FMFVSBL)
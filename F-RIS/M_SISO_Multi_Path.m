clear all
close all
 
f=10e9;
lambda=3e8/f;
N=4;
realization=3;
L=3;P=3;
r=lambda*2;
dm=lambda/2;
for reali=1:realization
    alpha=1/sqrt(2)*(normrnd(0,1,L,1)+1i*normrnd(0,1,L,1));
    beta=1/sqrt(2)*(normrnd(0,1,P,1)+1i*normrnd(0,1,P,1));
    gamma=1/sqrt(2)*(normrnd(0,1)+1i*normrnd(0,1));
    theta_B=unifrnd(-1,1,L,1);
    phi_B=unifrnd(-1,1,L,1);
    theta_R=unifrnd(-1,1,P,1);
    phi_R=unifrnd(-1,1,P,1);
    
    %% BO for joint {x,z,v} optimization
%     fun = @(x)MM_SNR(x.x1, x.x2, x.x3, x.x4, x.z1, x.z2, x.z3, x.z4, x.v1, x.v2, x.v3, x.v4,L,P,alpha,beta,gamma,lambda,theta_B,phi_B,theta_R,phi_R)
%     % 定义变量
%     vars = [
%         optimizableVariable('x1', [-r, r]),
%         optimizableVariable('x2', [-r, r]),
%         optimizableVariable('x3', [-r, r]),
%         optimizableVariable('x4', [-r, r]),
%         optimizableVariable('z1', [-r, r]),
%         optimizableVariable('z2', [-r, r]),
%         optimizableVariable('z3', [-r, r]),
%         optimizableVariable('z4', [-r, r]),
%         optimizableVariable('v1', [0, 2*pi]),
%         optimizableVariable('v2', [0, 2*pi]),
%         optimizableVariable('v3', [0, 2*pi]),
%         optimizableVariable('v4', [0, 2*pi])
%         ];
%     
%     con=@(x)constraint(x.x1, x.x2, x.x3, x.x4, x.z1, x.z2, x.z3, x.z4,dm)
%     results = bayesopt(fun,vars,'XConstraintFcn',con,'MaxObjectiveEvaluations',60,'ParallelMethod','min-observed')
    %results = bayesopt(fun, vars);
    
    %% BO for {x,z} optimization
    fun = @(x)MM_SNR(x.x1, x.x2, x.x3, x.x4, x.z1, x.z2, x.z3, x.z4, 0, 0,0,0,L,P,alpha,beta,gamma,lambda,theta_B,phi_B,theta_R,phi_R)
    % 定义变量
    vars = [
        optimizableVariable('x1', [-r, r]),
        optimizableVariable('x2', [-r, r]),
        optimizableVariable('x3', [-r, r]),
        optimizableVariable('x4', [-r, r]),
        optimizableVariable('z1', [-r, r]),
        optimizableVariable('z2', [-r, r]),
        optimizableVariable('z3', [-r, r]),
        optimizableVariable('z4', [-r, r])
        ];
    con=@(x)constraint(x.x1, x.x2, x.x3, x.x4, x.z1, x.z2, x.z3, x.z4,dm)
    results = bayesopt(fun,vars,'XConstraintFcn',con,'MaxObjectiveEvaluations',true,'ParallelMethod','min-observed','PlotFcn',[]);
    %% BO for v optimization
    xx=[-(2-1)/2:(2-1)/2]*lambda/2;
    zz=[-(2-1)/2:(2-1)/2]*lambda/2;
    xxn=[xx(1), xx(1), xx(2), xx(2)];
    zzn=[zz(1), zz(2), zz(1), zz(2)];
    fun = @(x)MM_SNR(xxn(1), xxn(2), xxn(3), xxn(4), zzn(1), zzn(2), zzn(3), zzn(4), x.v1, x.v2,x.v3,x.v4,L,P,alpha,beta,gamma,lambda,theta_B,phi_B,theta_R,phi_R)
    % 定义变量
    vars = [
        optimizableVariable('v1', [0, 2*pi]),
        optimizableVariable('v2', [0, 2*pi]),
        optimizableVariable('v3', [0, 2*pi]),
        optimizableVariable('v4', [0, 2*pi])
        ];
    results = bayesopt(fun,vars,'MaxObjectiveEvaluations',60,'ParallelMethod','min-observed')
    %% Closed form for v optimization
    for n=1:N
        asum=0;
        for l=1:L
            for p=1:P
                asum=asum+alpha(l)*conj(beta(p))*exp(1i*2*pi/lambda*((theta_B(l)-theta_R(p))*zzn(n)+(phi_B(l)-phi_R(p))*xxn(n)))
            end
        end
        v_cf(n)=angle(gamma)-angle(asum);
    end
    MM_SNR(xxn(1), xxn(2), xxn(3), xxn(4), zzn(1), zzn(2), zzn(3), zzn(4),v_cf(1),v_cf(2),v_cf(3),v_cf(4),L,P,alpha,beta,gamma,lambda,theta_B,phi_B,theta_R,phi_R)
end
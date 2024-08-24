function [mu] =  CMFV_SBL(y,PHI,PHI_co,K,II,D,sort_idx,idx)

Esig2=1;
[M,G]=size(PHI);
Lambda=eye(G);


a = 1e-6;
b = 1e-6;
c = 1e-6;
d = 1e-6;


Py=PHI'*y;
mu=zeros(G,1);
%mu=PHI'*inv(PHI*PHI')*y;
% D=G/K;
%idx
%  sort_idx=1:G;
%  D=ones(K,1)*G/K; %  number of each cluster
%  D(1)=12;D(2)=4;
SIG=zeros(G,G);
G_set=1:G;
for iter=1:II
    mu_old = mu;
    for k=1:K
        %         PHI_t=PHI;
        %        PHI_t(:,i)=[];
        if D(k)==0
            continue; % pass cluster k
        end
        mu_t=mu;

        cluter_k_idx=sort_idx(sum(D(1:k-1))+1:sum(D(1:k)));

        mu_t(cluter_k_idx)=[];
        PHI_t=PHI_co(cluter_k_idx,:);
        PHI_t(:,cluter_k_idx)=[];
        SIG(cluter_k_idx,cluter_k_idx)=inv(Esig2*PHI_co(cluter_k_idx,cluter_k_idx)+Lambda(cluter_k_idx,cluter_k_idx));
        mu(cluter_k_idx)=Esig2*SIG(cluter_k_idx,cluter_k_idx)*(Py(cluter_k_idx)-PHI_t*mu_t);
        % mu(i)=Esig2*SIG(i,i)*PHI(:,i)'*(y-PHI_t*mu_t)

    end
    for i=1:length(G_set)
        rho(G_set(i))=(1/2+a)/((abs(mu(G_set(i)))^2+SIG(G_set(i),G_set(i)))/2+b);
    end
    % Esig2=(c+G/2)/((norm(y)^2-2*real(Py'*mu)+trace(PHI_co*(mu*mu'+SIG)))/2+d);

    idx_all_K=sort_idx(1:sum(D(1:K)));
    Esig2=(c+M/2)/((norm(y)^2-2*real(Py(idx_all_K)'*mu(idx_all_K))+mu(idx_all_K)'*PHI_co(idx_all_K,idx_all_K)*mu(idx_all_K)+sum(diag(PHI_co(idx_all_K,idx_all_K)).*diag(SIG(idx_all_K,idx_all_K))))/2+d);

    Lambda=diag(rho);
    ind_remove = find(rho > 5e4);

    G_set=setdiff(G_set,ind_remove);
    rho(ind_remove)=0;
    mu(ind_remove)=0;
    for iq=1:length(ind_remove)
        sort_idx(find(sort_idx==ind_remove(iq)))=[];
        if D(idx(ind_remove(iq)))>0
            D(idx(ind_remove(iq)))=D(idx(ind_remove(iq)))-1;
        end
    end
    if norm(mu-mu_old)/norm(mu)<1e-8
        break
    end
end
end


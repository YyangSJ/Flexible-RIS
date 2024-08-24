function Result=IFSBL(A,y,Max_iter)
    [M,N]=size(A);
    alpha=ones(N,1);
  mu=A'*y;
   %  mu=A'*inv(A*A')*y;
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    S = svd(A,'econ');
    L = 2*S(1)^2;%+0*0.00001; % Lipschitz constant
    delta = 1;
    theta = mu;
    ATY = A'*y;
    ATA = A'*A;
    Ath = A*theta;
    for iter = 1 : Max_iter
        mu_old = mu;
        sigma = 1./(delta*L/2+alpha);% Calculate the disgnoal entries of the covariance matrix
        mu = (delta*theta*L+2*delta*(ATY-ATA*theta)).*sigma/2;% Update mu
        alpha = (0.5+a)./(0.5*(abs(mu).^2+sigma)+b);% Update alpha
%         Ath = A*theta;
        Amu = A*mu;
        %delta=(c+0.5*M)/(d+0.5*(norm(y-Amu)^2+2*real((mu-theta)'*A'*(A*theta-y))+ 0.5*L*norm(mu-theta)^2+0.5*L*sum(sigma)));
         delta = (c+0.5*M)/(d+0.5*((y-2*Amu+Ath)'*(y-Ath)+0.5*L*sum(abs((mu-theta)).^2)+0.5*L*sum(sigma)));% Update noise precision
        theta = mu;% Update theta
        Ath = Amu;
%         if norm(mu-mu_old)/norm(mu)<1e-8
%             break
%         end
    end
    Result.x = mu;
end
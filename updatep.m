function [p q lam] = updatep(t, p, q, K0, alpha, beta, mu)

N = size(t, 2);
lam = zeros(N,1);

for i=1:N
    for j=1:(i-1)
        %% Calculate p( i, j), i > j stands the probability of ONE SIGLE event at day i triggered by "ALL" events at day j
        p(i,j) = K0(j)*wblpdf((i-j), alpha, beta)*t(j);
        
        %% Calculate q( i, j), i > j stands the probability of ONE SIGLE event at day i triggered by "ONE SIGLE" event at day j
        % Note that (t(j)>0) is to see whether there is an event in day j or not
        % If so, we calculate single event at day j. If not, it remains zero.
        q(i,j) = K0(j)*wblpdf((i-j), alpha, beta)*(t(j)>0);
        
    end

    %% probablity ii is background event proportional to mu background rate
    p(i,i)=mu;
    q(i,i)=mu;
    
    %% save intensity at each event for analysis
    lam(i)=sum(p(i,1:i));
    
    %% normalize probabilities
    % Note that p will sum to 1 since all events at day j are considered, q won't.
    p(i,1:i) = p(i,1:i)/lam(i);
    q(i,1:i) = q(i,1:i)/lam(i);
    if sum(sum(isnan(p)))>0
        disp("temp")
    end
end

end
  



%%%%%%%%%%%% Legacy %%%%%%%%%%%%%%%%%%

            % probability i triggered by j is proportional to triggering
            % kernel evaluated at inter-point times and distances
%             p(i,j)=K0(ceil(9*j/T))*wblpdf((i-j),alpha,beta)*t(j);
         %   p(i,j)=K0(ceil(bins*j/T))*wbl_new(i,j,alpha,beta)*t(j);
%            d=i-j;
%             p(i,j)=K0(ceil(bins*j/T))*(exp(-((d-1)/alpha)^beta)-exp(-(d/alpha)^beta))*t(j);
%            p(i,j)=K0(ceil(bins*j/T))*alpha/beta*(gammainc(((d+1)/alpha)^beta,1/beta,'upper')+gammainc(((d-1)/alpha)^beta,1/beta,'upper')-2*gammainc((d/alpha)^beta,1/beta,'upper'))*t(j);
%lam = reshape(lam, 1, N);      
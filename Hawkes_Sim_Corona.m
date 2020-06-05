function [times ]=Hawkes_Sim_Corona( mu, alpha, beta, T, fK0, Tlow, itr, tr_in)
%this program simulates event times, called "times", according
%to a self-exciting point process with exponential triggering kernel and
%parameters mu (const. background rate)
%k0 (branching ratio) and w (exp parameter) on the time interval [0,T]

times = zeros(5000,1);

%first simulate "background" events
%this is done by picking p points where p is Poisson with parameter mu*T
%and then distributing the points uniformly in the interval [0,T]
n_K0cal = T-Tlow;
base = zeros(1, n_K0cal);
n_exh = zeros(1, n_K0cal);
%
ct = 1;
T_offset = (Tlow:T-1);
%
for t_s = 1 : n_K0cal
    t_stamps = (t_s+Tlow) - (1:Tlow);
    intense = wblpdf( t_stamps, alpha, beta).*fK0(1:Tlow).*tr_in;
    intense = sum(intense) + mu;
    base(t_s) = intense;
    %
    n_exh(t_s) = pois( intense );
    %
    times(ct:ct+n_exh(t_s)-1) = rand( n_exh(t_s), 1) + T_offset(t_s);
    ct = ct+n_exh(t_s);
    %
end

n_tol_exh = sum(n_exh);
times(1:n_tol_exh) = sort(times(1:n_tol_exh));
% 
% p=pois(mu*T);
% times(1:p,1)=rand(p,1)*T;
%counts = 1;
%countf = p;
counts = 1;
countf = n_tol_exh;


% Modify the code to drop the simulation that is before T-delta
% Start with the nearest one
%times(1:p,1) = sort(times(1:p,1));

%Next loop through every event and simulate the "offspring"
%even the offspring events can generate their own offspring

%dbstop if error

while((countf-counts)>-1)
    
    p = pois( fK0( ceil( times(counts) ) ) );
    
    %each event generates p offspring according to a Poisson r.v. with parameter k0
    % this generates an exponential r.v. on [t_counts,infty]
    temp = times(counts) + wblrnd( alpha, beta, p, 1); 
    temp(find(temp > T)) = [];
    
    times(1:countf + numel(temp)) = [times(1:countf);temp];
    
    countf = countf + numel(temp);

    counts=counts+1;

    % Check and remove the finished one 
    if( times(counts-1) < Tlow)
        times(counts-1) = [];
        counts = counts-1;
        countf = countf-1;
    end
    
end
data=[times(1:countf)];

times=sort(data);





end


function p=pois(S)

if(S<=100)
    temp=-S;
    L=exp(temp);
    k=0;
    p=1;
    while(p > L)
        k=k+1;
        p=p*rand();
    end
    p=k-1;
else
    p=floor(S+S^.5*randn());
end
end

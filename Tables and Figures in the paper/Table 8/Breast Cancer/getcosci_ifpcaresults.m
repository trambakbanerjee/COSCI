
load('selected.mat');
load('breast_y.mat');

trueclass = breasty;
CER = zeros(30,1);
K = 2;

parfor rep=1:30
    
    G = selected*(selected');
    [V, ~] = eigs(G, K - 1); 
    label = kmeans(V, K, 'replicates', 30);
    
    predclass = label;
    
    % Calculate CER
    
    n = size(selected,1);
    num = 0;
    den = n*(n-1)/2;
    for i= 1:(n-1)
        
        tt = predclass(i)-predclass(i+1:n);
        tt(tt~=0)=1;
        tt = 1-tt;
        dd = trueclass(i)-trueclass(i+1:n);
        dd(dd~=0)=1;
        dd = 1-dd;
        num = sum(abs(tt-dd))+num;
        
    end
    CER(rep,1) = num/den;
end
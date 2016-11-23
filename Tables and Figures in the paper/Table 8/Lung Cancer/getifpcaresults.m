
load('lunggordon_x.mat');

trueclass = [ones(150,1);2*ones(31,1)];
CER = zeros(30,2);

parfor rep=1:30
    
    [label, stats, numselect,signals] = ifpca(lunggordonx, 2);
    
    predclass = label;
    
    % Calculate CER
    
    n = size(lunggordonx,2);
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
    CER(rep,:) = [num/den numselect];
end
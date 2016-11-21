
wd = 'Enter your working directory here';
%%----------------------------------------------------------------------------

reps = 50;
p = 50;
s = [1 2 3 4 5];
noise = 6:p;

%%--------------------------------------------------------------------------------
n = 200;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

for r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,2:(p+1));
   [label, stats, numselect,signals] = ifpca(d',32);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

%%--------------------------------------------------------------------------------
n = 1000;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

for r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,2:(p+1));
   [label, stats, numselect,signals] = ifpca(d',32);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

%%--------------------------------------------------------------------------------
n = 2500;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

for r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,2:(p+1));
   [label, stats, numselect,signals] = ifpca(d',32);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

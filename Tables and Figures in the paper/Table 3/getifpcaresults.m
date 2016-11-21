
wd = 'Enter your working directory here';
%%------------------------------------------------------------------------------
reps = 50;
p = 100;
s = [1 2 3 4 5 6];
noise = 7:p;

%%--------------------------------------------------------------------------------
n = 200;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

parfor r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,1:p);
   [label, stats, numselect,signals] = ifpca(d',96);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

%%--------------------------------------------------------------------------------
n = 1000;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

parfor r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,1:p);
   [label, stats, numselect,signals] = ifpca(d',96);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

%%--------------------------------------------------------------------------------
n = 2500;
loc = strcat(wd,'\',num2str(n),'\');
out = zeros(reps,2);

parfor r = 1:reps
    
   d = readdata(loc, r);
   d = d(:,1:p);
   [label, stats, numselect,signals] = ifpca(d',96);
   out(r,:) = [length(s)-sum(ismember(s,signals)) sum(ismember(noise,signals))]; 
   
end

[sum(out(:,1))/reps std(out(:,1))/sqrt(reps) sum(out(:,2))/reps std(out(:,2))/sqrt(reps)]

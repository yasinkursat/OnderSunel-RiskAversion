clear all
clc
indicator_beowulf =1;
if indicator_beowulf>0
   load graphs\data_sim.txt
   load graphs\def_per.txt
   load graphs\param.txt
   load graphs\delta.txt
   
%    data_sim = graphs_data_sim;
%    def_per = graphs_def_per;
%    param = graphs_param;
%    delta = graphs_delta;
else
    load data_sim.txt
    load def_per.txt
    load param.txt
    load delta.txt
end
r=0.01;
coupon = (r+delta)/(1 + r);

per_num = param(1);  %HOW MANY PERIODS IN EACH SAMPLE
n = param(2);        %HOW MANY SAMPLES

y = zeros(per_num, n);
b = zeros(per_num, n);
q = zeros(per_num, n);
c = zeros(per_num, n);      
tb = zeros(per_num, n);
d = zeros(per_num, n);
e = zeros(per_num, n);
q_rn = zeros(per_num, n);
b_next = zeros(per_num-1, n);
duration = zeros(per_num,n);
recovery = zeros(per_num, n);
b_adj = zeros(per_num, n);

for i=1:n
   y(:,i) = data_sim((i-1)*per_num+1:i*per_num,1); 
   b(:,i) = data_sim((i-1)*per_num+1:i*per_num,2); 
   q(:,i) = data_sim((i-1)*per_num+1:i*per_num,3); 
   c(:,i) = data_sim((i-1)*per_num+1:i*per_num,4); 
   tb(:,i) = data_sim((i-1)*per_num+1:i*per_num,5);
   d(:,i) = data_sim((i-1)*per_num+1:i*per_num,6)-1;
   e(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
   q_rn(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
   yield(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
   b_adj(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
   recovery(:,i) = 0;%data_sim((i-1)*per_num+1:i*per_num,8);
%    duration(:,i) = (((1-delta).*q(:,i) +1)./(q(:,i)))./(((1-delta).*q(:,i) +1)./(q(:,i))-1+delta);
%    duration(:,i) = (((1-delta).*q(:,i)+ coupon)./coupon);
    duration(:,i) = (1+yield(:,i))./(yield(:,i)+delta);
end

b_next = b(2:per_num,:);
spread = (((1+yield)/(1+r)).^4) - 1;
%g = zeros(per_num, n);
%g(1,:) = 1+y(1,:);
%g(2:per_num,:) = 1+(y(2:per_num,:)-y(1:per_num-1,:));

%CREATE OUTPUT SERIES
num = 500; %HOW MANY PERIODS IN EACH SUBSAMPLE
%b1 = zeros(per_num, n);
%for j=1:n
%    y(1,j) = 1;
 %   for i=2:per_num
 %      y(i,j) = y(i-1,j)*g(i,j);
 %      b1(i,j) = y(i-1,j)*b(i,j);
 %   end
%end

last_y = y(per_num - num+1:per_num,:);  %TRIM THE FIRST 9500 OBSERVATIONS
last_c = c(per_num - num+1:per_num,:);  %TRIM THE FIRST 9500 OBSERVATIONS
last_tb = tb(per_num - num+1:per_num,:);  
last_d = d(per_num - num+1:per_num,:);
last_b = b(per_num - num+1:per_num,:);
%last_spreadannual = (((1-delta)./q(per_num - num+1:per_num,:) +1-delta)/1.01).^4 - 1;
last_spreadannual = (((1-delta).*q(per_num - num+1:per_num,:) +coupon)./(1.01*q(per_num - num+1:per_num,:))).^4 - 1;
%last_spread = (((1-delta)./q(per_num - num+1:per_num,:) +1-delta)/1.01) - 1;
last_spread = (((1-delta).*q(per_num - num+1:per_num,:) + coupon)./(1.01*q(per_num - num+1:per_num,:))) - 1;
last_duration = duration(per_num - num+1:per_num,:);



%CREATE VECTORS OF STD OF RETURNS AND CORRELATION BETWEEN RETURNS, OUTPUT
%AND TB. NEED TO TRIM OBSERVATIONS WHILE THE COUNTRY IS IN DEFAULT.

lambda = 1600;
y_trend = hpfilter(last_y,lambda);            %FILTER log(output)
c_trend = hpfilter(last_c, lambda);           %FILTER log(consumption)
tb_trend = hpfilter(last_tb, lambda);         %FILTER tb/output
spread_trend = hpfilter(last_spread, lambda);  %FILTER quaterly spread
spreadannual_trend = hpfilter(last_spreadannual, lambda); %FILTER annualized spread
%COMPUTE DEVIATIONS FROM TREND  
y_dev = last_y - y_trend;
c_dev = last_c - c_trend;
tb_dev = last_tb - tb_trend;
spread_dev = last_spread - spread_trend;
spreadannual_dev = last_spreadannual - spreadannual_trend;

fprintf('Statistics when all sample periods are considered \n')
fprintf('std of output = %f5 \n',100*mean(std(y_dev)))
fprintf('std of cons   = %f5 \n',100*mean(std(c_dev)))
fprintf('std of TB/Y   = %f5 \n',100*mean(std(tb_dev)))

matrix_corr = zeros(n,3);
matrix_corr1 = zeros(n,1);
for i=1:n
    matrix = corrcoef([y_dev(:,i) c_dev(:,i) tb_dev(:,i) spreadannual_dev(:,i)]);
    matrix_corr(i,:) = matrix(1,2:4);

    %COMPUTE CORRELATION BETWEEN FILTERED SPREAD AND TRADE BALANCE
    matrix1 = corrcoef([spreadannual_dev(:,i) tb_dev(:,i)]);
    matrix_corr1(i) = matrix1(1,2);
end
fprintf('corr(c, y)   = %f5 \n',mean(matrix_corr(:,1)))
fprintf('corr(tb,y)   = %f5 \n',mean(matrix_corr(:,2)))
%fprintf('Mean default rate = %f5 \n', 100*mean(sum(last_d)/500))
fprintf('Mean default rate = %f5 \n', 400*mean(sum(d)/per_num))
fprintf('E(R_s)       = %f5 \n',100*mean(mean(last_spreadannual)))
fprintf('Max R_s      = %f5 \n',100*max(max(last_spreadannual)))
fprintf('Duration     = %f5 \n',mean(mean(last_duration))/4)



sample_size = 74;
max_num_def = max(sum(d));
num_observations_vector = zeros(n,1); %NUMBER OF OBSERVATIONS PER SAMPLE (AN OBSERVATION
                               % IS A DEFAULT EPISODE WITH SUFFICIENT
                               % PERIODS BEFORE THE DEFAULT AND NO
                               % EXCLUSION IN BETWEEN
def_per_matrix = zeros(max_num_def,n); %MATRIX OF DEFAULT PERIODS SATISFYING THE ABOVE RESTRICTION

index_acum = 0; %INDEX OF ACUMULATED NUMBER OF DEFAULT EPISODES
for i = 1:n
    def_num = sum(d(:,i));
    if def_num>0
       def_periods = find(d(:,i)>0); %def_per(index_acum+1:index_acum+def_num);  %VECTOR OF DEFAULT PERIODS
       interperiods = zeros(def_num,1);
       interperiods(1) = def_periods(1);             %SEPARATION BETWEEN DEFAULT PERIODS
       interperiods(2:def_num)= diff(def_periods);   %SEPARATION BETWEEN DEFAULT PERIODS
       excl_acum = zeros(def_num,1);   %cumulated number of periods before a default
       excl_acum(1) = 0;
       if def_num>1
           for k = 2:def_num
               excl_acum(k) = sum(e(def_periods(k-1)+1:def_periods(k),i));
           end
       end
       interperiods = interperiods - excl_acum;
       
       index_acum = index_acum + def_num;            %COUNT DEFAULT PERIODS IN i SAMPLE.
                                         %NEEDED TO KEEP TRACK OF SAMPLE
                                         %CHANGES IN def_per
       %FIND SAMPLES CONTAINING AT LEAST sample_size PERIODS BEFORE THE
       %DEFAULT AND THAT THE LAST EXCLUSION PERIOD HAPPENED sample_size + 2
       %PERIODS AGO (NEEDED TO CLEAR THE SAMPLE FROM OUTLIERS)
       indices = find(interperiods>sample_size+24);  
       
       num_observations_vector(i) = length(indices);        %STORE THE NUMBER OF SUCH SAMPLES
       if length(indices)>0
           %ONLY STORE RESULTS IF THERE ARE A POSITIVE NUMBER OF
           %OBSERVATIONS.
           def_per_matrix(1:length(indices),i) = def_periods(indices);  
       end
    end
end


% load data_mat1;



num_observations = sum(num_observations_vector);

std_y_vector = zeros(num_observations,1);
std_c_vector = zeros(num_observations,1);
std_tb_vector = zeros(num_observations,1);
std_ra_vector = zeros(num_observations,1);
std_ra_vector_rn = zeros(num_observations,1);
mean_b_vector = zeros(num_observations,1);
mean_ra_vector = zeros(num_observations,1);
mean_ra_vector_rn = zeros(num_observations,1);
default_b_vector = zeros(num_observations,1);
default_rec_vector = zeros(num_observations,1);

corr_ra_y_vector = zeros(num_observations,1);
corr_y_c_vector = zeros(num_observations,1);
corr_y_bnext_vector = zeros(num_observations,1);
corr_ra_tb_vector = zeros(num_observations,1);
corr_ra_ra_rn_vector = zeros(num_observations,1);

corr_y_tb_vector = zeros(num_observations,1);

b_next_before1 = zeros(5,num_observations);
ra_before1     = zeros(5,num_observations);
b_next_after1  = zeros(4,num_observations);
ra_after1      = zeros(4,num_observations);
new_issuance_before = zeros(5,num_observations);
new_issuance_after  = zeros(4,num_observations);

ra_matrix = zeros(sample_size, num_observations);
y_matrix = zeros(sample_size, num_observations);
c_matrix = zeros(sample_size, num_observations);
y_matrix1 = y_matrix;
ra_matrix1 = ra_matrix;
index_observation = 1;
fr=[];
fr2=[];
fr3=[];
for i=1:n
       for j=1:num_observations_vector(i)
%            ra_vector = (((1-delta)./q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i) +1-delta)/1.01).^4 - 1;
%             ra_vector = (((1-delta).*q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i) +coupon)./...
%                 (1.01*q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i))).^4 - 1;
            ra_vector = (((1-delta).*q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i) +coupon)./...
                (1.017*q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i))).^4 - 1;
            ra_vector_rn = (((1-delta).*q_rn(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i) + coupon)./...
                (1.01*q_rn(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i))).^4 - 1;
                        
            y_vector = y(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            tb_vector = tb(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            c_vector = c(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            b_vector = b(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            b_adj_vector = b_adj(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            by_vector = b_vector./(4*exp(y_vector));
            tau_vector = exp(y_vector) - exp(c_vector);
            b_adj_y_vector = b_adj_vector./(4*exp(y_vector));
            b_next_vector = b_next(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            inflow_vector = q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i).*(b_next(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i)...
                -(1-delta)*b(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i));
            %by_vector = b_vector/(delta+r)./(4*exp(y_vector));
            
            ra_matrix(:, index_observation) = ra_vector;
            y_matrix(:, index_observation) = exp(y_vector);
            c_matrix(:, index_observation) = exp(c_vector);
            by_matrix(:, index_observation) = (by_vector);
            b_adj_y_matrix(:, index_observation) = (b_adj_y_vector);
            b_matrix(:, index_observation) = (b_vector);
            b_next_matrix(:, index_observation) = (b_next_vector);
            b_adj_matrix(:, index_observation) = (b_adj_vector);
            tau_matrix(:, index_observation) =exp(y_vector)-exp(c_vector);
            tb_matrix(:, index_observation) = exp(tb_vector);
            inflow_matrix(:, index_observation) = inflow_vector;
            
            
            lambda=1600;
            hp_y = hpfilter(y_vector,lambda);    %FILTER log 
            hp_bnext = hpfilter(b_next_vector,lambda);    %FILTER log 
            hp_tb = hpfilter(tb_vector,lambda);    %FILTER log 
            hp_ra = hpfilter(ra_vector,lambda);    %FILTER log 
            hp_ra_rn = hpfilter(ra_vector_rn,lambda); 
            hp_c = hpfilter(c_vector,lambda);    %FILTER log 
    
            dev_y = y_vector - hp_y;
            dev_bnext = b_next_vector - hp_bnext;
            dev_tb = tb_vector - hp_tb;
            dev_ra = ra_vector - hp_ra;
            dev_ra_rn = ra_vector_rn - hp_ra_rn;
            dev_c = c_vector - hp_c;
    
            default_b_vector(index_observation) = b(def_per_matrix(j,i),i)/exp(y(def_per_matrix(j,i),i));
            default_rec_vector(index_observation) = recovery(def_per_matrix(j,i),i);
            %by2_vector = b_vector ./ exp(y_vector);
            
            std_ra_vector(index_observation) = std(dev_ra);
            std_ra_vector_rn(index_observation) = std(dev_ra_rn);
            std_c_vector(index_observation) = std(dev_c);
            std_y_vector(index_observation) = std(dev_y);
            std_tb_vector(index_observation) = std(dev_tb);

            matrix = corrcoef(dev_y, dev_c);
            corr_y_c_vector(index_observation) = matrix(1,2);

            matrix = corrcoef(dev_y, dev_tb);
            corr_y_tb_vector(index_observation) = matrix(1,2);
            
            matrix = corrcoef(dev_y, dev_bnext);
            corr_y_bnext_vector(index_observation) = matrix(1,2);

            matrix = corrcoef(dev_ra, dev_y);
            corr_ra_y_vector(index_observation) = matrix(1,2);
    
            matrix = corrcoef(dev_ra, dev_ra_rn);
            corr_ra_ra_rn_vector(index_observation) = matrix(1,2);
            
            matrix = corrcoef(dev_ra, dev_tb);
            corr_ra_tb_vector(index_observation) = matrix(1,2);
    
            mean_ra_vector(index_observation) = mean(ra_vector);
            mean_ra_vector_rn(index_observation) = mean(ra_vector_rn);
            mean_y_vector(index_observation) = y_vector(sample_size);
            mean_b_vector(index_observation) = b_vector(sample_size);
            mean_cy_vector(index_observation) = tau_vector(sample_size);
            mean_by_vector(index_observation) = by_vector(sample_size);
            mean_tb_vector(index_observation) = tb_vector(sample_size);
            mean_b_adj_y_vector(index_observation) = b_adj_y_vector(sample_size);
            mean_tau_vector(index_observation) = exp(tau_vector(sample_size));
            max_ra_vector(index_observation) = max(ra_vector);

           if def_per_matrix(j,i)+4 < per_num
                b_next_before1(:,index_observation) = b(def_per_matrix(j,i) - 3: def_per_matrix(j,i) + 1,i)*1.01/(delta+.01);
                new_issuance_before(:,index_observation) = (b(def_per_matrix(j,i) - 3: def_per_matrix(j,i) + 1,i) - (1-delta)*(1-d(def_per_matrix(j,i) - 4: def_per_matrix(j,i) ,i)).*b(def_per_matrix(j,i) - 4: def_per_matrix(j,i) ,i))*(1-delta)*1.01/(delta+.01);
                %ra_before1(:,index_observation)     = (((1-delta)./q(def_per_matrix(j,i) - 4: def_per_matrix(j,i),i) +1-delta)/1.01).^4 - 1;
                ra_before1(:,index_observation) = (((1-delta).*q(def_per_matrix(j,i) - 4: def_per_matrix(j,i),i) +1)./(1.01*q(def_per_matrix(j,i) - 4: def_per_matrix(j,i),i))).^4 - 1;
                             
                b_next_after1(:,index_observation)  = b(def_per_matrix(j,i) + 2: def_per_matrix(j,i) + 5,i)*1.01/(delta+.01);
%                ra_after1(:,index_observation)      = (((1-delta)./q(def_per_matrix(j,i) + 1: def_per_matrix(j,i) + 4,i) +1-delta)/1.01).^4 - 1;
                ra_after1(:,index_observation) = (((1-delta).*q(def_per_matrix(j,i) +1: def_per_matrix(j,i)+4,i) +1)./...
                                                 (1.01*q(def_per_matrix(j,i) +1: def_per_matrix(j,i)+4,i)) ).^4 - 1;
                new_issuance_after(:,index_observation)  = (b(def_per_matrix(j,i) + 2: def_per_matrix(j,i) + 5,i) - (1-delta)*b(def_per_matrix(j,i) + 1: def_per_matrix(j,i) + 4,i) ) *(1-delta)*1.01/(delta+.01);
            end 

            index_observation = index_observation+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_simulations = min(500,num_observations);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' \n')
fprintf(['Compute statistics of ',num2str(num_simulations),' samples of ',num2str(sample_size),' observations before a default and without outliers. \n'])
fprintf('std(y)         = %f5 \n',100*mean(std_y_vector(1:num_simulations)))
fprintf('std(c)         = %f5 \n',100*mean(std_c_vector(1:num_simulations)))
fprintf('std(tb)        = %f5 \n',100*mean(std_tb_vector(1:num_simulations)))
fprintf('std(R_s)       = %f5 \n',100*mean(std_ra_vector(1:num_simulations)))
fprintf('corr(y,c)      = %f5 \n',mean(corr_y_c_vector(1:num_simulations)))
fprintf('corr(y,bnext)  = %f5 \n',mean(corr_y_bnext_vector(1:num_simulations)))
fprintf('corr(y,tb)     = %f5 \n',mean(corr_y_tb_vector(1:num_simulations)))
fprintf('corr(y,R_s)    = %f5 \n',mean(corr_ra_y_vector(1:num_simulations)))
fprintf('corr(R_s,tb)   = %f5 \n',mean(corr_ra_tb_vector(1:num_simulations)))
fprintf('Mean debt/y    = %f5 \n',100*mean(mean_by_vector(1:num_simulations)))
fprintf('E(R_s)         = %f5 \n',100*mean(mean_ra_vector(1:num_simulations)))
fprintf('Max R_s        = %f5 \n',100*max(max_ra_vector(1:num_simulations)))



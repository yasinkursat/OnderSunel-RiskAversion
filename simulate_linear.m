clear all

load data_sim.txt
load def_per.txt
load param.txt


per_num = param(1);  %HOW MANY PERIODS IN EACH SAMPLE
n = param(2);        %HOW MANY SAMPLES

y = zeros(per_num, n);
b = zeros(per_num, n);
q = zeros(per_num, n);
qdef = zeros(per_num, n);
c = zeros(per_num, n);
tb = zeros(per_num, n);
e = zeros(per_num, n);
d = zeros(per_num, n);
b_next = zeros(per_num-1, n);

for i=1:n
   y(:,i) = data_sim((i-1)*per_num+1:i*per_num,1); 
   b(:,i) = data_sim((i-1)*per_num+1:i*per_num,2); 
   q(:,i) = data_sim((i-1)*per_num+1:i*per_num,3);
   c(:,i) = data_sim((i-1)*per_num+1:i*per_num,4); 
   tb(:,i) = data_sim((i-1)*per_num+1:i*per_num,5);
   d(:,i) = data_sim((i-1)*per_num+1:i*per_num,6)-1;
   e(:,i) = data_sim((i-1)*per_num+1:i*per_num,7);
end

excl_periods = find(data_sim(:,7)>1);

y_wholesample = data_sim(:,1);% zeros(length(data_sim),1);
c_wholesample = data_sim(:,4);% zeros(length(data_sim),1);
y_wholesample(excl_periods) = [];
c_wholesample(excl_periods) = [];

mean_y = mean(y_wholesample);
mean_output = mean(exp(y_wholesample));
mean_c = mean(c_wholesample);

b_next = b(2:per_num,:);

%CREATE OUTPUT SERIES
num = 500; %HOW MANY PERIODS IN EACH SUBSAMPLE

sample_size = 74;
max_num_def = max(sum(d));
num_observations_vector = zeros(n,1); %NUMBER OF OBSERVATIONS PER SAMPLE (AN OBSERVATION
                               % IS A DEFAULT EPISODE WITH SUFFICIENT
                               % PERIODS BEFORE THE DEFAULT AND NO
                               % EXCLUSION IN BETWEEN
def_per_matrix = zeros(max_num_def,n); %MATRIX OF DEFAULT PERIODS SATISFYING TWO RESTRICTIONS:
                                       %1) NO OTHER DEFAULT OR EXCLUSION OVER THE LAST 75 PERIODS
default_matrix = zeros(max_num_def,n); %MATRIX OF DEFAULT PERIODS 
index_acum = 0; %INDEX OF ACUMULATED NUMBER OF DEFAULT EPISODES
for i = 1:n
    def_num = sum(d(:,i));
    if def_num>0
       default_matrix(1:def_num,i) = def_per(index_acum+1:index_acum+def_num);
       def_periods = def_per(index_acum+1:index_acum+def_num);  %VECTOR OF DEFAULT PERIODS
       interperiods = zeros(def_num,1);
       interperiods(1) = def_periods(1);             %SEPARATION BETWEEN DEFAULT PERIODS
       interperiods(2:def_num)= diff(def_periods);   %SEPARATION BETWEEN DEFAULT PERIODS
       index_acum = index_acum + def_num;            %COUNT DEFAULT PERIODS IN i SAMPLE.
                                         %NEEDED TO KEEP TRACK OF SAMPLE

       %FIND SAMPLES CONTAINING AT LEAST sample_size PERIODS BEFORE THE
       %DEFAULT AND THAT THE LAST EXCLUSION PERIOD HAPPENED sample_size +1
       %PERIODS AGO (NO EXCLUSION IN THE FIRST PERIOD OF THE SAMPLE)
       indices = find(interperiods>sample_size+0 & e(max(def_periods-sample_size-1,1),i)<2);  
       
       num_observations_vector(i) = length(indices);        %STORE THE NUMBER OF SUCH SAMPLES
       if length(indices)>0
           %ONLY STORE RESULTS IF THERE ARE A POSITIVE NUMBER OF
           %OBSERVATIONS.
           def_per_matrix(1:length(indices),i) = def_periods(indices);  
       end
    end
end
num_observations = sum(num_observations_vector);

std_y_vector = zeros(num_observations,1);
std_c_vector = zeros(num_observations,1);
std_tb_vector = zeros(num_observations,1);
std_ra_vector = zeros(num_observations,1);
mean_b_vector = zeros(num_observations,1);
corr_ra_y_vector = zeros(num_observations,1);
corr_y_c_vector = zeros(num_observations,1);
corr_ra_tb_vector = zeros(num_observations,1);
corr_y_tb_vector = zeros(num_observations,1);
corr_ra_c_vector = zeros(num_observations,1);


matrix_y = zeros(sample_size, num_observations);
matrix_ra = zeros(sample_size, num_observations);
matrix_b = zeros(sample_size, num_observations);

index_observation = 1;
time=[ones(sample_size,1)];
for i=1:n
       for j=1:num_observations_vector(i)
            ra_vector = (1./q(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i).^4) -1.017^4;
            y_vector = y(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            tb_vector = tb(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            c_vector = c(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i);
            b_vector = b(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i)./exp(y(def_per_matrix(j,i) - sample_size:def_per_matrix(j,i) - 1,i));
            
            matrix_y(:,index_observation) = y_vector;
            matrix_ra(:,index_observation) = ra_vector;
            matrix_b(:,index_observation) = b_vector;
            lambda=1600;
            hp_y = hpfilter(y_vector,lambda);    %FILTER log 
            hp_tb = hpfilter(tb_vector,lambda);    %FILTER log 
            hp_ra = hpfilter(ra_vector,lambda);    %FILTER log 
            hp_c = hpfilter(c_vector,lambda);    %FILTER log 
            
            dev_y = y_vector - hp_y;
            dev_tb = tb_vector - hp_tb;
            dev_ra = ra_vector - hp_ra;
            dev_c = c_vector - hp_c;
    
            std_ra_vector(index_observation) = std(dev_ra);
            std_c_vector(index_observation) = std(dev_c);
            std_y_vector(index_observation) = std(dev_y);
            std_tb_vector(index_observation) = std(dev_tb);

            matrix = corrcoef(dev_y, dev_c);
            corr_y_c_vector(index_observation) = matrix(1,2);

            matrix = corrcoef(dev_y, dev_tb);
            corr_y_tb_vector(index_observation) = matrix(1,2);

            matrix = corrcoef(dev_ra, dev_y);
            corr_ra_y_vector(index_observation) = matrix(1,2);
    
            matrix = corrcoef(dev_ra, dev_tb);
            corr_ra_tb_vector(index_observation) = matrix(1,2);
            
            matrix = corrcoef(dev_ra, dev_c);
            corr_ra_c_vector(index_observation) = matrix(1,2);
    
            mean_ra_vector(index_observation) = mean(ra_vector);
            mean_b_vector(index_observation) = mean(b_vector);
            max_ra_vector(index_observation) = max(ra_vector);
                      
            index_observation = index_observation+1;            
    end
end

            y_default = y(def_per) - mean_y;
            c_default = c(def_per)- mean_c;
            tb_default = tb(def_per-ones(length(def_per),1));
            b_default = b(def_per)/mean_output;%/exp(y(def_per_matrix(j,i),i));
            ra_default = (1./qdef(def_per).^4) - 1.017^4;


num_simulations = min(length(std_y_vector), 100);

fprintf(' \n')
fprintf(['Compute statistics of ',num2str(num_simulations),' samples of ',num2str(sample_size),' observations before the ast period before a default. \n'])
fprintf('std(y)       = %8.4f \n',100*mean(std_y_vector(1:num_simulations)))
fprintf('std(c)       = %8.4f \n',100*mean(std_c_vector(1:num_simulations)))
fprintf('std(tb)      = %8.4f \n',100*mean(std_tb_vector(1:num_simulations)))
fprintf('std(R_s)     = %8.4f \n',100*mean(std_ra_vector(1:num_simulations)))
fprintf('corr(y,c)    = %8.4f \n',mean(corr_y_c_vector(1:num_simulations)))
fprintf('corr(y,tb)   = %8.4f \n',mean(corr_y_tb_vector(1:num_simulations)))
fprintf('corr(y,R_s)  = %8.4f \n',mean(corr_ra_y_vector(1:num_simulations)))
fprintf('corr(c,R_s)  = %8.4f \n',mean(corr_ra_c_vector(1:num_simulations)))
fprintf('corr(R_s,tb) = %8.4f \n',mean(corr_ra_tb_vector(1:num_simulations)))
fprintf('Mean debt/y  = %8.4f \n',100*mean(mean_b_vector(1:num_simulations)))
fprintf('E(R_s)       = %8.4f \n',100*mean(mean_ra_vector(1:num_simulations)))
fprintf('Max R_s      = %8.4f \n',100*max(max_ra_vector(1:num_simulations)))
fprintf('Default rate   = %8.4f \n',400*sum(sum(d))/(n*per_num))


%COMPUTE STATISTICS IN THE DEFAULT PERIOD.
fprintf(' \n')
fprintf(['Period before a default: average over ',num2str(num_simulations),' default episodes. \n'])
fprintf('y       = %8.4f \n',100* mean(y_default))
fprintf('c       = %8.4f \n',100* mean(c_default))
fprintf('tb      = %8.4f \n', 100* mean(-tb_default))
fprintf('b / y   = %8.4f \n', 100* mean(b_default))
fprintf('Spread  = %8.4f \n',100 * mean(ra_default))
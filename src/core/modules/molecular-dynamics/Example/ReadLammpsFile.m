clear all
close all
clc
ts = 0.001;

total_iter = 2000;

tt = ts*total_iter;

dump_step = 100;

n_iter = total_iter/dump_step;
for k = 1:n_iter

  str = NanoObj;
  str.import('atom.lmp',(k-1)*dump_step,{'Si'});
  
  tmp1 = str.get_property('v_AVETEMP');
  temp_ave(k) = mean(tmp1);
  
  tmp2 = str.get_property('v_mytemp');
  temp_inst(k) = mean(tmp2);
  %temp_inst(k) = mean(temp); 
  
  
end

plot(linspace(0,ts*total_iter,n_iter)/1e3,temp_ave,'r-')
hold on
plot(linspace(0,ts*total_iter,n_iter)/1e3,temp_inst,'b-')
legend('Average','Instantaneus')
xlabel('Time [ns]')
ylabel('Temperature [K]')

Makenice

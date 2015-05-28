load EB.mat
load EBresults.mat
ngenes=length(c1mtEBA);
[void, order] = sort(llr);

nt=length(times2);

printpath='outeb';
printfileall=fullfile(printpath,'all_EB');
genes_vect(1,find(strcmp(genes_vect, 'OCT3/4')))={'OCT3'};
tmu = [0 times2]';

for i = order(end:-1:end-10)'
   
    
  model1=models{i,1};
  model2=models{i,2};
  [mu21, varsigma21] = gpPosteriorMeanVar(model2.comp{1}, tmu);
  [mu22, varsigma22] = gpPosteriorMeanVar(model2.comp{2}, tmu);
  
  rt=mu21-mu22;
  
  times2 = model2.comp{1}.X;
  
  h=figure
  
  %plot(tmu, rt,'b'), hold on
  plot(tmu, mu21, 'g'), hold on
  plot(tmu, mu21+sqrt(varsigma21), 'g:');
  plot(tmu, mu21-sqrt(varsigma21), 'g:');
  plot(tmu, mu22, 'r');
  plot(tmu, mu22+sqrt(varsigma22), 'r:');
  plot(tmu, mu22-sqrt(varsigma22), 'r:');
  plot(model2.comp{1}.X, model2.comp{1}.y, 'k+')
  plot(model2.comp{2}.X, model2.comp{2}.y, 'ko')
  plot(times2(1:7), c1ptEBA(:,i), 'r+')
  plot(times2(8:end), c3ptEBA(:,i), 'r*')
  plot(times2(1:7), c1mtEBA(:,i), 'g+')
  plot(times2(8:end), c3mtEBA(:,i), 'g*')
  
      
  
  dt=sprintf('%s llr=%2.5f',genes_vect{1, i},llr(i));
  title(dt)
  
  
  %legend('mt-pt gp','mt gp','pt gp','c1pt','c3pt','c1mt','c3mt')
  %printfile=fullfile(printpath,genes_vect{1, i});
  %print_str=sprintf('print -depsc2 %s',printfile);
  %print_str_all=sprintf('print -dpsc2 -append %s',printfileall);
  %eval(print_str)
  %eval(print_str_all)
  %close(h)
end

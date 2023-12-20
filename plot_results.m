clc; clear; close all;

% This function run the stochastic_model function and plot the pYAP ratio
% in case of small and large adheisons

%### totalSteps = 
totalSteps = [10]*10^6;

%### diffusion constant (um^2/s )
diffusion = 0.8;

%### number of YAP
particleNumber = 250;

%### number of simulations for each case
iteration = 3; 

%### adhesion properties
adhesionNumber = [1 9];
adhesionSites = [81 81];

%### different cases: [Rb Ra Ru,YAP Ru,pYAP Rdep] (1/s)
rates = {[50 200 0.1 0.1 0.035],[50 200 0.1 0.1 0.035]} ;

%### lattice size (um)
latticeSize = [0.2];

results = struct();
folder = [datestr(now, 30) '_results_test\'];
if ~exist(folder, 'dir')
       mkdir(folder);
end

folder_snap = ['snapshots\'];
if ~exist(folder_snap, 'dir')
       mkdir(folder_snap);
end

%### simulation and plot each case separately

for ls = 1:length(latticeSize) 
    for as = 1:length(adhesionSites)
        results(ls).latticeSize = latticeSize(ls);
        for r = 1:1 
               
                for i = 1:iteration   

                      [results(ls).time{r,as}(i,:), results(ls).pYAPAll{r,as}(i,:), results(ls).pYAPOut{r,as}(i,:)] =  ...
                      stochastic_model(totalSteps(ls), adhesionNumber(as), adhesionSites(as), rates{r,as}, latticeSize(ls), diffusion);

                      results(ls).avgpYAPAll{r,as}(i) = mean(results(ls).pYAPAll{r,as}(i,end-0.4*totalSteps:end));
                      results(ls).avgpYAPOut{r,as}(i) = mean(results(ls).pYAPOut{r,as}(i,end-0.4*totalSteps:end));
                            
                end

                    results(ls).avg_iterations_pYAPAll{r,as} = mean( results(ls).avgpYAPAll{r,as});
                    results(ls).std_iterations_pYAPAll{r,as} = std( results(ls).avgpYAPAll{r,as});

                    results(ls).avg_iterations_pYAPOut{r,as} = mean( results(ls).avgpYAPOut{r,as});
                    results(ls).std_iterations_pYAPOut{r,as} = std( results(ls).avgpYAPOut{r,as});  

                    % plot pYAP all  
                    figure('Position', [10 10 1500 500]); 
                    plot(results(ls).time{r,as}', results(ls).pYAPAll{r,as}'/particleNumber,'LineWidth',2);
                    title(['adhesion-number:' num2str(adhesionNumber(as)) ', adhesion-sites:' num2str(adhesionSites(as))]);
                    xlabel('Time (s)'); ylabel('pYAP ratio');
                    set(gca,'Fontname', 'Times New Roman','fontsize',14);
                    %saveas(gcf,[folder '_pYAPAll' '_lattice_size_' num2str(latticeSize(ls))  '_adhesion_number_' num2str(adhesionNumber(as)) '_Rb_' num2str(rates{r,as}(1)) '_Rdep_' num2str(rates{r,as}(end)) '.png']);        
                                   
        end
    end
end

%save([folder 'sspYAP'], '-v7.3');


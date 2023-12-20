
function [T, activatedNumber, activatedNumber_outside] = stochastic_model(totalSteps, adhesionNumber, adhesionSites, rates, latticeSize, diffusion) 

%### Inputs:
% totalSteps: total number of iterations
% adhesionNumber: number of adhesions
% adhesionSites: number of binding sites
% rates: rates for diffusion and reation events
% latticeSize: lattice size
% diffusion: diffusion constant 

%### Outputs:
% T is time
% activatedNumber: total number of pYAP
% activatedNumber_outside: total number of pYAP outside the adhesions

% Start Timing
tic

%### Define the "physical" system
dimensions = 3; boxSize = [40 40 20];

%### system parameters
particleNumber = 250; 
initialpYAP = 0; 

%### Define constant rates 
diffusionConstant = diffusion; % diffusion rate of YAP/pYAP
rateDiffusion = (2*dimensions*diffusionConstant)/latticeSize^2; 
rateBinding = rates(1);        % Rb: YAP binding rate 
rateActivation = rates(2);     % Ra: YAP phosphorylation rate
rateUnbinding = rates(3);      % Ru,YAP: YAP-adhesion unbinding rate 
rateUnbinding2 = rates(4);     % Ru,pYAP: pYAP-adhesion unbinding rate
rateDeactivation = rates(5);   % Rdep: dephosphorylation rate 

%### initial Random positions  
domain = zeros(boxSize(1),boxSize(2),boxSize(3));

%### particles random positions (no overlap)
particleRandomPosition  = [0 0 0; 0 0 0];
while any(pdist(particleRandomPosition) < 0.1)
    particleRandomPosition = ceil((boxSize(1)-4)*rand(particleNumber,2)) + 2;        
    particleRandomPosition(:,3) = ceil((boxSize(3)-4)*rand(particleNumber,1)) + 2;  
end
for m = 1:length(particleRandomPosition)
    domain(particleRandomPosition(m,1),particleRandomPosition(m,2),particleRandomPosition(m,3)) = 1;
end

%### initial adhesion positions on the membrane
[integrinPosition,adhesionCMX,adhesionCMY] = integrinDistribution(boxSize, adhesionSites, adhesionNumber); %adhesionSize --> integrinNumber

%### Generate Position Arrays
[particlePosition(:,1), particlePosition(:,2), particlePosition(:,3)] = ind2sub(boxSize,find(domain));
particleIndex = sub2ind( boxSize, particlePosition(:,1), particlePosition(:,2), particlePosition(:,3));

%### activation state of particles (YAP)
particleState = zeros(particleNumber,1); 
particleState(randperm(particleNumber,initialpYAP)) = 1;

%### possible directions
dR = {[-1 0 0],[+1 0 0],[0 -1 0],[0 +1 0],[0 0 -1],[0 0 +1]}; %3D     

%### number of pYAP 
activatedNumber = zeros(totalSteps,1);
activatedNumber_outside = zeros(totalSteps,1);

%### no pYAP in the beginning
T(1) = 0; 
activatedNumber(1) = length(find(particleState == 1));
activatedNumber_outside(1) = length( find(particleState == 1 & particlePosition(:,3) ~= 1 ) );

%### generating first set of random numbers are generated 
% for the first step the only possible event is diffusion 

tau = -log(rand(1,particleNumber))./rateDiffusion;        
[sortedTau, indexTau] = sort(tau);

%### time stepping

for step = 1:totalSteps  

    T(step+1) = sortedTau(1);   
    particle = indexTau(1);

    %### decide which reaction, and update  
    [particlePosition, particleIndex, particleState] = ...
        newEvent(particle,particlePosition,particleIndex,particleState,dR,dimensions,boxSize);

    sortedTau(1) = nextEvent(particle,particlePosition,particleState,particleIndex,sortedTau);

    [sortedTau, indexTauNew] = sort(sortedTau);
    indexTau = indexTau(indexTauNew);
    
    % Update Domain Image
    domain = zeros(boxSize(1),boxSize(2),boxSize(3));
    domain(particleIndex) = 1;
    
    if  mod(step,200000) == 0  
        %### plot the position every 200000 steps
         %plotPosition(particlePosition, particleState, integrinPosition, boxSize, step, T);
    end

    %### print to screen 
    if mod(step,10000000) == 0
        fprintf(1,'Step is %d \n', step )
    end
    
    activatedNumber(step+1) = length(find(particleState == 1));
    activatedNumber_outside (step+1) = length( find(particleState == 1 & particlePosition(:,3) ~= 1 ) );

end

toc

%### periodic boundary conditions
function [outputPosition, moveIndex] = PeriodicBoundary(tempPosition,boxSize)

    if tempPosition(3) <= boxSize(3) && tempPosition(3) >= 2
        tempPosition(tempPosition>boxSize(1)) = tempPosition(tempPosition>boxSize(1)) - boxSize(1);
        tempPosition(tempPosition<1) = tempPosition(tempPosition<1) + boxSize(1);
    else
        if tempPosition(3) < 2  %particles cannot diffuse to membrane where there is no adhesion 
            tempPosition(3) = tempPosition(3) + 1;
        elseif tempPosition(3) > boxSize(3)
            tempPosition(3) = tempPosition(3) - 1;
        end
    end

    outputPosition = tempPosition;
    % Convert to indicies for speedy search
    moveIndex = sub2ind( boxSize, tempPosition(1), tempPosition(2), tempPosition(3));

end

%### plot the positions
function plotPosition(particlePosition, particleState, integrinPosition, boxSize, step, T)

figure('Position', [10 10 700 700]); %'visible','off'

X = particlePosition(particleState == 0, 1);
Y = particlePosition(particleState == 0, 2);
Z = particlePosition(particleState == 0, 3);
scatter3(X,Y,Z, 14, 'b','filled'); hold on;

X = particlePosition(particleState == 1, 1);
Y = particlePosition(particleState == 1, 2);
Z = particlePosition(particleState == 1, 3);
scatter3(X,Y,Z, 164,'pentagram', 'r','filled'); hold on;

[X,Y] = meshgrid(1:boxSize(1),1:boxSize(2));
scatter3(X,Y,ones(boxSize(1),boxSize(2)),6,'MarkerEdgeColor',[220 220 220]/255,'MarkerFaceColor',[220 220 220]/255); hold on;

scatter3(integrinPosition(:,1), integrinPosition(:,2), integrinPosition(:,3), 28,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 1 0]); 

view(40,35)

xlim([0, boxSize(1)]); ylim([0, boxSize(2)]); zlim([1, boxSize(3)]);
title(['Step = ' num2str(step) ' Time = ' num2str(T(end))]); xlabel('X'); ylabel('Y'); zlabel('Z');
set(gca,'Fontname', 'Times New Roman','fontsize',14);   
saveas(gcf,['./snapshots/step_' num2str(step) '.png']);
    
end

%### random distribution of integrins on the membrane
function [integrinPosition, adhesionCMX, adhesionCMY] = integrinDistribution(boxSize, adhesionSites, adhesionNumber) 

integrinPeradhesion = adhesionSites/adhesionNumber; 
integrinPosition = ones(adhesionSites,4); 
d = (ceil(sqrt(integrinPeradhesion)) - 1)/2; 
r = 0;

attempt_counter = 1 ; 
Distances = 0;
    while any(Distances < 2*ceil(d)+1)
        attempt_counter = attempt_counter + 1; 
        adhesionCMX = ceil((boxSize(1)-2*ceil(d))*rand(adhesionNumber,1)) + ceil(d) ; %randi([1+d boxSize(1)-d], adhesionNumber, 1);
        Distances = pdist(adhesionCMX);
    end    

attempt_counter = 1 ; 
Distances = 0;
    while any(Distances < 2*ceil(d)+1)
        attempt_counter = attempt_counter + 1; 
        adhesionCMY = ceil((boxSize(1)-2*ceil(d))*rand(adhesionNumber,1)) + ceil(d) ; %randi([1+d boxSize(2)-d], adhesionNumber, 1);
        Distances = pdist(adhesionCMY);
    end  

range = -ceil(d):ceil(d);

    for n = 1:adhesionNumber        
        for j = range
            for l = range
                r=r+1;
                if r <= adhesionSites
                integrinPosition(r,:) = [adhesionCMX(n)+l adhesionCMY(n)+j 1 n];
                % use the upper squred number
                end
                if  mod(r,ceil(adhesionSites/adhesionNumber)) == 0 
                break
                end
            end
            if  mod(r,ceil(adhesionSites/adhesionNumber)) == 0 
                break
            end
        end
        
     end 

end

%### diffusion
function [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize)

        randDirection = ceil( 2 * dimensions * rand );
        tempPosition = particlePosition(particle,:) + dR{randDirection};
        [tempPosition, moveIndex] = PeriodicBoundary(tempPosition,boxSize);
            
        % If where we want to move is empty, move
        if ~sum(ones(length(particleIndex),1)*moveIndex == particleIndex)
            particlePosition(particle,:) = tempPosition;
            particleIndex(particle) = moveIndex;
        end

end

%### updating YAP state/position based on possible events
function [particlePosition, particleIndex, particleState] = ...
        newEvent(particle,particlePosition,particleIndex,particleState,dR,dimensions,boxSize)

    if particlePosition(particle,3) > 2
    
        if particleState(particle) == 0
           [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize);
            
        else
           sumRates = rateDiffusion + rateDeactivation;  
           outputEvent = randsrc(1,1,[1 2; rateDiffusion/sumRates rateDeactivation/sumRates]); 
           switch outputEvent
               case 1
             [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize);            
               case 2 
             particleState(particle) = 0;    
           end
        end
    
    elseif particlePosition(particle,3) == 2
    
          if  particleState(particle) == 0      
              if  intersect( find(particlePosition(particle,1) == integrinPosition(:,1)), ...
                             find(particlePosition(particle,2) == integrinPosition(:,2))) 
                  
                  sumRates = rateDiffusion + rateBinding;
                  outputEvent = randsrc(1,1,[1 2; rateDiffusion/sumRates rateBinding/sumRates]); 
                  switch outputEvent
                      case 1
                     [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize);                    
                      case 2   %here we also check if the binding site on integrin is not occupied
                                moveIndex = sub2ind( boxSize, particlePosition(particle,1), particlePosition(particle,2), 1);
    
                            if  ~sum(ones(length(particleIndex),1)*moveIndex == particleIndex)
                                particlePosition(particle,3) = 1;  
                                particleIndex(particle) = moveIndex;
                            end
                  end
              else                        
                 [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize);
    
              end
    
          else                           
              sumRates = rateDiffusion + rateDeactivation;
              outputEvent = randsrc(1,1,[1 2; rateDiffusion/sumRates rateDeactivation/sumRates]);
              switch outputEvent
                  case 1
                 [particlePosition, particleIndex] = diffusionEvent(particle,particlePosition,particleIndex,dR,dimensions,boxSize);
                  case 2
                 particleState(particle) = 0; 
              end
          end
    
    elseif particlePosition(particle,3) == 1

        if particleState(particle) == 0          
            sumRates = rateActivation + rateUnbinding;
            outputEvent = randsrc(1,1,[1 2; rateActivation/sumRates rateUnbinding/sumRates]);
            switch outputEvent
                  case 1   
                   particleState(particle) = 1;
                  case 2
                   particlePosition(particle,3) = 2;    
                   moveIndex = sub2ind( boxSize, particlePosition(particle,1), particlePosition(particle,2), particlePosition(particle,3));
                   particleIndex(particle) = moveIndex;
            end       
        else                           
                   particlePosition(particle,3) = 2;      
                   moveIndex = sub2ind( boxSize, particlePosition(particle,1), particlePosition(particle,2), particlePosition(particle,3));
                   particleIndex(particle) = moveIndex;       
        end  

    end    
end

%### calculating time to the next even
function nextEventTau = nextEvent(particle,particlePosition,particleState,particleIndex,sortedTau)
    if particlePosition(particle,3) > 2

        if particleState(particle) == 0
           nextEventTau = -log(rand)/rateDiffusion + sortedTau(1); 
        else
           sumRates = rateDiffusion + rateDeactivation; 
           nextEventTau = -log(rand)/sumRates + sortedTau(1);
        end

    elseif particlePosition(particle,3) == 2

          if  particleState(particle) == 0       
              if  intersect( find(particlePosition(particle,1) == integrinPosition(:,1)), ...
                             find(particlePosition(particle,2) == integrinPosition(:,2)))             
                  sumRates = rateDiffusion + rateBinding;
                  nextEventTau = -log(rand)/sumRates + sortedTau(1);
              else                       
                  nextEventTau = -log(rand)/rateDiffusion + sortedTau(1);  
              end
          else                          
              sumRates = rateDiffusion + rateDeactivation;
              nextEventTau = -log(rand)/sumRates + sortedTau(1); 
          end

    elseif particlePosition(particle,3) == 1     
        if particleState(particle) == 0          
            sumRates = rateActivation + rateUnbinding;
            nextEventTau = -log(rand)/sumRates + sortedTau(1);   
        else                            
            sumRates = rateUnbinding2; 
            nextEventTau = -log(rand)/sumRates + sortedTau(1);
        end
    end

end

end

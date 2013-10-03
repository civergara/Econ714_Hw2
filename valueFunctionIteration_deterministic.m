function[value]=valueFunctionIteration_deterministic( vGridCapital , nGridLabor, tolerance )

%%  1. Calibration

aalpha = 1/3;           % Elasticity of output w.r.t. capital
bbeta  = 0.95;          % Discount factor
ddelta = .09;          % Depreciation
laborSteadyState=1/3;   % Labor in steady state

% Productivity values
vProductivity = [0.9792; 0.9896; 1.0000; 1.0106; 1.0212]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
                 0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
                 0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
                 0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
                 0.0000, 0.0000, 0.0000, 0.0273, 0.9727];
             
%% 2. Steady State

capitalSteadyState = ((1/bbeta-1+ddelta)/(aalpha*laborSteadyState^(1-aalpha)))^(1/(aalpha-1));
outputSteadyState = capitalSteadyState^aalpha*laborSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-capitalSteadyState*ddelta;

psi=outputSteadyState*(1-aalpha)/((outputSteadyState-ddelta*capitalSteadyState)*laborSteadyState^2);

fprintf(' Output = %2.6f, Capital = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState); 
fprintf('\n')
fprintf(' Value of Leisure=%2.6f\n',psi)
fprintf('\n')

% We generate the grid of capital

vGridLabor=linspace(0.5*laborSteadyState,1.5*laborSteadyState,nGridLabor);

nGridCapital=length(vGridCapital);
nGridProductivity = length(vProductivity);

%% 2. Required matrices and vectors for First Iteration

mOutput           = zeros(nGridCapital,nGridLabor);
mValueFunction    = ones(nGridCapital,1)*(log(consumptionSteadyState)-psi*laborSteadyState^2/2);
mValueFunctionNew = zeros(nGridCapital,1);
mPolicyFunction   = zeros(nGridCapital,1);
mLaborFunction    = zeros(nGridCapital,1);


%% 3. We pre-build output for each point in the grid for First Iteration

for nCapital=1:nGridCapital
    for nLabor=1:nGridLabor
            mOutput(nCapital,nLabor) = (vGridCapital(nCapital))^aalpha*(vGridLabor(nLabor))^(1-aalpha)+...
                vGridCapital(nCapital)*(1-ddelta);
        end
end

%% 4. First Iteration

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

while (maxDifference>tolerance)  
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            gridLabor=1;
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                valueHighSoFar1=-1000.0;
                
                for nLabor = gridLabor:nGridCapital
                    
                consumption = mOutput(nCapital,nLabor)-vGridCapital(nCapitalNextPeriod);
                valueProvisional =(1-bbeta)*(log(consumption)-psi*vGridLabor(nLabor)^2/2)+...
                bbeta*mValueFunction(nCapitalNextPeriod);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    gridLabor=nLabor;
                    laborChoice=vGridLabor(nLabor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital) = valueHighSoFar;
            mPolicyFunction(nCapital) = capitalChoice;
            mLaborFunction(nCapital)=laborChoice;
        end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')

fprintf(' My chek = %2.6f\n', mPolicyFunction(10)); 
fprintf('\n')

%% 5. Required matrices and vectors for Main Iteration

mValueFunction1=zeros(nGridCapital,nGridProductivity);

for nCapital=1:nGridCapital
    for nProductivity=1:nGridProductivity
        mValueFunction1(nCapital,:)=mValueFunction(nCapital);
    end
end

mValueFunction=mValueFunction1;
clear mValueFunction1 capitalChoice consumption gridCapitalNextPeriod gridLabor iteration...
    laborChoice mLaborFunction mOutput mPolicyFunction mValueFunctionNew maxDifference nCapital...
    nCapitalNextPeriod tolerance valueHighSoFar valueHighSoFar1 valueProvisional nLabor nProductivity;

mOutput           = zeros(nGridCapital,nGridLabor,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunction   = zeros(nGridCapital,nGridProductivity);
mLaborFunction    = zeros(nGridCapital,nGridProductivity);
 
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 6. We pre-build output for each point in the grid for Main Iteration

for nCapital=1:nGridCapital
    for nLabor=1:nGridLabor
        for nProductivity=1:nGridProductivity
            mOutput(nCapital,nLabor,nProductivity) = (vGridCapital(nCapital))^aalpha*...
                vProductivity(nProductivity)*(vGridLabor(nLabor))^(1-aalpha);
        end
    end
end

%% 7. Main iteration

maxDifference = 10.0;
iteration = 0;

%while (maxDifference>tolerance)  
while(iteration<240)   
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for nCapital = 1:nGridCapital
                        
            valueHighSoFar = -1000.0;
            capitalChoice  = vGridCapital(1);
            gridLabor=1;
            
            for nCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                valueHighSoFar1=-1000.0;
                
                for nLabor = gridLabor:nGridLabor
                    
                consumption = mOutput(nCapital,nLabor,nProductivity)-vGridCapital(nCapitalNextPeriod)+vGridCapital(nCapital)*(1-ddelta);
                valueProvisional =(1-bbeta)*(log(consumption)-psi*vGridLabor(nLabor)^2/2)+...
                bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);              
                
                if (valueProvisional>valueHighSoFar1)
                    valueHighSoFar1=valueProvisional;
                    gridLabor=nLabor;
                    laborChoice=vGridLabor(nLabor);
                else
                    break;
                end
                
                end
                      
                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    capitalChoice = vGridCapital(nCapitalNextPeriod);
                    gridCapitalNextPeriod = nCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    
                  
            end
            
            mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar;
            mPolicyFunction(nCapital,nProductivity) = capitalChoice;
            mLaborFunction(nCapital,nProductivity)=laborChoice;
        end
        
    end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')

fprintf(' My chek = %2.6f\n', mPolicyFunction(10,3)); 
fprintf('\n')

toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy Function')

subplot(3,1,3)
plot(vGridCapital,mLaborFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Labor Function')

value=mValueFunction;

end

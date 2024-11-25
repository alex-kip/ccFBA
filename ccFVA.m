function [carbonConst] = ccFVA(model,bnds,condition,carbonCount,rlxBnds,rlxInt,cofPairs,minmax,runInitMM,runFinalMinMax,limEX,totalCarbon,skipRxns,CBBRxns,delimiter,SearchProb_rxns)
%% This script is designed to constrain the solution space based on the 
%   number of carbon atoms made available to the model by the exchange
%   reactions. This is done in order to avoid physiologically meaningless
%   futile cycles.
%
%   CREATED 2017-09-30 by Maximilian Lularevic
%
% INPUT
%   model       =   cobra model (must contain model.metFormulas)
%
%   
%   
%                   
% OPTIONAL
%   
%   bnds        =   bounds or constraints as cell or double. for double 
%                   first row corresponds to rxnIDs, second row to lower 
%                   and the third row to upper bnds. These bounds should 
%                   only contain are only uptake and secretion rates.
%                   Otherwise intracellular constraints are potentially 
%                   taken into the calculation of totalCarbon uptaken
%
%   condition   =   string containing information as to which condition is 
%                   run (e.g. lactate producing conditon)
%
%   carbonCount =   carbonCount Structure as received form function 
%                   findMetCarbon.m (needs to contain
%                   carbonCount.carbonAtoms field). If this is not passed
%                   into the function it will create said structure by
%                   running the findMetCarbon function.
%
%   relaxBnds   =   If true and carbon constraints do not solve due to too 
%                   tight constraintsalgorithm will increase total carbon 
%                   consumed based on relaxInt (relaxation interval) which 
%                   is a percentage. If flase algorithm will throw an error
%                   warning the user that the model does not solve under
%                   the carbon constraints and exit the script. Default =
%                   true
%
%   relaxInt    =   percentage step totalCarbon will be increased/ relaxed
%                   in case the model does not solve after carbon
%                   constraining. Default = 0.1 (i.e. 10%)
%
%   cofParis    =   cofParis structure as received from function
%                   findCofactorPairs.m. If not passed into the function it
%                   will create said structure by running cofPairs
%                   function.
%
%   minmax      =   FVA for model under experimental constraints. (needs 
%                   to be processed with fixMinMax.m). If not passed into
%                   the function it will create said structure by running
%                   runMinMax.m (alternative to fluxvariability.m). 
%
%   runInitMM   =   Logical input true or false (1 or 0). Default = 1. This
%                   variable get only evaluated if minmax variable is empty
%                   or not input into function. If no initial FVA was
%                   input the experimenter has the option to perform a
%                   FVA with the runMinMax.m function (supplied with the
%                   toolbox). If not the the original bounds of the model
%                   (lb and ub vecotrs will be used to creat the minmax
%                   variable).
%
%   runFinalMM  =   logical true or false. if true a final MinMax (FVA)
%                   will be run
%
%   limEX       =   value to be used to constrain drain rxns. This is
%                   supposed to be a single carbon constraint (e.g. left
%                   over carbon can be set to CO2 rate which is a single
%                   carbon. If flux for CO2 (closing the carbon by mass balance)
%                   is 0.5 mmol/gDCW*h this would be the rate to feed into
%                   the limEX variable).
%
%   totalCarbon =   total single carbon flux over-write. This value is
%                   generally calculated based on the input bounds for
%                   exchange reactions but can be input here manually. 
%
%
%   skipRxns    =   vector of doubles that denote reactions that are to be
%                   skipped for carbon constraining (could be ussed for
%                   reactions without a chemical formula for example or the
%                   biomass reaction (skipped by deafault))
%  
%   CBBRxns     =   Extension for Calvin Cycle under autotrophic conditions. 
%                   Vector of doubles that denote reactions involved in the
%                   Calvin Cycle.
%
%   delimiter   =   This variables default is '_' meaning metaboites in the
%                   input model are contain the compartment information
%                   after the underscroe (e.g. glc_D_c, or nadh_m). Often
%                   models come with other delimiter such as square
%                   brackets (e.g. glc_D[c], nadh[m]) and therfore this 
%                   input should be carefully checked.
%        
% SearchProb_rxns = This Variables default is true, meaning that IF there
%                   is infeasibility or reduced Objective Value created by the post_CCFVA bounds,
%                   it searches for possible Rxns that create the problem.
%
% 
% OUTPUT
%   carbonConst =   structure containing the following fields
%                       - Title
%                       - Date and Time of job
%                       - constraintsName (Name of bnds variable: e.g. 'LacProd')
%                       - model.description
%                       - total carbon consumed
%               *** optional *** (start)
%                       - actual carbon used for constraining (only when
%                         passed into function)
%               *** optional *** (finish)
%                       - carbonRelaxed
%               *** optional *** (start)
%                       - rlxVal
%                       - percentageRelaxed
%                       - relaxationIterations
%                       - relaxedTotalCarbonConsumed
%               *** optional *** (finish)
%                       - report
%                           - cell with the lenght of model.rxns which 
%                             states whther reaction was balanced or unbalanced  
%                           - number of balanced rxns
%                           - number of unbalanced rxns
%                           - number of set rxns (aka constraints)
%                           - number of rxns with no flux
%                       - bndAdjusted
%                           - rawData (logical 2 colum matrix)
%                           - number of lower bounds adjusted by cc
%                             algorithm
%                           - number of upper bnds adjusted by cc algorithm
%                       - carbonCount
%                           - dateAndTime
%                           - modelDescription
%                           - metName
%                           - metFormula
%                           - carbon atoms for each metabolite
%                           - table with all the above information compiled
%                       - cofactorPairs (cell with mets and met IDs and
%                         occurance of these mets with corresponding
%                         rxnIDs)
%                       - original minmax (either the one passed into the
%                         function or herein calculated one)
%                       - newMinMax passed on carbon balance
%                       - solSpaceReduction size of new solution space
%                         compared to old one (calculatons based on minmax
%                         range comparison)
%                       - solutionNew (optimizeCbModel output with new bnds)
%                       - EXrxnsAdjusted -> title
%                       - limEX (single carbon limit passed into fucntion
%                         to constrain the Exchange reactions output)
%                       - number of EX bnds adjusted
%                       - solSpaceReductionEX in % compared to original
%                       - Problematic_rxns: Reactions that may drop the
%                         objective value when constrained
%                       - EXspaceReduction (how much has the space shrunk
%                         just comparing old EX ranges with new ones
%                       - post New MinMax -> title
%                       - postNewMinMax (new minmax run with runMinMax
%                         and cc bnds)
%                       - solSpacereductionPostCC
%                       - solutionPostNew
%
%                            
%                           -----------------------------------------------------------
%                            This code is publicly available and is the
%                            property of the Laboratory of Biochemical
%                            Engineering©  https://lbe.cheng.auth.gr/
%                            
%
% CHANGE LOG
%
%   Date:       Comment:                            Change made by
%
%   28-Feb-19   Added automatic relaxation          Maximilian Lularevic
%               option for totalCarbon and
%               changed function name to
%               ccFVA from carbonConstraints
%   
%   01-Mar-19   Added runInitMM variable and        Maximilian Lularevic
%               adjusted the evaluation section
%               of minmax variable where if 
%               runInitMM is false the initial
%               minmax variable will be set to
%               the mdels lower and upper bounds
%
%   24-Mar-19   Addtion of delimiter variable       Maximilian Lularevic
%               was added to the function
%
%   Jul-2022    Addition of Calvin Cycle            Giannis Vourvachakis
%               constraint methodology.
%
%   July-2023   Change in findLrgMlcs. Now the
%               function can find also antiport
%               reactions.                          Nikolaos Stratis
%

%% check if minmax was passed into function and, run minmax if it hasn't


if exist('runInitMM','var')
    if isempty(runInitMM)
        runInitMM = 1;
    end
else
    runInitMM = 1;
end

[~,num_rxns] = size(model.S);
% Prepare model for ccFVA (add the given bounds and find directionality of
% the Rxns
[~,model] = PrepareModel(model);

% if Variable bnds is passed into the function  and its empty then find
% the bnds

if exist('bnds','var')
    if isempty(bnds)
     [bnds,model] = PrepareModel(model);    
    end
else 
     [bnds,model] = PrepareModel(model);   
end


% load those constraints 
model = loadConstraints(model,bnds,'FBA');
if exist('minmax','var')
    if isempty(minmax)
        if runInitMM
            minmax = runMinMax(model);
            minmax = fixMinMax(minmax);
        else
            minmax = [model.lb,model.ub];
        end
    end
else
    if runInitMM
        minmax = runMinMax(model);
        minmax = fixMinMax(minmax);
    else
        minmax = [model.lb,model.ub];
    end
end

%% check if condition was passed into function
if exist('condition','var')
    if isempty(condition)
        condition = 'n/a';
    end
else
    condition = 'n/a';
end

% check relaxBnds was passed into function
if exist('rlxBnds','var')
    if isempty(rlxBnds)
        rlxBnds = 1;
    end 
else
    rlxBnds = 1;
end

% check relaxBnds was passed into function
if exist('rlxInt','var')
    if isempty(rlxInt)
        rlxInt = 0.1;
    end 
else
    rlxInt = 0.1;
end

% set rlxVal
rlxVal = 1; % 1 = 100%

% check if final minmax should be run
if exist('runFinalMinMax','var')
    if isempty(runFinalMinMax)
        runFinalMinMax = 1;
    end 
else
    runFinalMinMax = 1;
end

% check if there are any rxns to be skipped
if exist('skipRxns','var')
    if isempty(skipRxns)
        skipRxns = find(model.c ==1);  %skip the biomass equation
    else
        [s,z] = size(skipRxns);
        if z > s
            skipRxns = skipRxns';
        end
    end
else
    skipRxns = find(model.c == 1);  
end
if exist('CBBRxns','var')
   if isempty(CBBRxns)
        CBBRxns = [];
    else
        [s,z] = size(CBBRxns);
        if z > s
            CBBRxns = CBBRxns';
        end
    end
else
    CBBRxns = [];
end  
% check if delimiter variable was put into function
if exist('delimiter','var')
    if isempty(delimiter)
        delimiter = '_';
    end
else
    delimiter = '_';
end

% check if SearchProb_rxns was passed into function
if exist('SearchProb_rxns','var')
    if isempty(SearchProb_rxns)
        SearchProb_rxns = 1;
    end 
else
    SearchProb_rxns = 1;
end

%% extract lower and upper bnds from bnds variable
if iscell(bnds)
    k = 0;
    [num_vars,~] = size(bnds);
    for i=1:num_vars
        varName = bnds{i,1};
        var_index = find(ismember(model.rxns,varName));
        if isempty(var_index)
            fprintf('%s not found\n',bnds{i,1});
        else
            k = k+1;
            varID(k,1) = var_index;
            varMin(k,1) = bnds{i,2};
            varMax(k,1) = bnds{i,3};
        end
    end
elseif isnumeric(bnds)
    varID = bnds(:,1);
    varMin = bnds(:,2);
    varMax = bnds(:,3);   
end
varID_tmp = varID;
%% find corresponding metIDs from varIDs
% first remove all bnd rxn that are not true drains (e.g. biomass, IgG formation)
isDrain = find(checkDrainRxns(model));
noExRxns = find(~ismember(varID,isDrain));
% adjust the previously retrieved variables
varID(noExRxns) = [];
varMin(noExRxns) = [];
varMax(noExRxns) = [];
% get all drains without the set bnds
isDrainNew = isDrain(~ismember(isDrain,varID));

%% check if number of carbon atoms per metabolite was passed into function
if exist('carbonCount','var') 
    if isempty(carbonCount)
        carbonCount = findMetCarbon(model);
    end
else
    carbonCount = findMetCarbon(model);
end
% extract the carbon vector (length of model.mets)

metCarb = carbonCount.carbonAtoms;

%% multiply the lower and upper bnds witht their respective number of carbon
% atoms
for i = 1:size(varID,1)
    metID = find(model.S(:,varID(i,1)));
    atomCount(i,1) = sum(carbonCount.carbonAtoms(metID,1));
end
varMinCarb = varMin.*atomCount;
% varMaxCarb = varMax.*metCarb(metID);
if exist('totalCarbon','var')
    if isempty(totalCarbon)
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    elseif ~isnumeric(totalCarbon)
        h = warndlg('Warning the totalCarbon input is not numeric. The totalCarbon will be calculated based on the lower bnds of the input variables. If you wish to continue press OK or else kill the code and input a numeric value for totalCarbon.','Problem with input');
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    else
        total_carbon_cons = totalCarbon;
        actual_total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    end
else
    total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    intotal_carbon_cons=total_carbon_cons;
end
total_carbon_cons=total_carbon_cons;
intotal_carbon_cons=total_carbon_cons; %initial total carbon consumed (might be followed by relaxation)
%% check if cofactorPairs have been passed on 

% check for cofPair variable and/or calculate new
if exist('cofPairs','var')
    if isempty(cofPairs)
        smallMets = findSmallMetsCarb(model);
        metPairs = findMPs(model,smallMets.smallMetIDs,false,delimiter);
        cofPairs = findCofactorPairs(model,metPairs);
    end
else
    smallMets = findSmallMetsCarb(model);
    metPairs = findMPs(model,smallMets.smallMetIDs,false,delimiter);
    cofPairs = findCofactorPairs(model,metPairs);
end

% initialize variable
cofMetIDs = [];
% retrieve all the metabolite IDs for cofactors
for i = 1:size(cofPairs,1)
    cof_tmp = cell2mat(cofPairs{i,6});
    cofMetIDs = [cofMetIDs;cof_tmp];
end
% retrieve all the metabolite IDs for cofactors
for k = 1:size(cofPairs, 1)
                          cofPairs{k, 5} = [cofPairs{k, 5}{:}];
                          cofPairs{k, 6} = [cofPairs{k, 6}{:}];
end

% set cofMets to 0 in model.S
model_tmp = model;
model_tmp.S(cofMetIDs,:) = 0;

% find all rxns containing coa mets on both sides of the rxn
coa = findCoAs(model_tmp);

% find all rxns containing cytochromes, quinones, and hemes
lrgMlc = findLrgMlcs(model,carbonCount);

%% Initial result structure 
carbonConst.title = '*** ccFVA ***';
carbonConst.dateAndTime = datetime();
carbonConst.constraintsName = condition;
if isfield(model,'description')
    carbonConst.model = model.description;
end
carbonConst.totalCarbonConsumed = total_carbon_cons;
if exist('actual_total_carbon_cons','var')
    if ~isempty(actual_total_carbon_cons)
        carbonConst.actualTotalCarbonConsumed = actual_total_carbon_cons;
    end
end

%% calculating new bnds algorithm
minmax_new = minmax;
mmAdj = zeros(num_rxns,2);
reRun = 1;
% counter = 0;
while reRun
    for i = 1:num_rxns
        % evaluate if rxn is part of the measured constraints 
        % rxns
        if ismember(i,varID_tmp)
            bericht{i,1} = 'set_bnd';
            a{i,1}='exchange';
            continue
        elseif ismember(i,skipRxns)
            bericht{i,1} = 'skipped rxn';
             a{i,1}='skipped';
            continue
        else %intracellurar reactions
             % check if reaction carries flux under normal conditions
            if sum(abs(minmax(i,:))) == 0
                bericht{i,1} = 'no_flux';
                a{i,1}='no flux';
                continue
            elseif ismember(i,isDrain)
                bericht{i,1} = 'balanced - drain';
                subIDs = find(model.S(:,i) < 0);
                prodIDs = find(model.S(:,i) > 0);

                % account for stoichiometric coefficients and add up carbon
                subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
                prodCarbon = sum(abs(full(model.S(prodIDs,i))).*metCarb(prodIDs,1));

                if subCarbon == 0 && prodCarbon == 0
                    bericht{i,1} = 'carbon unconstrained';
                    continue
                elseif subCarbon == 0
                    constraint = abs(total_carbon_cons/prodCarbon);
                else
                    constraint = abs(total_carbon_cons/subCarbon);
                end
                if model.qualDir(1,i) == 0 %bidirectional 
                        if abs(minmax(i,1)) > constraint
                            minmax_new(i,1) = -constraint;
                            mmAdj(i,1) = 1;
                        end
                        if abs(minmax(i,2)) > constraint
                           minmax_new(i,2) = constraint;
                           mmAdj(i,2) = 1;
                        end
                    
                elseif model.qualDir(1,i) == 1 
                    %unidirectional forward [0 - 100]
                           if abs(minmax(i,2)) > constraint
                            if constraint > abs(minmax(i,1))
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                          end
                    
                
                elseif model.qualDir(1,i) == -1 
                    %unidirectional backwords [-100 - 0]
                        if abs(minmax(i,1)) > constraint
                            if constraint > abs(minmax(i,2))
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                        end
                 end
            else
                % find substrates and products 
                subIDs = find(model.S(:,i) < 0);
                prodIDs = find(model.S(:,i) > 0);

                % account for stoichiometric coefficients and add up carbon
                subCarbon(i) = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
                prodCarbon(i) = sum(abs(full(model.S(prodIDs,i))).*metCarb(prodIDs,1));


                          %detection of cofactor pairs in the reaction
                           found = false(size(cofPairs, 1), 1);
                           for k = 1:size(cofPairs, 1)
                               found(k) = any(cofPairs{k, 5} == i);
                           end
                           cofrow=find(found);
                           
                           
                           round=0;
                 if ~isempty(cofrow)
                 a{i,1}='cofactor pairs identified';
                 subcofids=[];
                 prodcofids=[];
                       for m = 1:size(cofrow, 1)  
                         round = round + 1;
                         arr{round,1} = subIDs(ismember(subIDs,cofPairs{cofrow(m,1), 6} )); %cofactor reactants
                         arr{round,2} =prodIDs(ismember(prodIDs,cofPairs{cofrow(m,1), 6})); %cofactor products
                                
                        for j = 1:size(arr{round,1},1)
                           subCofids=arr{round,1};
                        end
                        subcofids = [subcofids;subCofids];
                       
                        for j = 1:size(arr{round,2},1)
                           prodCofids=arr{round,2};
                        end
                        prodcofids = [prodcofids;prodCofids];
                        end
                             
                       subcofids=unique(subcofids);
                       prodcofids=unique(prodcofids);
                       
                       cofsubCarbon=sum(abs(full(model.S(subcofids,i))).*metCarb(subcofids,1));
                       cofprodCarbon=sum(abs(full(model.S(prodcofids,i))).*metCarb(prodcofids,1));
                       
                        %check if subCarbon of cofactor reactants equal to prodCarbon of cofactor reactants
                            if cofsubCarbon > cofprodCarbon %more cofactor reactants  than cofactor products, constraining less
                                subCarbon(i) = subCarbon(i) - cofprodCarbon;
                                prodCarbon(i) = prodCarbon(i) - cofprodCarbon;
                           else
                                subCarbon(i) = subCarbon(i) - cofsubCarbon;
                                prodCarbon(i) = prodCarbon(i) - cofsubCarbon;
                            end
                    
                  else
                     a{i,1}='cofactor pairs not identified';
                  end
            

                % subtract 21 carbons form total sub and prd carbon
                % CoA = C21H32N7O16P3S
                if ismember(i,coa.coaRxnIDs)
                     %check if subCarbon of coa reactants equal to prodCarbon of coa reactants
                     %Smaller number is subtracted from carbon balance
                        if coa.prdCarbon(ismember(coa.coaRxnIDs,i)) > coa.subCarbon(ismember(coa.coaRxnIDs,i))
                            subCarbon(i) = subCarbon(i) - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                            prodCarbon(i) = prodCarbon(i) - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                        elseif coa.prdCarbon(ismember(coa.coaRxnIDs,i)) < coa.subCarbon(ismember(coa.coaRxnIDs,i))
                            subCarbon(i) = subCarbon(i) - coa.prdCarbon(ismember(coa.coaRxnIDs,i));
                            prodCarbon(i) = prodCarbon(i) - coa.prdCarbon(ismember(coa.coaRxnIDs,i));
                       
                    else %subcarbon=prodcarbon 
                        subCarbon(i) = subCarbon(i) - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                        prodCarbon(i) = prodCarbon(i) - coa.subCarbon(ismember(coa.coaRxnIDs,i));
                        end
                end
                          

               % subtrac lrgMlcs   
                if ismember(i,lrgMlc.lmRxnIDs)
                    %check if subCarbon of lrgMlc reactants equal to prodCarbon of lrgMlc reactants
                     %Smaller number is subtracted from carbon balance
                        if lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i)) > lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i))
                            subCarbon(i) = subCarbon(i) - lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i));
                            prodCarbon(i) = prodCarbon(i) - lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i));
                        elseif lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i)) < lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i))
                               subCarbon(i) = subCarbon(i) - lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i));
                            prodCarbon(i) = prodCarbon(i) - lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i));
                       
                    else
                    subCarbon(i) = subCarbon(i) - (lrgMlc.subCarbon(ismember(lrgMlc.lmRxnIDs,i)));
                    prodCarbon(i) = prodCarbon(i) - lrgMlc.prdCarbon(ismember(lrgMlc.lmRxnIDs,i));
                        end
                   
                end

              

                if subCarbon(i) == 0 && prodCarbon(i) == 0
                    bericht{i,1} = 'carbon unconstrained';
                    continue
                end

                   % check carbon balancing
                if subCarbon(i) == prodCarbon(i)
                    bericht{i,1} = 'balanced';
                     
                  if ismember(i,CBBRxns) 
                      %if a Calvin cycle reaction under autotrophic conditions
                       constraint = 6*abs(total_carbon_cons/subCarbon(i));
                  else 
                    constraint = abs(total_carbon_cons/subCarbon(i));
                  end
                  
                    if model.qualDir(1,i) == 0 %bidirectional
                        if abs(minmax(i,1)) > constraint
                            minmax_new(i,1) = -constraint;
                            bericht{i,1} = 'constrained-balanced';
                            mmAdj(i,1) = 1;
                        end
                        if abs(minmax(i,2)) > constraint
                           minmax_new(i,2) = constraint;
                           bericht{i,1} = 'constrained-balanced';
                           mmAdj(i,2) = 1;
                        end
                    
                    elseif model.qualDir(1,i) == 1 
                        %unidirectional forward [0 - 100]
                        if abs(minmax(i,2)) > constraint
                            if constraint > abs(minmax(i,1))
                                minmax_new(i,2) = constraint;
                                bericht{i,1} = 'constrained-balanced';
                                mmAdj(i,2) = 1;
                            end
                        end
                    
                
                    elseif model.qualDir(1,i) == -1 
                        %unidirectional backwords [-100 - 0]
                        if abs(minmax(i,1)) > constraint
                            if constraint > abs(minmax(i,2))
                                minmax_new(i,1) = -constraint;
                                bericht{i,1} = 'constrained-balanced';
                                mmAdj(i,1) = 1;
                            end
                        end
                    end
                        
                else % if carbon is not balanced          
                        bericht{i,1} = 'unbalanced';
                        
                        if subCarbon(i) < prodCarbon(i) % if unbalanced then go with the smaller one (constraining less)
                               
                            if ismember(i,CBBRxns)
                                %if a Calvin cycle reaction under autotrophic conditions
                                constraint = 6*abs(total_carbon_cons/subCarbon(i));
                     
                            else 
                            constraint = abs(total_carbon_cons/subCarbon(i));
                            end
                            
                        else
                            if ismember(i,CBBRxns)
                                %if a Calvin cycle reaction under autotrophic conditions
                                constraint = 6*abs(total_carbon_cons/prodCarbon(i));

                            else
                            constraint = abs(total_carbon_cons/prodCarbon(i));
                             end
                         end
                    
                       if model.qualDir(1,i) == 0 %bidirectional
                        if abs(minmax(i,1)) > constraint
                            minmax_new(i,1) = -constraint;
                            bericht{i,1} = 'constrained-unbalanced';
                            mmAdj(i,1) = 1;
                        end
                        if abs(minmax(i,2)) > constraint
                           minmax_new(i,2) = constraint;
                           bericht{i,1} = 'constrained-unbalanced';
                           mmAdj(i,2) = 1;
                        end
                    
                       elseif model.qualDir(1,i) == 1 
                           %unidirectional forward [0 - 100]
                            if abs(minmax(i,2)) > constraint
                            if constraint > abs(minmax(i,1))
                                minmax_new(i,2) = constraint;
                                bericht{i,1} = 'constrained-unbalanced';
                                mmAdj(i,2) = 1;
                            end
                        end
                    
                
                       elseif model.qualDir(1,i) == -1 
                           %unidirectional backwords [-100 - 0]
                            if abs(minmax(i,1)) > constraint
                            if constraint > abs(minmax(i,2))
                                minmax_new(i,1) = -constraint;
                                bericht{i,1} = 'constrained-unbalanced';
                                mmAdj(i,1) = 1;
                            
                            end
                            end 
                       end  
                end 
            end    
        end
        progressbar(i/num_rxns) 
    end
    

    % test if model solves under new conditions
    model.lb = minmax_new(:,1);
    model.ub = minmax_new(:,2);
    try 
    sol = optimizeCbModel(model,'max');
    catch 
        error('The drain reactions bounds that you imposed are possibly unbalanced');
    end 
    sol_cc = sol;
    
    
      % flag to check for infeasibillity
      flag = isnan(sol.f); %infeasible solution 
          
      
      if flag == true  
       if rlxBnds           
                  
          %If true and carbon constraints do not solve due to too      
          %tight constraints, algorithm will increase total carbon       
          %consumed based on relaxInt (relaxation interval) which        
          %is a percentage. 
          
          warning('Model does not solve!')
          rlxVal = rlxVal + rlxInt; 
           %rlxVal=1 (1 = 100%)
           %rlxInt = percentage step totalCarbon will be increased/ relaxed 
           %in case the model does not solve after carbon constraining. Default = 0.1 (i.e. 10%)
           
          total_carbon_cons = intotal_carbon_cons * rlxVal;
          reRun = 1;
       else 
           error('Model does not solve!. Think about relaxing the bounds.')
       end
     else 
        %feasible solution
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
        reRun = 0;
     end 
end 


%% qunatify how much solution space was constraint compared to original minmax
mm_diff = minmax(:,2)-minmax(:,1);
mm_new_diff = minmax_new(:,2)-minmax_new(:,1);
diff = 100*(1-(sum(mm_new_diff)/sum(mm_diff)));
        
%% variable complier
unbRxns = sum(ismember(bericht,'unbalanced'))+sum(ismember(bericht,'unbalanced - CofPair'));
balRxns = sum(ismember(bericht,'balanced'))+sum(ismember(bericht,'balanced - X'))...
    +sum(ismember(bericht,'balanced - drain'))+sum(ismember(bericht,'unconstr - oxPhos'));
setRxns = sum(ismember(bericht,'set_bnd'));
nfRxns = sum(ismember(bericht,'no_flux'));

carbonConst.relaxation = '*** carbon relaxation ***';
if rlxVal > 1
    carbonConst.carbonRelaxed = 'YES';
    carbonConst.rlxVal = rlxVal;
    carbonConst.percentageRelaxed = strcat(num2str(rlxVal*100),'%');
    carbonConst.relaxationIterations = (rlxVal-1)/rlxInt;
    carbonConst.relaxedTotalCarbonConsumed = total_carbon_cons;
else
    carbonConst.carbonRelaxed = 'NO';
end
carbonConst.summary = '*** carbon relaxation ***';
carbonConst.report.summary = bericht;
carbonConst.report.unbalanced = unbRxns;
carbonConst.report.balanced = balRxns;
carbonConst.report.setRxns = setRxns;
carbonConst.report.noFlux = nfRxns;
carbonConst.bndAjusted.rawData = mmAdj;
carbonConst.bndAjusted.numberOf_LB_adj = sum(mmAdj(:,1));
carbonConst.bndAjusted.numberOf_UB_adj = sum(mmAdj(:,2));
carbonConst.CoAevaluation = coa;
carbonConst.electrTransferMlcs = lrgMlc;
carbonConst.carbonCount = carbonCount;
carbonConst.cofactorPairs = cofPairs;
carbonConst.originalMinMax = minmax;
carbonConst.newMinMax = minmax_new;
carbonConst.solSpaceReduction = strcat(num2str(diff),'%');
carbonConst.solutionNew = sol;

%% constrain uptakes based on passed in limit (e.g. lost/closed carbon)

if exist('limEX','var')
    if ~isempty(limEX)
        carbonConst.EXrxnAdjust = '*** EX rxns adjustment ***';
        carbonConst.limEX = limEX;
        
        % get all the max values for the EX rxns not set
        varMaxN = minmax_new(isDrainNew,2);
        % get met indices
        metID_ind = find(model.S(:,isDrainNew));
        % get metIDs from Met indices
        [metID_EX,~] = ind2sub(size(model.S),metID_ind);
        % get rxnIDs where rxn actually has a flux and ignore any fluxes
        % that are set to 0
        tmpRxnID = isDrainNew(find(varMaxN));
        % how many rxns have been set to 0?
        numEXzero = size(isDrainNew,1)-size(tmpRxnID,1);
        % get metIDs from EX rxn IDs so carbon per met can be determined
        tmpMetID = metID_EX(find(varMaxN));
        % calculate new bnds
        limEXnew = limEX./metCarb(tmpMetID);
        % get rid of Inf
        limEXnew(~isfinite(limEXnew)) = 0;
        % find which reactions needf adjustment
        adjEX = tmpRxnID(limEXnew < varMaxN(ismember(isDrainNew,tmpRxnID)));
        %  check if there is any bnds to be adjusted
        if ~isempty(adjEX)
            % make calculations as to how much solution space shrinks
            % compared to original
            diffLim = 100*(1-(sum(limEXnew(ismember(tmpRxnID,adjEX)))/sum(minmax_new(adjEX,2))));
            limEXadj = limEXnew(ismember(tmpRxnID,adjEX));
            minmax_new(adjEX,2) = limEXnew(ismember(tmpRxnID,adjEX));
            
            mm_lim_diff = minmax_new(:,2)-minmax_new(:,1);
            diff_n = 100*(1-(sum(mm_lim_diff)/sum(mm_diff)));
            
            model.lb = minmax_new(:,1);
            model.ub = minmax_new(:,2);
            sol = optimizeCbModel(model,'max');

            if isempty(sol.f)
                warning('model does not solve. think about relaxing the bounds')
            else
                disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
            end
            
            carbonConst.numAdjusted = size(adjEX,1);
            carbonConst.solSpaceReductionEX  = strcat(num2str(diff_n),'%');
            carbonConst.EXspaceReduction = strcat(num2str(diffLim),'%');
        else
            warning('no adjustment of EXrxns was done as already constrained more than limEX')
            carbonConst.warning = 'no adjustment was done';
        end
    end
end

%% check if carbon constraining solves



%% run minmax on new minmax bnds and see if solution space is reduced further
if runFinalMinMax
    model.lb = minmax_new(:,1);
    model.ub = minmax_new(:,2);
    minmax_post_new = runMinMax(model);
    minmax_post_new = fixMinMax(minmax_post_new);
    carbonConst.postNmm = '*** post new minmax ***';
    carbonConst.postNewMinMax = minmax_post_new;
    % qunatify how much solution space was constraint compared to original minmax
    mm_post_diff = minmax_post_new(:,2)-minmax_post_new(:,1);
    diff = 100*(1-(sum(mm_post_diff)/sum(mm_diff)));
    carbonConst.solSpaceReductionPostCC = strcat(num2str(diff),'%');

    % test if new minmax solves
    model.lb = minmax_post_new(:,1);
    model.ub = minmax_post_new(:,2);
    sol = optimizeCbModel(model,'max');
    % find the rxnID that changed after final minmax
    id_lower = find((minmax_new(:,1) - minmax_post_new(:,1)>0));
    id_upper = find((minmax_new(:,2) - minmax_post_new(:,2)>0));
    id_tmp = vertcat(id_lower,id_upper);
    id = sort(unique(id_tmp));
    %initialize vars
    infeas = {};
    infeas_ID = [];
    if isempty(sol.f)
        warning('model does not solve. think about relaxing the bounds')
    end 

    if SearchProb_rxns  % Here we search for reactions that create infeasibillity after constraining with CCFBA.  We do that by constraining reactions one by one.
     if sol.f < sol_cc.f
        warning('Objective Value decreased: searching for reactions that might cause the problem')
        c =0;
    
     for i = 1:length(id)
        model.lb = minmax_new(:,1);
        model.ub = minmax_new(:,2);
        model.lb(i) = minmax_post_new(i,1);
        model.ub(i) = minmax_post_new(i,2);
    
        sol_post = optimizeCbModel(model);
        if sol_post.f < sol_cc.f
         c=c+1;
         infeas(c,1) = model.rxnNames(i);
         infeas_ID(c,1) = i;
      
        end 
        progressbar(i/length(id))
     end
     end 
    carbonConst.Problematic_rxns = infeas;
    carbonConst.Problematic_rxns_ID = infeas_ID;
    end 
     else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
 end

    carbonConst.solutionPostNew = sol;
    
end
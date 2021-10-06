%% DYNAMIC COMPARTMENTAL TJ PERMEABILITY MODEL
% (c) Aapo Tervonen, April 2018
% aapo.tervonen@tut.fi
% Tampere University of Technology
%
% The model simulates diffusion through network of dynamic TJ strands
% using a multi-compartmental approach. An approximate example system
% geometry is presented below.
%      _______________________________________________________
%      |                                                     |
%      |                          40                         |
%      |                                                     |
%      |_____________________________________________________|
%      |30|  31 |  32 |  33 |  34 |  35 |  36 |  37 |  38 |39|
%      |__|_____|_____|_____|_____|_____|_____|_____|_____|__|
%      |  21 |  22 |  23 |  24 |  25 |  26 |  27 |  28 |  29 |
%      |_____|_____|_____|_____|_____|_____|_____|_____|_____|
%      |11|  12 |  13 |  14 |  15 |  16 |  17 |  18 |  19 |20|
%      |__|_____|_____|_____|_____|_____|_____|_____|_____|__|
%      |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |  10 |
%      |_____|_____|_____|_____|_____|_____|_____|_____|_____|
%      |                                                     |
%      |                          1                          |
%      |                                                     |
%      |_____________________________________________________|
%
% The small compartments in the middle are lined by the TJ strands and
% the large compartments above and below the TJ describe the space above
% and below TJs. The model calculates the molecular diffusion throgh the
% strand network by assigning random strand breaks into the strands between
% the compartments. The simulation is run multiple times to average the
% stochaistic strand behavior.    


function TJ_permeability_model(varargin)
% TJ_PERMEABILITY_MODEL Simulate the TJ permeability model
%   The functions takes PARAMETER_NAME, PARAMETER_VALUE pairs as an
%   optional input
%
%   Syntax: TJ_permeability_model(varargin)
%
%   The possible input parameters are (default value in parenthesis):
%       N_SIMULATIONS (256): how many time to run the simulation for
%           averaging
%       N_CORES (the number of physical cores): how many cores to use to
%           run the simulations
%       SIMULATION_TIME (7200): simulated time in seconds
%       N_STRANDS (4): the number of horizontal strands in the system
%       N_COMPARTMENTS 5(50): the horizontal number of the small
%           compartments in the model
%       OPENING_PROBABILITY (0.03): strand break forming probability per µm
%           of strand per second
%       CLOSING_PROBABILITY (0.033): break sealing probability per second
%       JUNCTION_LENGTH (0.180): cell boundary length per area
%       TRICELLULAR_DENSITY (0.014): the density of the tricellular
%           junctions per area
%       SAME_ROW_INTERACTION_LENGTH (100e-9): the length of the strecth of
%           strand between the compartments in the same compartment row
%       OTHER_ROW_INTERACTION_LENGTH (50e-9): the length of the strecth of
%           strand between the compartments in different compartment rows
%       BREAK_WIDTH (20e-9): width of the strand break
%       LARGE_COMPARTMENT_HEIGHT (200e-9): height of the large top and
%           bottom compartments
%       TJ_WIDTH (4e-9): half-width of the TJs between the cells
%       TJ_HEIGHT (6e-9): height of the TJ strand
%       TTJ_RADIUS (5e-9): radius of the tricellular junction pores
%       TTJ_HEIGHT (1e-6): height of the tricellular junctino pores
%       BOTTOM_CONCENTRATION (1): concentration in the bottom compartment
%       MOLECULAR_MASS (457): mass of the permeating PEG molecule
%       INITIAL_VALUES_CASE (2): defines if the simulation begins from the
%       zero initial states (1) or from the equilibrium state (2)
%       EXPORT ([1 1 1 0 0]): matrix that defines what is exported (1 =
%           yes, 0 = no) [RawData CalculatedResults Parameters
%           AllMoleculeQuantities CopyModelFile]
%       EXPORT_TEXT ('no'): text to be added into the exported filenames
%
%   The exported files are names in the following method:
%       STARTINGDATE_STARTINGTIME_FILE_TYPE_SIMULATEDTIME_NSIMULATIONS_
%       NSTRANDS_OPENINGPROBABILITY_CLOSINGPROBABILITY_EXPORTTEXT
%
%   Example:
%       TJ_permeability_model('N_SIMULATIONS', 16)
%       TJ_permeability_model('N_STRANDS', 3, 'JUNCTION_LENGTH', 0.15)
%
%   Other m-files required: none
%   Other files required: initial_values_multipliers.csv (if equilibrium
%       initial conditions are used)
%   Subfunctions: input_check_function, permeability_coefficients_function,
%       interaction_function, compartment_area_function, 
%       unique_interaction_function, rate_constant_function, ode_functions,
%       print_simulation_time_function, calculate_permeability_function,
%       export_function
%
% (c) Aapo Tervonen, April 2018

% Timing
tic

% Print the simulation start time
startingTime = clock;
fprintf('Starting time: %s \r',datestr(startingTime,'HH:MM:SS dd/mmm/yyyy'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign the default parameter values
PARAMETERS = struct(...
    'N_SIMULATIONS', 512,...
    'N_CORES', feature('numCores'),...
    'SIMULATION_TIME', 7200,...
    'N_STRANDS', 4,...
    'N_COMPARTMENTS', 50,...
    'OPENING_PROBABILITY', 0.03,...
    'CLOSING_PROBABILITY', 0.033,...
    'JUNCTION_LENGTH', 0.180,...
    'TRICELLULAR_DENSITY', 0.014,...
    'SAME_ROW_INTERACTION_LENGTH', 100e-9,...
    'OTHER_ROW_INTERACTION_LENGTH', 50e-9,...
    'BREAK_WIDTH', 20e-9,...
    'LARGE_COMPARTMENT_HEIGHT', 200e-9,...
    'TJ_WIDTH', 4e-9,...
    'TJ_HEIGHT', 6e-9,...
    'TTJ_RADIUS', 5e-9,...
    'TTJ_HEIGHT', 1e-6,...
    'BOTTOM_CONCENTRATION', 1,...
    'MOLECULAR_MASS', 547,...
    'INITIAL_VALUES_CASE', 2,...
    'EXPORT', [1 1 1 0 0],...
    'EXPORT_TEXT', 'no');

% Checks the input parameters and assigns them
PARAMETERS = input_check_function(PARAMETERS, varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DERIVED PARAMETERS AND NATURAL CONSTANTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Width of the model system
PARAMETERS.SYSTEM_WIDTH = PARAMETERS.N_COMPARTMENTS*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH*2;

% Area of the bottom compartment
PARAMETERS.LARGE_COMPARTMENT_AREA = PARAMETERS.SYSTEM_WIDTH*PARAMETERS.LARGE_COMPARTMENT_HEIGHT;

% The probability of strand being broken in the beginning of the simulation
% (calculated for strands between different rows and scaled later for the
% same row strands)
strandStates = [0.5 0.5]*[1-PARAMETERS.OPENING_PROBABILITY*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH/1e-6 PARAMETERS.OPENING_PROBABILITY*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH/1e-6; PARAMETERS.CLOSING_PROBABILITY 1-PARAMETERS.CLOSING_PROBABILITY]^1e6;
PARAMETERS.OPEN_PROBABILITY = strandStates(2);

% Avogadros constant
PARAMETERS.N_A = 6.022e23;

% Calculate the permeability coefficients for the breaks and tricellular
% pores
[PARAMETERS.BREAK_PERMEABILITY, PARAMETERS.TRICELLULAR_PERMEABILITY] = permeability_coefficients_function(PARAMETERS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FIND THE INTERACTIONS BETWEEN THE COMPARTMENTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[interactionMatrix, interactionLengths, nFunctions] = interaction_function(PARAMETERS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CREATE A MATRIX FOR THE COMPARTMENT AREAS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compartmentAreas = compartment_area_function(interactionMatrix, nFunctions, PARAMETERS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FIND THE INDICES OF THE INTERACTIONS 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unique linear inds (only including interactions between compartments i
% and j once)
linInteractionInds = unique_interaction_function(interactionMatrix);

% Row and Column indices (for all of the interactions)
[interactionRows, interactionColumns] = find(interactionMatrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SETUP PARALLEL LOOP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Argument for option running parfor with one core
parforArgument = 0;

% Start the parallel pool
if PARAMETERS.N_CORES > 1
    % Create a pool object
    poolObject = gcp('nocreate');
    % Check that if there is a pool running already
    if isempty(poolObject)
        % if no, create a new one
        parpool('local', PARAMETERS.N_CORES);
    else
        % if yes, delete the old one and create a new one
        delete(poolObject);
        parpool('local', PARAMETERS.N_CORES);
    end
    % Set the pargor argument to normal parfor mode
    parforArgument = Inf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   BICELLULAR SIMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the results cell array for the concentrations in the top
% compartments
cData = cell(1,PARAMETERS.N_SIMULATIONS);

% Initialize the results cell array for the final molecule quantities in
% all the compartments
qData = cell(1,PARAMETERS.N_SIMULATIONS);

% Loop for the simulations (run parallel if cores is more than one)
parfor (k = 1:PARAMETERS.N_SIMULATIONS, parforArgument)

    % For timing
    tic
    
    % Short pause so that the random seeds between simulations are
    % different
    pause(k*0.1);
    rng shuffle
    
    % Define the random rate constants for the simulation
    [rateConstants, rateConstantsLinInd] = rate_constant_function(compartmentAreas, interactionLengths, nFunctions, linInteractionInds, interactionRows, interactionColumns, PARAMETERS);
    
    % Initial values for the simulation
    initialValues = initial_values_function(PARAMETERS,nFunctions);
    
    % Simulation is divided into part of 15 minutes to save memory
    simulationParts = ceil(PARAMETERS.SIMULATION_TIME/900);
    
    % Initialization of the result vectors
    qSim = [];
    tSim = [];
    
    % Running the simulation in parts
    for i = 1:simulationParts
        
        % solve the ODE system
        [tTemp, qTemp] = ode23(@(t,Q) ode_functions(t, Q, rateConstants, rateConstantsLinInd), [(i-1)/simulationParts*PARAMETERS.SIMULATION_TIME i/simulationParts*PARAMETERS.SIMULATION_TIME], initialValues );
        
        % Take the last state of the system for the initial state of the
        % next part
        initialValues = qTemp(end,:);
        
        % Save the results
        qSim = [qSim; qTemp(:,end)]; %#ok<*AGROW>
        tSim = [tSim; tTemp(:,end)];
    end
    
    % Create a time series to resample the results
    qTimeSeries = timeseries(qSim,tSim);
    % Resample the results
    qResampled = resample(qTimeSeries,0:PARAMETERS.SIMULATION_TIME);
    
    % Save the resample data that is also converted to concentration
    cData{k} = qResampled.data/PARAMETERS.N_A/PARAMETERS.LARGE_COMPARTMENT_AREA;
    
    % Save the compartments' molecule quantities
    qData{k} = qTemp(end,:)';
    
    % Time of the single simulation run
    tSimulation = toc;
    
    % Print a message the the simulation run has been finished (can be
    % commented out if the number of simulations is high)
    print_simulation_time_function(k, tSimulation, 1)
    
end

% If multiple cores were used, shut down the parfor
if PARAMETERS.N_CORES > 1
    delete(gcp('nocreate'));
end

% Transform the results data saved in the cell into a matrix
cData = cell2mat(cData);

qData = cell2mat(qData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CALCULATE THE RESULTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the bicellular TJ permeability
calculatedPermeabilities = calculate_permeability_function(cData, PARAMETERS);

% Save the tricellular permeability into the same vector
calculatedPermeabilities(2) = PARAMETERS.TRICELLULAR_PERMEABILITY;

% Save the total permeability into the same vector
calculatedPermeabilities(3) = calculatedPermeabilities(1) + calculatedPermeabilities(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FINALIZING THE SIMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export the data
export_function(startingTime, cData, qData, calculatedPermeabilities, PARAMETERS);

% Time of the whole run
tRun = toc;
% Print the duration of the whole simulation run
print_simulation_time_function(1, tRun, 2)

% Beeps for the end :)
beep

end

%% INPUT_CHECK_FUNCTION

function PARAMETERS = input_check_function(PARAMETERS, inputs)
% INPUT_CHECK_FUNCTION Checks the input and assigns the values
%   Takes the input parameter name, parameter value pairs, checks that they
%   are valid and assigns their values to the PARAMETERS struct
%
%   Inputs:
%       PARAMETERS: a struct that holds the default parameter values
%       inputs: the user input parameter name, value pairs
%
%   Output:
%       PARAMETERS: the modified PARAMETERS struct with the new values
%
% (c) Aapo Tervonen, April 2018

% Gets the field names of the PARAMETERS struct
parameterNames = fieldnames(PARAMETERS);

% Gets the number of inputted 
nInputs = length(inputs);

% Checks that the inputs are given in pairs
if mod(nInputs,2) ~= 0
    error('The inputs have to be propertyName/propertyValue pairs');
end

% Go through the input pairs
for inputPair = reshape(inputs,2,[])
    
    % Checks that the property is a valid parameter
    if any(strcmpi(inputPair{1}, parameterNames))
        
        % Checks that the parameter fullfills the requirements for the
        % parameter, if not, gives an error
        if strcmpi(inputPair{1}, 'N_SIMULATIONS')
            assert(isscalar(inputPair{2}) && mod(inputPair{2}, 1.0) == 0 && inputPair{2} > 0, 'Number of simulations must be a positive integer scalar');
        elseif strcmpi(inputPair{1}, 'N_CORES')
            assert(isscalar(inputPair{2}) && mod(inputPair{2}, 1.0) == 0 && inputPair{2} > 0 && inputPair{2} <= feature('numCores'), 'Number of cores must be a positive integer scalar and smaller or equal than the maximum number of physical cores');
        elseif strcmpi(inputPair{1}, 'SIMULATION_TIME')
            assert(isscalar(inputPair{2}) && mod(inputPair{2}, 1.0) == 0 && inputPair{2} > 0, 'Simulation time must be a positive integer scalar');
        elseif strcmpi(inputPair{1}, 'N_STRANDS')
            assert(isscalar(inputPair{2}) && mod(inputPair{2}, 1.0) == 0 && inputPair{2} >= 2, 'Number of strands must be a scalar integer larger or equal to 2');
        elseif strcmpi(inputPair{1}, 'N_COMPARTMENTS')
            assert(isscalar(inputPair{2}) && mod(inputPair{2}, 1.0) == 0 && inputPair{2} > 2, 'Number of compartments must be a scalar integer larger than 2');
        elseif strcmpi(inputPair{1}, 'OPENING_PROBABILITY')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0 && inputPair{2} <= 1, 'Break opening probability must be a scalar in the range [0,1]');
        elseif strcmpi(inputPair{1}, 'CLOSING_PROBABILITY')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0 && inputPair{2} <= 1, 'Break closing probability must be a scalar in the range [0,1]');
        elseif strcmpi(inputPair{1}, 'JUNCTION_LENGTH')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'junction length must be a positive scalar');
        elseif strcmpi(inputPair{1}, 'SECTION_LENGTH')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'Section length must be a positive scalar');
        elseif strcmpi(inputPair{1}, 'BREAK_WIDTH')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0 && inputPair{2} <= 1, 'Break width must be a scalar in the range [0,1]');
        elseif strcmpi(inputPair{1}, 'LARGE_COMPARTMENT_HEIGHT')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'Large compartment height must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'TJ_WIDTH')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'TJ width must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'TJ_HEIGHT')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'TJ heigth must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'TJ_RADIUS')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'tTJ radius must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'TTJ_HEIGHT')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'tTJ height must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'BOTTOM_CONCENTRATION')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'Bottom compartment concentration must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'INITIAL_VALUES_CASE')
            assert(inputPair{2} == 1  || inputPair{2} == 2, 'Initial values case is not 1 or 2');
        elseif strcmpi(inputPair{1}, 'MOLECULAR_MASS')
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'Molecular mass must be a postive scalar');
        elseif strcmpi(inputPair{1}, 'EXPORT')
            isBinary = true;
            for i = inputPair{2}
               isBinary = isBinary && (i == 1 || i == 0); 
            end
            assert(isvector(inputPair{2}) && length(inputPair{2}) == 4 && isBinary, 'Export must be a binary vector of 4 elements');
        elseif strcmpi(inputPair{1}, 'EXPORT_DATA')
            assert(ischar(inputPair{2}), 'Export text must be a string');
        end
        
        % If the parameter passes the checks, assign the value
        PARAMETERS.(upper(inputPair{1})) = inputPair{2};
    else
        
        % The parameter name is not recognized
        error('%s is not a recoqnized parameter name', inputPair{1});
    end
end

end

%% PERMEABILITY_COEFFICIENTS_FUNCTION

function [breakPermeability, tricellularPermeability] = permeability_coefficients_function(PARAMETERS)
% PERMEABILITY_COEFFICIENTS_FUNCTION Calculate break and tTJ permeabilities
%   Calculates the permeability coefficients for the strand breaks and the
%   tTJ pores based on the given molecule size and break and tTJ dimensions
%
%   Input:
%       PARAMETERS: a struct holding the parameters
%   Outputs:
%       breakPermeability: the permeability coefficient of the strand
%           breaks
%       tricellularPermeability: the permeability coefficient of the tTJ
%           pores
%
% (c) Aapo Tervonen, April 2018

% Diffusion coefficient
diffusionCoefficient = 9.9e-9*PARAMETERS.MOLECULAR_MASS^(-0.453);

% Molecular radius (for PEGs)
molecularRadius = 0.29*PARAMETERS.MOLECULAR_MASS^0.454*1e-10;

% Cell boundary length per area (multiplication to convert the unit to m)
cellBoundaryLength = PARAMETERS.JUNCTION_LENGTH/1e-6;

% Relative area of the lateral space
epsilonLateralSpace = 2*PARAMETERS.TJ_WIDTH*cellBoundaryLength;

% Relative area of the tricellular pores (division to change the unit to m)
epsilonTricellularPores = pi*PARAMETERS.TTJ_RADIUS^2*PARAMETERS.TRICELLULAR_DENSITY/1e-12;

% Break permeability coefficient
lambdaBreak = molecularRadius/PARAMETERS.TJ_WIDTH;
breakPermeability = epsilonLateralSpace*diffusionCoefficient/PARAMETERS.TJ_HEIGHT*(1 + 9/16*lambdaBreak*log(lambdaBreak) - 1.19358*lambdaBreak + 0.4285*lambdaBreak^3 - 0.3192*lambdaBreak^4 + 0.08428*lambdaBreak^5);

% Tricellular pore permeability coefficient
lambdaTricellular = molecularRadius/PARAMETERS.TTJ_RADIUS;
tricellularPermeability = epsilonTricellularPores*diffusionCoefficient/PARAMETERS.TTJ_HEIGHT*(1 + 9/8*lambdaTricellular*log(lambdaTricellular) - 1.56034*lambdaTricellular + 0.528155*lambdaTricellular^2 + 1.91521*lambdaTricellular^3 - 2.81903*lambdaTricellular^4 + 0.270788*lambdaTricellular^5 + 1.10115*lambdaTricellular^6 - 0.435933*lambdaTricellular^7);

end

%% INTERACTION_FUNCTION

function [interactionMatrix,interactionLengths, nFunctions] = interaction_function(PARAMETERS)
% INTERACTION_FUNCTION Defines interactions between the compartments
%   Creates a matrix that describes the interactions between the
%   compartments in the model. Each compartment is described both in the
%   rows and columns, and if they are connected in the geometry, value 1 is
%   assigned. If there is no interactions, 0 is assigned.
%
%   Input:
%       PARAMETERS: a struct holding the parameters
%
%   Outputs:
%       interactionMatrix: a matrix describing the interactions
%       interactionLengths: a matrix that has the length of the strand
%           between the compartments for each of the interactions
%       nFunctions: the number of functions for the small compartments in
%           the model
%
% (c) Aapo Tervonen, April 2018

% Calculate the number of functions in the model (the number of strand rows
% times the the number of compartments per row plus one compartment for
% each even row 
nFunctions = (PARAMETERS.N_STRANDS-1)*PARAMETERS.N_COMPARTMENTS + floor((PARAMETERS.N_STRANDS-1)/2);

% Creates a vector that will form the diagonal of the interaction matrix
% (interactions with self)
selfInteractions = ones(nFunctions,1);

% Initialize the vector for the interactions and interaction lengths with
% the next compartments in the same row
nextInteractions = [];
nextInteractionLengths = [];

% Go through the strand rows
for i = 1:PARAMETERS.N_STRANDS-1
    % Even strand rows
    if mod(i,2) == 0
        % Add the interactions for the row and zero in to end for the last
        % compartment
        nextInteractions = [nextInteractions ones(1,PARAMETERS.N_COMPARTMENTS) 0];
        % Add the interaction lenghts similarly to the interactions
        nextInteractionLengths = [nextInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS).*PARAMETERS.SAME_ROW_INTERACTION_LENGTH 0];
    % Odd strand rows
    else
        % Add the interactions for the row and zero in to end for the last
        % compartment
        nextInteractions = [nextInteractions ones(1,PARAMETERS.N_COMPARTMENTS-1) 0];
        % Add the interaction lenghts similarly to the interactions
        nextInteractionLengths = [nextInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS-1).*PARAMETERS.SAME_ROW_INTERACTION_LENGTH 0];
    end
end

% Remove the final extra zeros from the vectors
nextInteractions(end) = [];
nextInteractionLengths(end) = [];

% Create the interactionMatrix including the self and next interactions
% using the diagonal vectors
interactionMatrix = diag(selfInteractions) + diag(nextInteractions,1);
% Create the interactionLenghts matrix including the next interaction
% lengths using the diagonal vector (no length for self interaction)
interactionLengths = diag(nextInteractionLengths,1);

% If there are more than 2 strands
if PARAMETERS.N_STRANDS > 2
    
    % Each compartment (excluding the top row) has two compartments above
    % them, initializing a vector for each and for the lengths
    aboveFirstInteractions = [];
    aboveSecondInteractions = [];
    aboveFirstInteractionLengths = [];
    aboveSecondInteractionLengths = [];
    
    % Go through the compartmentrows excluding the top one
    for i = 1:PARAMETERS.N_STRANDS-2
        % Even row
        if mod(i,2) == 0
            
            % Left above interactions and lengths
            aboveFirstInteractions = [aboveFirstInteractions ones(1,PARAMETERS.N_COMPARTMENTS)];
            aboveFirstInteractionLengths = [aboveFirstInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS).*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH];
            % Rigth above interactions and lengths
            aboveSecondInteractions = [aboveSecondInteractions ones(1,PARAMETERS.N_COMPARTMENTS) 0];
            aboveSecondInteractionLengths = [aboveSecondInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS).*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH 0];
        
            % Odd row 
        else
            % Left above interactions and lengths
            aboveFirstInteractions = [aboveFirstInteractions ones(1,PARAMETERS.N_COMPARTMENTS) 0];
            aboveFirstInteractionLengths = [aboveFirstInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS).*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH 0];
            % Rigth above interactions and lengths
            aboveSecondInteractions = [aboveSecondInteractions ones(1,PARAMETERS.N_COMPARTMENTS)];
            aboveSecondInteractionLengths = [aboveSecondInteractionLengths ones(1,PARAMETERS.N_COMPARTMENTS).*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH];
        end
    end
    
    % If last compartment row is odd, remove the last zeros of the second
    % interactions and lengths
    if mod(PARAMETERS.N_STRANDS-1,2) ~= 0
        aboveSecondInteractions(end) = [];
        aboveSecondInteractionLengths(end) = [];
    end
    
    % Add the above interactions to the interactionMatrix
    interactionMatrix = interactionMatrix + diag(aboveFirstInteractions,PARAMETERS.N_COMPARTMENTS) + diag(aboveSecondInteractions,PARAMETERS.N_COMPARTMENTS+1);
    % Add the above interaction lengths to the interactionLengths matrix
    interactionLengths = interactionLengths + diag(aboveFirstInteractionLengths,PARAMETERS.N_COMPARTMENTS) + diag(aboveSecondInteractionLengths,PARAMETERS.N_COMPARTMENTS+1);
end

% Copy the upper triangle of the matrices to the bottom triangle to make
% them symmetric
interactionMatrix = interactionMatrix + triu(interactionMatrix,1)';
interactionLengths = interactionLengths + triu(interactionLengths,1)';
   
% Interactions of the bottom compartments with other compartments (no
% gain or loss of molecules)
firstRow = zeros(1,nFunctions+2);

% Interactions of the small compartments with the bottom compartment
firstColumn = [ones(PARAMETERS.N_COMPARTMENTS,1);zeros(nFunctions-PARAMETERS.N_COMPARTMENTS,1)];

% Interaction lengths of the small compartment with the large bottom
% compartmnt
firstColumnLengths = [ones(PARAMETERS.N_COMPARTMENTS,1);zeros(nFunctions-PARAMETERS.N_COMPARTMENTS,1)].*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Interactions of the top compartment with other compartments

% Even number of strands
if mod(PARAMETERS.N_STRANDS,2) == 0
    lastRow = [zeros(1,nFunctions-PARAMETERS.N_COMPARTMENTS+1) ones(1,PARAMETERS.N_COMPARTMENTS) 0];
    lastRowLengths = [zeros(1,nFunctions-PARAMETERS.N_COMPARTMENTS+1) ones(1,PARAMETERS.N_COMPARTMENTS) 0].*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
% Odd number of strands
else
    lastRow = [zeros(1,nFunctions-(PARAMETERS.N_COMPARTMENTS+1)+1) ones(1,PARAMETERS.N_COMPARTMENTS+1) 0];
    lastRowLengths = [zeros(1,nFunctions-(PARAMETERS.N_COMPARTMENTS+1)+1) 1 2.*ones(1,PARAMETERS.N_COMPARTMENTS-1) 1 0].*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
end

% Interactions of the small compartments with the top compartment (no
% backflow since the top compartment is considered large)
lastColumn = zeros(nFunctions,1);

% Full interactionMatrix and interactionLengths matrix
interactionMatrix = [firstRow; firstColumn interactionMatrix lastColumn; lastRow];
interactionLengths = [firstRow; firstColumnLengths interactionLengths lastColumn; lastRowLengths];

end

%% COMPARTMENT_AREA_FUNCTION

function compartmentAreas = compartment_area_function(interactionMatrix, nFunctions, PARAMETERS)
% COMPARTMENT_AREA_FUNCTION Calculates the areas for each interaction
%   Creates a matrix that contains the area of each of the donor
%   compartements for the interactions similar
% 
%   Inputs:
%       interactionMatrix - a matrix of the compartment interactions
%       nFunctions: number of small compartment functions
%       PARAMETERS: a struct holding the parameters
%
%   Output:
%       compartmentArea: matrix of the compartment ares for each
%           interaction
%
% (c) Aapo Tervonen, April 2018

% Initialize the compartment areas matrix
compartmentAreas = interactionMatrix;

% Interactions from the bottom compartment (first column)
compartmentAreas(:,1) = compartmentAreas(:,1).*PARAMETERS.LARGE_COMPARTMENT_AREA;

% Interactions from the normal compartments (all the other compartments)
compartmentAreas(2:end,2:end-1) = compartmentAreas(2:end,2:end-1).*PARAMETERS.SAME_ROW_INTERACTION_LENGTH.*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Indices from the small compartments (first and last compartments of
% the even compartment rows)
smallCompartmentInds = 1 + [PARAMETERS.N_COMPARTMENTS+1:2*PARAMETERS.N_COMPARTMENTS+1:nFunctions 2*PARAMETERS.N_COMPARTMENTS+1:2*PARAMETERS.N_COMPARTMENTS+1:nFunctions];
% Interactions from the small compartments
compartmentAreas(:,smallCompartmentInds) = interactionMatrix(:,smallCompartmentInds).*PARAMETERS.SAME_ROW_INTERACTION_LENGTH.*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

end

%% UNIQUE_INTERACTION_FUNCTION

function linearInteractionInds = unique_interaction_function(interactionMatrix)
% UNIQUE_INTERACTION_FUNCTION Finds the indices of the unique interactions
%   Finds the linear indices of the unique interactions (top triangle of
%   the interacionMatrix excluding the diagonal), so that interactions
%   between compartments i and j are only taken into account once
%
%   Input:
%       interactionMatrix: matrix of the compartment interactions
%
%   Output: 
%       linInteractionInds: vector of the linear indices of the
%           interactions
%
% (c) Aapo Tervonen, April 2018

% Gets the rows and columns of all interactions
[interactionInds(1,:), interactionInds(2,:)] = find(interactionMatrix == 1);

% Finds the duplicate interactions (in one half of the matrix)
duplicateInteractions = interactionInds(1,:) <= interactionInds(2,:);

% Removes the duplicates
interactionInds(:,duplicateInteractions) = [];

% Gets the linear indices from the rows and columns
linearInteractionInds = sub2ind(size(interactionMatrix),interactionInds(1,:), interactionInds(2,:));

end

%% RATE_CONSTANT_FUNCTION

function [rateConstants, rateConstantsLinInd] = rate_constant_function(compartmentAreas, interactionLengths, nFunctions, linearInteractionInds, interactionRows, interactionColumns, PARAMETERS)
% RATE_CONSTANT_FUNCTION Defines the rate constants for the interaction
%   Defines the dynamic rate constants for each of the interaction for
%   each time step in the simulation. This is done by changing the break
%   state of the strands with given probabilities on each time step
%
%   Inputs:
%       compartmentAreas: a matrix of the compartment ares for each
%           interactio
%       interactionLengths: a matrix that has the length of the strand
%           between the compartments for each of the interactions
%       nFunctions: number of small compartment functions
%       linInteractionInds: vector of the linear indices of the
%           interactions
%       interactionRows: rows of all the interactions
%       interactionColumnds: columns of all the interactoins
%       PARAMETERS: a struct holding the parameters
%
%   Outputs:
%       rateConstants: a cell of sparse matrices that includes the rate
%			constant matrix for each time step)
%       rateConstantLinInds: a cell of vectors that include the rate
%			constants of the unique inteactions for faster comparison
%			between time points in the ODE function)
%
% (c) Aapo Tervonen, April 2018

% Initializes the rateConstants cell
rateConstants = cell(1,PARAMETERS.SIMULATION_TIME);

% Creates a vector to keep track which interactions have broken strands
% (OPEN_PROBABILITY is scaled for the different strand lengths)
breaks = rand(1,length(linearInteractionInds)) < PARAMETERS.OPEN_PROBABILITY.*interactionLengths(linearInteractionInds)./PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Go through the time
for k = 1:PARAMETERS.SIMULATION_TIME+1
    
    % Create a temporary vector for the rate constant matrix
    rateConstantsTemp = zeros(nFunctions+2,nFunctions+2);
    
    % Assign the rate constants using the linear indices for the broken
    % strands
    rateConstantsTemp(linearInteractionInds) = breaks.*PARAMETERS.BREAK_PERMEABILITY.*PARAMETERS.BREAK_WIDTH;
    
    % Copy the lower half of the matrix to the upper half (without the
    % diagonal)
    rateConstantsTemp = tril(rateConstantsTemp,-1) + tril(rateConstantsTemp,-1)';
    
    % Divide by compartment areas to obtain the rate constants from the
    % permeabilities
    rateConstantsTemp = rateConstantsTemp./compartmentAreas;
    
    % Set the interactions of the bottom compartment to zeros
    rateConstantsTemp(1,:) = 0;
    
    % Set the interactoins of the top row of compartments with the top
    % compartments to zeros (suppress backflow)
    rateConstantsTemp(:,end) = 0;
    
    % Set NaNs to zeros
    rateConstantsTemp(isnan(rateConstantsTemp)) = 0;
    
    % Go through the diagonal values ("self interactions")
    for i = 2:nFunctions+1
        % Sum over all the neighboring compartments and scale with the
        % compartment area to find the rate of loss for each compartment
        rateConstantsTemp(i,i) = -sum(rateConstantsTemp(i,2:end-1).*compartmentAreas(i,2:end-1)./compartmentAreas(i,i)) - rateConstantsTemp(i,1)*PARAMETERS.LARGE_COMPARTMENT_AREA/(PARAMETERS.SAME_ROW_INTERACTION_LENGTH*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH) - rateConstantsTemp(end,i);
    end
    
    % Find the new state of the strands (break and seal depending on the
    % state)
    breaks = breaks + (breaks == 0).*(rand(1,length(linearInteractionInds)) < PARAMETERS.OPENING_PROBABILITY.*interactionLengths(linearInteractionInds)./1e-6) - (breaks == 1).*(rand(1,length(linearInteractionInds)) < PARAMETERS.CLOSING_PROBABILITY);
    
    % Get the linear indices of all the interactions
    linearInds = sub2ind(size(rateConstantsTemp),  interactionRows, interactionColumns);
    
    % Save the rate constants as linear indices into the
    % rateConstantsLinInd cell
    rateConstantsLinInd{k} = rateConstantsTemp(linearInds);
    
    % Save the rate constant matrix as a sparse matrix into the
    % rateConstants cell
    rateConstants{k} = sparse(rateConstantsTemp);
    
end
end

%% INITIAL_CONDITIONS_FUNCTION

function initialValues = initial_values_function(PARAMETERS,nFunctions)
% INITIAL_VALUES_FUNCTION Defines the initial conditions for the
%   simulation either as zero (except for the basal compartment) or as
%   the equilibrium state (loaded from a file). The file to be loaded
%   (initial_value_multipliers.csv) includes the multipliers for the
%   compartment rows, compared to the basal molecular quantity. The file
%   may contain multipliers for systems with different numbers of strands. 
%   In this case, the strand number is first on each row
%   (See the initial_value_multiplier.csv for an example).
%   
%   Input:
%       PARAMETERS: a struct holding the parameters
%       nFunctions: number of small compartment functions
%
%   Output:
%       initialValues: vector of the initial values for the simulation
%
% (c) Aapo Tervonen, April 2018

% If the user wants to simulate from zero initial conditions
if PARAMETERS.INITIAL_VALUES_CASE == 1
    
    % Set the initial conditions of all the compartments to zero
    initialValues = zeros(1,nFunctions+2);
    
    % Set the basal compartment initial state
    initialValues(1) = PARAMETERS.BOTTOM_CONCENTRATION*PARAMETERS.N_A*PARAMETERS.LARGE_COMPARTMENT_AREA;

% If the user wants to simulate from equilibrium state
elseif PARAMETERS.INITIAL_VALUES_CASE == 2
    
    % Load the initial_value_multipliers file
    initialValuesFile = csvread('initial_value_multipliers.csv');
    
    % If the file contains multipliers for only one model system
    if size(initialValuesFile,1) == 1
        
        % Get the multipliers
        multipliers = initialValuesFile;
    
    % If the file contains multipliers for systems with different number of
    % strands
    else
    
        % Get the row index for the simulated strand number
        tempInd = initialValuesFile(:,1) == PARAMETERS.N_STRANDS;
        
        % Check that the strand number is found in the file
        assert(~isempty(tempInd),'No suitable strand number found from the initial_value_multipliers.csv file');
        
        % Get the multipliers
        multipliers = initialValuesFile(tempInd,2:PARAMETERS.N_STRANDS);
        
        % Remove zeros from the vector
        multipliers(multipliers == 0) = [];
        
        % Check that there is correct number of initial value multipliers
        assert(length(multipliers) == PARAMETERS.N_STRANDS-1,'Wrong number of multipliers for the compartment rows in the initial_value_multipliers.csv file');
        
    end
    
    % Calculate the constant value of the basal compartment
    initialValues = PARAMETERS.BOTTOM_CONCENTRATION*PARAMETERS.N_A*PARAMETERS.LARGE_COMPARTMENT_AREA;
    
    % Scales the multipliers based on the area difference of the basal and
    % the small TJ compartments
    multipliers = multipliers*(PARAMETERS.SAME_ROW_INTERACTION_LENGTH^2/PARAMETERS.LARGE_COMPARTMENT_AREA);
    
    % Go through the compartment rows
    for i = 1:PARAMETERS.N_STRANDS-1
       
        % With odd rows
        if mod(i,2) == 1
           
            % Calculate the initial values with the multipliers
            initialValues = [initialValues initialValues(1).*multipliers(i).*ones(1,PARAMETERS.N_COMPARTMENTS)];
        
        % With even rows
        else
            
            % Calculate the initial values with the multipliers (taking
            % into account that the first and last compartment with half of
            % the size compared to others (and thus having half of the
            % molecular quantity)
            initialValues = [initialValues initialValues(1).*multipliers(i)/2 initialValues(1).*multipliers(i).*ones(1,PARAMETERS.N_COMPARTMENTS-1) initialValues(1).*multipliers(i)/2];
            
        end
    end
    
    % Set the initial value of the apical compartment to zero
    initialValues(end+1) = 0;
    
end

end

%% ODEA_FUNCTIONS

function dQ = ode_functions(t,Q, rateConstants, rateConstantsLinInd)
% ODEA_FUNCTIONS The ODE function system solved during the simulation
%	Defines the system of ODEs for each time step by interpolating the
%	strand break behavior between the given time step of 1 seconds
%	
%	Inputs:
%		rateConstants: a cell of sparse matrices that includes the rate
%			constant matrix for each time step)
%       rateConstantLinInds: a cell of vectors that include the rate
%			constants of the unique inteactions for faster comparison
%			between time points in the ODE function)
%
%
% (c) Aapo Tervonen, April 2018

% Gets the floor and ceiling of time for interpolation
tFloor = floor(t)+1;
tCeil = ceil(t)+1;

% Checks if ther are any changes in the system during this time step
% (rateConstantsLinInd is considerably faster for the comparision than the
% full or sparse matrix)
if rateConstantsLinInd{tFloor} == rateConstantsLinInd{tCeil}
    % If not, take the floor rateConstants
    interpolatedReactionConstants = rateConstants{tFloor};
else
    % If yes, do linear interpolation between the two time points
    interpolatedReactionConstants = rateConstants{tFloor} + (rateConstants{tCeil} - rateConstants{tFloor}).*(t - floor(t));
end

% Solve the linear system
dQ = interpolatedReactionConstants*Q;

end

%% PRINT_SIMULATION_TIME_FUNCTION

function print_simulation_time_function(k, tSimulation, option)
% PRINT_SIMULATION_TIME_FUNCTION Prints duration and/or time
% 	Prints either the simulation time and finishing time for one single
%	simulation or the duration for the whole set of simulations
%
%   Inputs:
%		k: simulation index
%		tSimulation: simulation/simulation set duration
%		option: to set what is printed (1 =  single simulation,
%			2 = simulation set)
%
% (c) Aapo Tervonen, April 2018

% Gets the current time
dateVector = clock;

% Formats the time
timeSimulationEnd = datestr(dateVector, 'HH:MM:SS');
dataSimulationEnd = datestr(dateVector, 'dd/mmm/yyyy');

% Different plotting options for different durations

% Over a minute
    if tSimulation >= 60
        % Over an hour
        if tSimulation >= 3600
            % Single simulation
            if option == 1
                runMessage = sprintf('Simulation %.0f finished after %d h, %d min %.2f s at %s, %s', k, floor(tSimulation/3600), floor((tSimulation - 3600*floor(tSimulation/3600))/60), rem(tSimulation - 3600*floor(tSimulation/3600),60), timeSimulationEnd, dataSimulationEnd);
            % Simulation set
            elseif option == 2
                runMessage = sprintf('Calculation time: %d h %d min %.2f s', floor(tSimulation/3600), floor((tSimulation - 3600*floor(tSimulation/3600))/60), rem(tSimulation - 3600*floor(tSimulation/3600),60));
            end
        % Under an hour
        else
            % Single simulation
            if option == 1
                runMessage = sprintf('Simulation %.0f finished after %d min %.2f s at %s, %s', k, floor(tSimulation/60), rem(tSimulation,60), timeSimulationEnd, dataSimulationEnd);
            % Simulation set
            elseif option == 2
                runMessage = sprintf('Calculation time: %d min %.2f s', floor(tSimulation/60), rem(tSimulation,60));
            end
            
        end
    % Under a minute
    else
        % Single simulation
        if option == 1
            runMessage = sprintf('Simulation %.0f finished after %.2f s at %s, %s', k, rem(tSimulation,60), timeSimulationEnd, dataSimulationEnd);
        % Simulation set
        elseif option == 2
            runMessage = sprintf('Calculation time: %.2f s', rem(tSimulation,60));
        end
    end
    
    % Display the message
    disp(runMessage);

end

%% CALCULATE_PERMEABILITY_FUNCTION

function calculatedPermeability = calculate_permeability_function(cData, PARAMETERS)
% CALCULATE_PERMEABILITY_FUNCTION Calculates the bicellular permeability
%   Calculate the bicellular permeability by fitting linear equation to the
%   average of the simulated data
%
%   Inputs:
%       cData: concentration data from the top compartment for all
%           simualations
%       PARAMETERS: a struct holding the parameters
%
%   Output:
%       calculatePermeability: the calculated bicellular permeability
%
% (c) Aapo Tervonen, April 2018

% Create a time vector
t = 0:PARAMETERS.SIMULATION_TIME;

% Convert concentration into molecule quantity
qData = cData.*PARAMETERS.LARGE_COMPARTMENT_AREA;

% Calculate the mean molecule quantity from the multiple simulations
qDataMean = mean(qData,2);

% Find the linear polynomic fit for the data
qFit = polyfit(t', qDataMean,1);

% Get the slope of the fit
permeabilitySlope = qFit(1);

% Calculate the permeability (P = flux/(C_0*A_measurement)
calculatedPermeability = permeabilitySlope/(PARAMETERS.BOTTOM_CONCENTRATION*PARAMETERS.SYSTEM_WIDTH);

end

%% EXPORT_FUNCTION

function export_function(startingTime, cData, qData, calculatedPermeabilities, PARAMETERS)
% EXPORT_FUNCTION Exports the data that the user wants
%   Export the required data, options are raw simulation result data,
%   calculated permeability values, and simulation parameters. Also, the
%   user can copy the simulation files so that the results can be traced
%   back to the code that produced them.
%
%   Inputs:
%       startingTime: starting time of the simulation
%       cData: concentration data from the top compartment for all
%           simualations
%       calculatedPermeabilities: bicellular, tricellular and total
%           permeabilities
%       PARAMETERS: a struct holding the parameters
%
% (c) Aapo Tervonen, April 2018

% Check if something is to be exported and if Results folder exists
if exist('./Results') ~= 7 && sum(PARAMETERS.EXPORT) > 0 %#ok<EXIST>
    % If no folder but something is to be exported, create the Results
    % folder
    mkdir('./', 'Results');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DEFINE THE FILENAMES FOR THE EXPORTED FILENAMES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if something is to be exported
if sum(PARAMETERS.EXPORT) > 0
    
    % Create the date and time strings for the starting time
    dateFilename = datestr(startingTime, 'yyyymmdd');
    timeFilename = datestr(startingTime, 'HHMMSS');
    
    % If no EXPORT_TEXT
    if ~strcmp(PARAMETERS.EXPORT_TEXT,'no')
        
        % RAW file consists of the matrix cData
        filenameRawData = strcat('./Results/RawData/', dateFilename, '_', timeFilename, '_RawData_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', PARAMETERS.EXPORT_TEXT, '.dat');
        
        % Calculated results file consists of the calculated final results
        filenameCalculated = strcat('./Results/CalculatedResults/', dateFilename, '_', timeFilename, '_CalculatedResults_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', PARAMETERS.EXPORT_TEXT, '.dat');
        
        % Parameters file consists of the parameters values in the
        % simulation
        filenameParameters = strcat('./Results/Parameters/', dateFilename, '_', timeFilename, '_Parameters_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', PARAMETERS.EXPORT_TEXT, '.txt');
        
        % AllMoleculeQuantities file consists of the matrix
        % qData
        filenameAllMoleculeQuantities = strcat('./Results/AllMoleculeQuantities/', dateFilename, '_', timeFilename, '_AllMoleculeQuantities_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', PARAMETERS.EXPORT_TEXT, '.dat');
        
        % Name for the folder to input the simulation file copy
        foldernameCode = strcat(dateFilename, '_', timeFilename,'_simulationCode_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', PARAMETERS.EXPORT_TEXT);
    
    % IF EXPORT_TEXT
    else
        
        % RAW file consists of the matrix M_dataOutput
        filenameRawData = strcat('./Results/RawData/', dateFilename, '_', timeFilename, '_RawData_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '.dat');
        
        % Calculated results file consists of the calculated final results
        filenameCalculated = strcat('./Results/CalculatedResults/', dateFilename, '_', timeFilename, '_CalculatedResults_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '.dat');
        
        % Parameters file consists of the parameters values in the
        % simulation
        filenameParameters = strcat('./Results/Parameters/', dateFilename, '_', timeFilename, '_Parameters_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '.txt');
        
        % AllMoleculeQuantities file consists of the matrix
        % qData
        filenameAllMoleculeQuantities = strcat('./Results/AllMoleculeQuantities/', dateFilename, '_', timeFilename, '_AllMoleculeQuantities_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '.dat');
        
        % Name for the folder to input the simulation file copy
        foldernameCode = strcat(dateFilename, '_', timeFilename,'_simulationCode_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_SIMULATIONS), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY));
    end
end

% If RAW data is to be exported
if PARAMETERS.EXPORT(1)
    % Check if the RawData folder already exists
    if exist('./Results/RawData') ~= 7 %#ok<EXIST>
        % If no, make one
        mkdir('./Results/', 'RawData');
    end
    % Write the RawData file
    csvwrite(filenameRawData,[(0:PARAMETERS.SIMULATION_TIME)' cData]);
end

% If calculated results are to be exported
if PARAMETERS.EXPORT(2)
    % Check if the CalculatedResults folder already exist
    if exist('./Results/CalculatedResults') ~= 7 %#ok<EXIST>
        % If no, make one
        mkdir('./Results/', 'CalculatedResults');
    end
    % Write the CalculatedResults file
    fid = fopen(filenameCalculated,'w');
    fprintf(fid,'P_bTJ,P_tTJ,P_total\n');
    fclose(fid);
    dlmwrite(filenameCalculated,calculatedPermeabilities, '-append');
end

% Get the finishing time
finishTime = clock;
% Get the simulation set duration
simulationDuration = toc;

% If parameter file is to be exported
if PARAMETERS.EXPORT(3) == 1
    % Check if the Parameters folder already exist
    if exist('./Results/Parameters') ~= 7 %#ok<EXIST>
        % If no, make one
        mkdir('./Results/', 'Parameters');
    end
    
    % Remove nonnumeric data from struct for writing
    fields2Remove = {'EXPORT', 'EXPORT_TEXT'};
    parametersExport = rmfield(PARAMETERS, fields2Remove);
    
    % Get the field names of the numeric parameter data
    fieldNames = fieldnames(parametersExport);
    
    % Create a matrix of the the PARAMETERS struct
    parameterMatrix = cell2mat(struct2cell(parametersExport));
    
    % Open the parameters file for writing
    fileID = fopen(filenameParameters, 'w');
    
    % Write the starting time, finishing time and calculation duration
    fprintf(fileID, 'Starting time: %s\r\nFinishing time: %s\r\nCalculation time: %d seconds\r\n\r\n', datestr(startingTime),datestr(finishTime),round(simulationDuration));
    
    % Write each of the numeric parameters
    for i = 1:length(parameterMatrix)
        fprintf(fileID,'%s: %s\r\n', fieldNames{i}, num2str(parameterMatrix(i)));
    end
    
    % Write EXPORT matrix
    fprintf(fileID,'%s: %s\r\n', 'EXPORT', num2str(PARAMETERS.EXPORT));

    % Write EXPORT_TEXT
    fprintf(fileID,'%s: %s\r\n', 'EXPORT_TEXT', num2str(PARAMETERS.EXPORT_TEXT));
    
    % Close the file
    fclose(fileID);
end

% If AllMoleculeQuantities data is to be exported
if PARAMETERS.EXPORT(4)
    % Check if the RawData folder already exists
    if exist('./Results/AllMoleculeQuantities') ~= 7 %#ok<EXIST>
        % If no, make one
        mkdir('./Results/', 'AllMoleculeQuantities');
    end
    % Write the RawData file
    csvwrite(filenameAllMoleculeQuantities,qData);
end

% If the simulation code is to be exported (if e.g. changes have been made)
if PARAMETERS.EXPORT(5)
    % Check if the Codes folder already exist
    if exist('./Results/Codes') ~= 7 %#ok<EXIST>
        % If no, make one
        mkdir('./Results/', 'Codes');
    end
    
    % Makes the directory
    mkdir('./Results/Codes', foldernameCode);
    % Creates a string of the pathway
    copyFolder = strcat('./Results/Codes/', foldernameCode);
    % Copies each of the model files to the created folder
    copyfile(strcat(mfilename, '.m'), copyFolder);
end

end
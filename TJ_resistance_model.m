%% DYNAMIC TJ RESISTANCE MODEL
% (c) Aapo Tervonen, April 2018
% aapo.tervonen@tut.fi
% Tampere University of Technology
%
% The model simulates resistance of the dynamic TJ strands using nodal
%   analysis. Example system geometry showing the TJ compartments is
%   presented below (4 compartments and 3 horizontal strands.
%      
%       _______________________________________________
%      |     |           |           |           |     |
%      |  5  |     6     |     7     |     8     |  9  |
%      |_____|___________|___________|___________|_____|
%      |           |           |           |           |
%      |     1     |     2     |     3     |     4     |
%      |___________|___________|___________|___________|
%  
%   The current loops are formed into this system are further presented
%   below (compartments not shown):
%
%                               _____________________________
%                              |                             |
%          ____________________|____________________         |
%         |        |           |           |        |        |
%         R   11   R    12     R    13     R   14   R        |
%         |__R_____|_____R_____|_____R_____|_____R__|        |
%         |     |     |     |     |     |     |     |        +
%         R  4  R  5  R  6  R  7  R  8  R  9  R 10  R   15   V
%         |_____|__R__|_____|__R__|_____|__R__|_____|        -
%            |           |           |           |           |
%            R     1     R     2     R     3     R           |
%            |___________|___________|___________|           |
%                               |                            |
%                               |____________________________|
%
%   The model calculates the resistance of the system by defining the
%   current in the outer loop (15 in the example). The resistors (R) are
%   dynamics, changing based on the random strand breaks in the strands
%   between the compartments. The simulation is run for a long time to
%   average the results.

function TJ_resistance_model(varargin)
% TJ_RESISTANCE_MODEL Simulate the TJ resistance model
%   The functions takes PARAMETER_NAME, PARAMETER_VALUE pairs as an
%   optional input
%   
%   Syntax: TJ_resistance_model(varargin)
%
%   The possible input parameters are (default value in parenthesis):
%       SIMULATION_TIME (1000000): simulated time in seconds
%       N_STRANDS (4): the number of horizontal strands in the system
%       N_COMPARTMENTS (50): the horizontal number of the small
%           compartments
%       OPENING_PROBABILITY (0.03): strand break forming probability per
%           µm of strand per second
%       CLOSING_PROBABILITY (0.033): break sealing probability per second
%       JUNCTION_LENGTH (0.180): cell boundary length per area in µm/µm^2
%       TRICELLULAR_DENSITY (0.0140): the density of the tTJs per area
%           1/µm^2
%       SAME_ROW_INTERACTION_LENGTH (100e-9): the length of the strecth of
%           strand between the compartments in the same compartment row in
%           m
%       OTHER_ROW_INTERACTION_LENGTH (50e-9): the length of the strecth of
%           strand between the compartments in different compartment rows
%           in m
%       BREAK_WIDTH (20e-9): width of the strand break in m
%       TJ_WIDTH (4e-9): half-width of the TJs between the cells in m
%       TJ_HEIGHT (6e-9): height of the TJ strand in m
%       TTJ_RADIUS (5e-9): radius of the tricellular junction pores in m
%       TTJ_HEIGHT (1e-6): height of the tricellular junction pores in m
%       EXTERNAL_VOLTAGE (1): voltage of the voltage source, used for
%           scaling, the actual values does not affect the results
%       FLUID_RESISTIVITY (0.537): resistivity of the physiological fluid
%       STRAND_RESISTANCE (0.3e9): resistance of TJ strand per µm
%       EXPORT([1 1 1 0]) - matrix that defines what is exported (1 = yes,
%           0 = no) [RawData CalculatedResults Parameters CopyModelFile] 
%       EXPORT_TEXT ('no'): text to be added into the exported filenames
%
%   The exported files are names in the following method:
%       STARTINGDATE_STARTINGTIME_SIMULATEDTIME_NSTRANDS_OPENINGPROBABILITY
%       _CLOSINGPROBABILITY_STRANDRESISTANCE_EXPORTTEXT
%
%   Example:
%       TJ_resistance_model('N_SIMULATIONS', 16)
%       TJ_resistance_model('N_STRANDS', 3, 'JUNCTION_LENGTH', 0.15)
%
%   Other m-files required: none
%   Subfunctions: input_check_function, interaction_function, 
%       bicellular_resistance_function, unique_interactions_function,
%       export_function, print_simulation_time_function
%
% (c) Aapo Tervonen, April 2018

% Timing
tic

% Print the simulation start time
startingTime = clock;
fprintf('Starting time: %s \r', datestr(startingTime,'HH:MM:SS dd/mmm/yyyy'));

% Initialize the random number generator
rng shuffle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign the default parameter values
PARAMETERS = struct(...
    'SIMULATION_TIME', 1e6,...
    'N_STRANDS', 4,...
    'N_COMPARTMENTS', 50,...
    'OPENING_PROBABILITY', 0.005,...
    'CLOSING_PROBABILITY', 0.033,...
    'JUNCTION_LENGTH', 0.05,...
    'TRICELLULAR_DENSITY', 0.001,...
    'SAME_ROW_INTERACTION_LENGTH', 100e-9,...
    'OTHER_ROW_INTERACTION_LENGTH', 50e-9,...
    'BREAK_WIDTH', 20e-9,...
    'TJ_WIDTH', 5e-9,...
    'TJ_HEIGHT', 6e-9,...
    'TTJ_RADIUS', 5e-9,...
    'TTJ_HEIGHT', 1e-6,...
    'EXTERNAL_VOLTAGE', 1,...
    'FLUID_RESISTIVITY', 0.537,...
    'STRAND_RESISTANCE', 10.62e9,...
    'EXPORT', [1 1 1 0],...
    'EXPORT_TEXT', 'no');

% Checks the input parameters and assigns them
PARAMETERS = input_check_function(PARAMETERS, varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   DERIVED PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Width of the model system (for the tiled geometry)
PARAMETERS.SYSTEM_WIDTH = PARAMETERS.N_COMPARTMENTS*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH*2;

% The probability of strand being broken in the beginning of the simulation
% (calculated for strands between different rows and scaled later for the
% same row strands)
strandStates = [0.5 0.5]*[1-PARAMETERS.OPENING_PROBABILITY*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH/1e-6 PARAMETERS.OPENING_PROBABILITY*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH/1e-6; PARAMETERS.CLOSING_PROBABILITY 1-PARAMETERS.CLOSING_PROBABILITY]^1e6;
PARAMETERS.OPEN_PROBABILITY = strandStates(2);

% Calculate the resistance of the strand breaks
PARAMETERS.BREAK_RESISTANCE = PARAMETERS.FLUID_RESISTIVITY*PARAMETERS.TJ_HEIGHT/(2*PARAMETERS.TJ_WIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FIND THE INTERACTIONS BETWEEN THE CURRENT LOOPS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[interactionMatrix, interactionLengths, nFunctions] = interaction_function(PARAMETERS);

% find the unique interactions (removes duplicates)
interactionInds = unique_interaction_function(interactionMatrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   BICELLULAR SIMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bicellularResistance = bicellular_resistance_function(interactionLengths, nFunctions, interactionInds, PARAMETERS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TRICELLULAR AND TOTAL RESISTANCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates the tricellular resistance
tricellularResistance = PARAMETERS.FLUID_RESISTIVITY*PARAMETERS.TTJ_HEIGHT/(pi*PARAMETERS.TTJ_RADIUS^2)/(PARAMETERS.TRICELLULAR_DENSITY*1e12)*10000;

% Calculates the total resistance for each time point
TER = 1./(1./bicellularResistance + 1/tricellularResistance);

% Saves the bicellular resistance, tricellular resistance and total TER to
% a vector

% Saves the bicellular resistance, tricellular resistance and total TER to
% a vector
calculatedTERs = mean(bicellularResistance);
calculatedTERs(2) = tricellularResistance;
calculatedTERs(3) = mean(TER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FINALIZING THE SIMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export the data
export_function(startingTime, TER, calculatedTERs, PARAMETERS)

% Time of the simulation run
tSimulation = toc;

% Print the duration of the whole simulation run
print_simulation_time_function(tSimulation);

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
%   (c) Aapo Tervonen, April 2018

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
        % parameter
        if strcmpi(inputPair{1}, 'SIMULATION_TIME')
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
            assert(isscalar(inputPair{2}) && inputPair{2} >= 0, 'Break width must be a positive scalar');
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

%% INTERACTION_FUNCTION

function [interactionMatrix, interactionLengths, nFunctions] = interaction_function(PARAMETERS)
% INTERACTION_FUNCTION Defines interactions between the current loops
%   Finds the indices of the interactions between the current loops in a
%   matrix that describes the interactions between current loops. Each loop
%   is described both in the rows and columns.
%
%   Input:
%       PARAMETERS: a struct holding the parameters
%
%   Outputs:
%       interactionMatrix: a matrix that includes the the interactions
%           between the current loops
%       interactionLengths: a matrix that has the length of the strand
%           between the compartments for each of the interactions 
%           (symmetric, so only upper tringle is populated)
%       nFunctions: the number of functions for the small compartments in
%           the model
%
%   (c) Aapo Tervonen, April 2018

% Calculate the number of functions for the small current loops in the 
% model (excludes the outer current loop) 
nFunctions = (PARAMETERS.N_COMPARTMENTS - 1) + (PARAMETERS.N_STRANDS - 2)*(2*PARAMETERS.N_COMPARTMENTS - 1) + (mod(PARAMETERS.N_STRANDS,2) == 0)*(PARAMETERS.N_COMPARTMENTS - 1) + (mod(PARAMETERS.N_STRANDS,2) == 1)*PARAMETERS.N_COMPARTMENTS;

% Initializes the vector that will form the diagonal of the interaction
% matrix (interactions with self), includes only the strands that are not
% shared by any other current loop (e.g. left side of loops 1, 4 and 11 in 
% the example figure).
selfInteractions = zeros(1, nFunctions);

% The lengths of these strecthes of strand are similarly saved into a
% vector.
selfInteractionsLengths = zeros(1, nFunctions);

% The first current loop always has a resistor that is not shared
selfInteractions(1) = 1;
selfInteractionsLengths(1) = 2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Goes through the current loops rows excluding the first and last
for i = 1:PARAMETERS.N_STRANDS-2
    % Assigns the selfinteractions and lengths for the first loops of each
    % row
    selfInteractions(PARAMETERS.N_COMPARTMENTS + (i-1)*(2*PARAMETERS.N_COMPARTMENTS - 1)) = 1;
    selfInteractionsLengths(PARAMETERS.N_COMPARTMENTS + (i-1)*(2*PARAMETERS.N_COMPARTMENTS - 1)) = PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
end

% Set the selfinteraction for the first loop of the last current loop row
selfInteractions(PARAMETERS.N_COMPARTMENTS + (PARAMETERS.N_STRANDS-2)*(2*PARAMETERS.N_COMPARTMENTS - 1)) = 1;

% Checks if the strand length for the first loop of the last current loop
% row is normal or twice the normal other row interaction length
if mod(PARAMETERS.N_STRANDS,2) == 0
    % if twice the normal (even number of strands)
    selfInteractionsLengths(PARAMETERS.N_COMPARTMENTS + (PARAMETERS.N_STRANDS-2)*(2*PARAMETERS.N_COMPARTMENTS - 1)) = 2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
else
    % otherwise normal (odd number of strands)
    selfInteractionsLengths(PARAMETERS.N_COMPARTMENTS + (PARAMETERS.N_STRANDS-2)*(2*PARAMETERS.N_COMPARTMENTS - 1)) = PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
end

% Initialize the vector for the interactions and interaction lengths with
% the next loops in the same row
sameRowInteractions = [ones(1,PARAMETERS.N_COMPARTMENTS-2) 0];
sameRowInteractionsLengths = [ones(1,PARAMETERS.N_COMPARTMENTS-2) 0].*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Go through the current loop rows (excluding the last row)
for i = 1:PARAMETERS.N_STRANDS-2
    % Add the interactions for the row and zero in to end for the last
    % compartment (since last loop of one row does not share anything with
    % the first loop of the next row)
    sameRowInteractions = [sameRowInteractions ones(1,2*PARAMETERS.N_COMPARTMENTS-2) 0]; %#ok<*AGROW>
    % Add the interaction lenghts similarly to the interactions
    sameRowInteractionsLengths = [sameRowInteractionsLengths ones(1,2*PARAMETERS.N_COMPARTMENTS-2).*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH 0];
end

% Even number of strands
if mod(PARAMETERS.N_STRANDS,2) == 0
    % Add the interactions and interaction length for the row
    sameRowInteractions = [sameRowInteractions ones(1, PARAMETERS.N_COMPARTMENTS-2)];
    sameRowInteractionsLengths = [sameRowInteractionsLengths ones(1, PARAMETERS.N_COMPARTMENTS-2).*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH];
% Odd number of strands
else
    % Add the interactions and interaction length for the row
    sameRowInteractions = [sameRowInteractions ones(1, PARAMETERS.N_COMPARTMENTS-1)];
    sameRowInteractionsLengths = [sameRowInteractionsLengths ones(1, PARAMETERS.N_COMPARTMENTS-1).*2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH];
end

% Finds the indices for the interactions between the current loops of
% different loop rows

% Only two strands
if PARAMETERS.N_STRANDS == 2
    % Get the indices (every loop of the first row interacts with the loop
    % directly above them
    otherRowInteractionsIndices = [1:(PARAMETERS.N_COMPARTMENTS-1); PARAMETERS.N_COMPARTMENTS:(2*(PARAMETERS.N_COMPARTMENTS-1))];
% More than two strands
elseif PARAMETERS.N_STRANDS > 2
    % Loops of the bottom row interact with every other loop in the second
    % row (starting from the second loop)
    otherRowInteractionsIndices = [1:(PARAMETERS.N_COMPARTMENTS - 1) ; PARAMETERS.N_COMPARTMENTS+1:2:(PARAMETERS.N_COMPARTMENTS - 1 + 2*PARAMETERS.N_COMPARTMENTS - 1)];
    % Every other loop of the second to last row interact with the loops of
    % the last row
    otherRowInteractionsIndices = [otherRowInteractionsIndices(1,:)  nFunctions - ((PARAMETERS.N_COMPARTMENTS):2:PARAMETERS.N_COMPARTMENTS+2*PARAMETERS.N_COMPARTMENTS - 1) ; otherRowInteractionsIndices(2,:) nFunctions - (0:(PARAMETERS.N_COMPARTMENTS-1))];
    % If there are more than 3 strands
    if PARAMETERS.N_STRANDS > 3
        % Every other loop in a row interacts with every other loop in the
        % next row
        otherRowInteractionsIndices = [otherRowInteractionsIndices(1,:) PARAMETERS.N_COMPARTMENTS:2:(PARAMETERS.N_COMPARTMENTS + (PARAMETERS.N_STRANDS-3)*(2*PARAMETERS.N_COMPARTMENTS-1))  ; otherRowInteractionsIndices(2,:) (PARAMETERS.N_COMPARTMENTS + 2*PARAMETERS.N_COMPARTMENTS - 1):2:((PARAMETERS.N_COMPARTMENTS + 2*PARAMETERS.N_COMPARTMENTS - 1) + (PARAMETERS.N_STRANDS-3)*(2*PARAMETERS.N_COMPARTMENTS-1))];
    end
end

% Construct the interaction matrix with the diagonal self and next
% interactions
interactionMatrix = (diag(selfInteractions) + diag(sameRowInteractions,1));
% Add the other row interactions into the matrix
interactionMatrix(sub2ind([nFunctions nFunctions], otherRowInteractionsIndices(1,:), otherRowInteractionsIndices(2,:))) = 1;

% Initialize a vector for the last column describing the outer current loop
lastColumn = zeros(nFunctions + 1,1);
% Input the interactions between the outer loop and the last loop of each
% of the rows
lastColumn([(PARAMETERS.N_COMPARTMENTS-1):(2*PARAMETERS.N_COMPARTMENTS - 1):(PARAMETERS.N_COMPARTMENTS-1 + (PARAMETERS.N_STRANDS - 2)*(2*PARAMETERS.N_COMPARTMENTS - 1)) nFunctions]) = 1;

% Add the last column (and last row of zeros) to get the final form of the
% interaction matrix with just the upper triangle populated
interactionMatrix = [[interactionMatrix; zeros(1,nFunctions)] lastColumn];
% Copies the upper triangle to the bottom triangle
interactionMatrix = interactionMatrix + triu(interactionMatrix,1)';

% Construct the interaction length matrix with the diagonal self and next
% interactions
interactionLengths = diag(selfInteractionsLengths) + diag(sameRowInteractionsLengths,1);
% Add the other row interaction lengths into the matrix
interactionLengths(sub2ind([nFunctions nFunctions], otherRowInteractionsIndices(1,:), otherRowInteractionsIndices(2,:))) = PARAMETERS.SAME_ROW_INTERACTION_LENGTH;

% Initialize a vector for the last column in the interactions length matrix
lastColumnLengths = zeros(nFunctions + 1,1);

% If there are even number of strands
if mod(PARAMETERS.N_STRANDS,2) == 0
    % Last loop of the first row
    lastColumnLengths(PARAMETERS.N_COMPARTMENTS-1) = 2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
    % Last loop of the middle rows
    lastColumnLengths(((PARAMETERS.N_COMPARTMENTS-1)+(2*PARAMETERS.N_COMPARTMENTS - 1)):(2*PARAMETERS.N_COMPARTMENTS - 1):(PARAMETERS.N_COMPARTMENTS-1 + (PARAMETERS.N_STRANDS - 2)*(2*PARAMETERS.N_COMPARTMENTS - 1))) = PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
    % Last loop of the last row
    lastColumnLengths(nFunctions) = 2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
% Odd number of strands
else
    % Last loop of the first row
    lastColumnLengths(PARAMETERS.N_COMPARTMENTS-1) = 2*PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
    % Last loop of the middle rows and the last row
    lastColumnLengths([((PARAMETERS.N_COMPARTMENTS-1)+(2*PARAMETERS.N_COMPARTMENTS - 1)):(2*PARAMETERS.N_COMPARTMENTS - 1):(PARAMETERS.N_COMPARTMENTS-1 + (PARAMETERS.N_STRANDS - 2)*(2*PARAMETERS.N_COMPARTMENTS - 1)) nFunctions]) = PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;
end

% Add the last column (and last row of zeros) to get the final form of the
% interaction length matrix with just the upper trianle populated
interactionLengths = [[interactionLengths; zeros(1,nFunctions)] lastColumnLengths];

end

%% UNIQUE INTERACTIONS FUNCTION

function interactionInds = unique_interaction_function(interactionMatrix)
% UNIQUE_INTERACTION_FUNCTION Removes the duplicated interactions in the
% interactionMatrix (between compartments i and j as well as j and i)
%
%   Input:
%       interactionMatrix: a matrix that includes the the interactions
%           between the current loops
%
%   Outputs:
%       interactionInds: a vector that includes the linear indices of all
%       the unique interactions
%
%   (c) Aapo Tervonen, April 2018

% Finds the indices of all the interactions in the system
[interactionInds(1,:), interactionInds(2,:)] = find(interactionMatrix == 1);

% Finds the indices of the duplicate interactions (whose row is larger than
% column)
duplicateInteractions = interactionInds(1,:) > interactionInds(2,:);

% Removes the the duplicate interaction indices
interactionInds(:,duplicateInteractions) = [];

% The linear indices from the row-column indices
interactionInds = sub2ind(size(interactionMatrix),interactionInds(1,:), interactionInds(2,:));

end

%% BICELLULAR RESISTANCE FUNCTION

function bicellularTER = bicellular_resistance_function(interactionLengths, nFunctions, interactionInds, PARAMETERS)
% BICELLULAR_RESISTANCE_FUNCTION Solves the bicellular resistance model
%   Defines the linear system that solves the bicellular resistance model
%   and evolves the model in time by changing the strand state
%
%   Inputs:
%       interactionLengths: a matrix that has the length of the strand
%           shared by the current loop for each of the interactions
%       nFunctions: number of small compartment functions
%       linInteractionInds: vector of the linear indices of the
%           interactions
%       PARAMETERS: a struct holding the parameters
%
%   Outputs:
%       bicellularTER - a vector for the resistance in the outer loop (the
%       bicellular resistance) for each time step
%
%   (c) Aapo Tervonen, April 2018

% Defines the scaling constants for the resistances for later use (scales
% the simulated system width the epithelial-wide values and multiplies by
% 10000 to give the results as Ohm*cm^2 instead of Ohm*m^2)
scalingConstant = PARAMETERS.SYSTEM_WIDTH/(PARAMETERS.JUNCTION_LENGTH*1e6)*10000;

% Defines the voltage vector
V = zeros(nFunctions + 1,1);
% Assigns the outer loop voltage
V(end) = PARAMETERS.EXTERNAL_VOLTAGE;

% Initializes the vector for the results
bicellularTER = zeros(PARAMETERS.SIMULATION_TIME+1,1);

% Creates a vector to keep track which interactions have broken strands
% (OPEN_PROBABILITY is scaled for the different strand lengths)
breaks = rand(1,length(interactionInds)) < PARAMETERS.OPEN_PROBABILITY.*interactionLengths(interactionInds)./PARAMETERS.OTHER_ROW_INTERACTION_LENGTH;

% Gets the row-column indices from the linear indices
[nonLinearRowInds, nonLinearColumnsInds] = ind2sub([nFunctions+1 nFunctions+1], interactionInds);

% Run the simulation through the time steps
for i = 1:PARAMETERS.SIMULATION_TIME+1
    
    % Calculate the resistances of the resistors in the circuit (parallel
    % connection between intact strand and a break, if the strand is
    % broken, and only intact strand if the strand is intact
    resistances = (breaks == 0).*-PARAMETERS.STRAND_RESISTANCE/1e6./interactionLengths(interactionInds)...
        + (breaks == 1).*-1./((PARAMETERS.BREAK_WIDTH./interactionLengths(interactionInds))./(PARAMETERS.BREAK_RESISTANCE/1e6./interactionLengths(interactionInds)) + (1-PARAMETERS.BREAK_WIDTH./interactionLengths(interactionInds))./(PARAMETERS.STRAND_RESISTANCE)/1e6./interactionLengths(interactionInds));
    
    % Create a sparse matrix from the resistors    
    resistances = sparse(nonLinearRowInds, nonLinearColumnsInds, resistances, nFunctions+1, nFunctions+1);
    
    % Mirror the matrix to the bottom half and substract the diagonal twice
    % from the matrix (this basically flips the sign)
    resistances = resistances + triu(resistances,1)' - 2.*diag(diag(resistances));
    
    % Calculate the total diagonal values and adds them into the matrix
    % (since the self interactions were already defined in the last line of
    % code, they are removed from the summation by substracting them from
    % the sum with the last term)
    resistances = resistances + diag(diag(resistances)) - diag(sum(resistances,2));
    
    % Solve the linear system to get the current flowing in each loop
    currents = mldivide(resistances,V);
    
    % Calculate the bicellular resistance from the current in the outer
    % (last) loop with Ohm's law (and scale the results)
    bicellularTER(i) = PARAMETERS.EXTERNAL_VOLTAGE/currents(end)*scalingConstant;
    
    % Progress the state of the breaks for the strands based on the break
    % probabilities
    breaks = breaks + (breaks == 0).*(rand(1,length(interactionInds)) < (PARAMETERS.OPENING_PROBABILITY.*interactionLengths(interactionInds)./1e-6)) - (breaks == 1).*(rand(1,length(interactionInds)) < PARAMETERS.CLOSING_PROBABILITY);
    
end

end

%% EXPORT_FUNCTION

function export_function(startingTime, TER, calculatedTERs, PARAMETERS)
% EXPORT_FUNCTION Exports the data that the user wants
%   Export the required data, options are raw simulation result data,
%   calculated resistance values, and simulation parameters. Also, the
%   user can copy the simulation files so that the results can be traced
%   back to the code that produced them.
%
%   Inputs:
%       startingTime - starting time of the simulation
%       TER - Raw TER data from the simulation
%       calculatedResistances - bicellular, tricellular and total
%           permeabilities
%       PARAMETERS - a struct holding the parameters
%
%   (c) Aapo Tervonen, April 2018

% Check if something is to be exported and if Results folder exists
if exist('./Results_TER') ~= 7 && sum(PARAMETERS.EXPORT) > 0 %#ok<EXIST>
    % If no folder but something is to be exported, create the Results
    % folder
    mkdir('./', 'Results_TER');
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
        
        % RAW file consists of the vector TER
        filenameRawData = strcat('./Results_TER/RawData/', dateFilename, '_', timeFilename, '_RawData_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '_' , PARAMETERS.EXPORT_TEXT, '.dat');
        
        % Calculated results file consists of the calculated final results
        filenameCalculated = strcat('./Results_TER/CalculatedResults/', dateFilename, '_', timeFilename, '_CalculatedResults_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '_', PARAMETERS.EXPORT_TEXT, '.dat');
        
        % Parameters file consists of the parameters values in the simulation
        filenameParameters = strcat('./Results_TER/Parameters/', dateFilename, '_', timeFilename, '_Parameters_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '_', PARAMETERS.EXPORT_TEXT, '.txt');
        
        % Name for the folder to input the simulation file copies
        foldernameCode = strcat(dateFilename, '_', timeFilename,'_simulationCode_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '_', PARAMETERS.EXPORT_TEXT);

    % If EXPORT_TEXT
    else
        % RAW file consists of the vector TER
        filenameRawData = strcat('./Results_TER/RawData/', dateFilename, '_', timeFilename, '_RawData_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '.dat');
        
        % Calculated results file consists of the calculated final results
        filenameCalculated = strcat('./Results_TER/CalculatedResults/', dateFilename, '_', timeFilename, '_CalculatedResults_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '.dat');
        
        % Parameters file consists of the parameters values in the simulation
        filenameParameters = strcat('./Results_TER/Parameters/', dateFilename, '_', timeFilename, '_Parameters_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE), '.txt');
        
        % Name for the folder to input the simulation file copies
        foldernameCode = strcat(dateFilename, '_', timeFilename,'_simulationCode_', num2str(PARAMETERS.SIMULATION_TIME), '_', num2str(PARAMETERS.N_STRANDS),'_', num2str(PARAMETERS.OPENING_PROBABILITY), '_', num2str(PARAMETERS.CLOSING_PROBABILITY), '_', num2str(PARAMETERS.STRAND_RESISTANCE));
        
    end
    
    % If RAW data is to be exported
    if PARAMETERS.EXPORT(1) == 1
        % Check if the RawData folder already exists
        if exist('./Results_TER/RawData') ~= 7 %#ok<EXIST>
            % If no, make one
            mkdir('./Results_TER/', 'RawData');
        end
        % Write the RawData file
        csvwrite(filenameRawData,[(0:PARAMETERS.SIMULATION_TIME)' TER]);
    end
    
    % If calculated results are to be exported
    if PARAMETERS.EXPORT(2) == 1
        % Check if the CalculatedResults folder already exist
        if exist('./Results_TER/CalculatedResults') ~= 7 %#ok<EXIST>
            % If no, make one
            mkdir('./Results_TER/', 'CalculatedResults');
        end
        
        % Write the CalculatedResults file
        fid = fopen(filenameCalculated,'w');
        fprintf(fid,'R_bTJ,R_tTJ,TER\n');
        fclose(fid);
        dlmwrite(filenameCalculated,calculatedTERs, '-append');
    end
    
    % Get the finishing time
    finishTime = clock;
    % Get the simulation set duration
    simulationDuration = toc;
    
    % If parameter file is to be exported
    if PARAMETERS.EXPORT(3) == 1
        
        % Check if the Parameters folder already exist
        if exist('./Results_TER/Parameters') ~= 7 %#ok<EXIST>
            % If no, make one
            mkdir('./Results_TER/', 'Parameters');
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
    
    % If the simulation code is to be exported (if e.g. changes have been
    % made)
    if PARAMETERS.EXPORT(4) == 1
        
        % Check if the Codes folder already exist
        if exist('./Results_TER/Codes') ~= 7 %#ok<EXIST>
            % If no, make one
            mkdir('./Results_TER/', 'Codes');
        end
        
        % Makes the directory
        mkdir('./Results_TER/Codes', foldernameCode);
        % Creates a string of the pathway
        copyFolder = strcat('./Results_TER/Codes/', foldernameCode);
        % Copies each of the model files to the created folder
        copyfile(strcat(mfilename, '.m'), copyFolder);
    end
end

end

%% PRINT_SIMULATION_TIME_FUNCTION

function print_simulation_time_function(tSimulation)
% PRINT_SIMULATION_TIME_FUNCTION Prints duration and/or time
% 	Prints either the simulation time and finishing time for one single
%	simulation or the duration for the whole set of simulations
%
%   Input:
%		Simulation - simulation/simulation set duration
%
%   (c) Aapo Tervonen, April 2018

% Different plotting options for different durations

% Over a minute
if tSimulation >= 60
    % Over an hour
    if tSimulation >= 3600
        runMessage = sprintf('Calculation time: %d h %d min %.2f s', floor(tSimulation/3600), floor((tSimulation - 3600*floor(tSimulation/3600))/60), rem(tSimulation - 3600*floor(tSimulation/3600),60));
    % Under an hour
    else
        runMessage = sprintf('Calculation time: %d min %.2f s', floor(tSimulation/60), rem(tSimulation,60));
    end
% under a minute
else
    runMessage = sprintf('Calculation time: %.2f s', rem(tSimulation,60));
end

% Display the message
disp(runMessage);

end
% INITIALIZE MATLAB
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BANDGAP SEARCH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Desired frequency range (input by user)
desired_freq_min = 0.3; % Example value, change as needed
desired_freq_max = 0.4; % Example value, change as needed

% Tolerance for maximum frequency (slightly above the desired maximum)
max_freq_tolerance = 0.05;

% Search ranges for a, r, and er
a_values = 0.1:0.1:1; % Values from 0 to 1
er_values = 1:1:12; % Values from 1 to 12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEARCH FOR BANDGAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HIGH RESOLUTION GRID 
Nx = 512; % Mapping out accuracy of the space that your unit cell will be in
Ny = Nx;

% PWEM Parameters
P = 3; % Needs to be odd
Q = P;

% RECIPROCAL LATTICE VECTORS
T1 = (2*pi) * [1;0];
T2 = (2*pi) * [0;1];

% KEY POINT OF SYMMETRY
G = [0;0];
X = 0.5*T1;
M = 0.5*T1 + 0.5*T2;

% Generate List
L1 = norm(X-G);
L2 = norm(M-X);
L3 = norm(G-M);

N1 = 20;  % Number of points along G-X (adjust as needed)
N2 = round(N1 * L2 / L1);
N3 = round(N1 * L3 / L1);

BX1 = linspace(G(1), X(1), N1);
BY1 = linspace(G(2), X(2), N1);

BX2 = linspace(X(1), M(1), N2);
BY2 = linspace(X(2), M(2), N2);

BX3 = linspace(M(1), G(1), N3);
BY3 = linspace(M(2), G(2), N3);

% Concatenate the segments
BX = [BX1, BX2, BX3];
BY = [BY1, BY2, BY3];

BETA = [BX ; BY];
BETA(:, [N1+1, N1+N2+1]) = []; % Remove duplicate points

KP = [1, N1, N1+N2-1, N1+N2+N3-2];
KL = {'G', 'X', 'M', 'G'};

% Initialize variables to store the closest bandgap
best_freq_match = inf;
best_params = [];

% Loop over all combinations of a, r, and er
for a = a_values
    r_values = 0.1*a:0.1*a:0.5*a; % Values from 0.1*a to 0.5*a
    for r = r_values
        for er = er_values
            % Build unit cell
            dx = a / Nx;
            dy = a / Ny;
            xa = (1:Nx) * dx; 
            xa = xa - mean(xa);
            ya = (1:Nx) * dx; 
            ya = ya - mean(ya);
            [Y, X] = meshgrid(ya, xa);

            UR = ones(Nx, Ny);
            ER = (X.^2 + Y.^2) < r^2; % Currently a bool because all the stuff outside of the circle is 0 and inside is 1
            ER = 1.0 + (er - 1.0) * ER; % Now we get an array of dielectric cylinders

            % Compute convolution matrices
            URC = convmat(UR, P, Q); % Convolution of the uniform relative permeability UR over the given number of spatial harmonics P and Q
            ERC = convmat(ER, P, Q); % Constructing a matrix that, when multiplied with the spatial harmonic components of the wave vector, simulates the effect of the uniform relative permeability on the electromagnetic fields within the unit cell

            % Initialize band data
            NBETA = length(BETA(1,:));
            K0 = zeros(P*Q, NBETA); % Store the sorted eigenvalues (band energies) for each beta value

            % Main loop -- iterate over Bloch wave vector
            for nbeta = 1:NBETA
                % Get next beta
                bx = BETA(1, nbeta);
                by = BETA(2, nbeta);

                % Form K matrices
                KX = bx - 2*pi*linspace(-floor(P/2), floor(P/2), P)/a;
                KY = by - 2*pi*linspace(-floor(Q/2), floor(Q/2), Q)/a;
                [KY, KX] = meshgrid(KY, KX);
                KX = diag(KX(:)); % Creates a diagonal matrix where the elements of the vector are placed along the diagonal
                KY = diag(KY(:));

                % Compute A matrix
                A = (KX / URC) * KX + (KY / URC) * KY; % Represents the Hamiltonian matrix in the eigenvalue problem, where KX and KY contribute to the wave vector components adjusted by the dielectric properties (URC)
                B = ERC;
                D = eig(full(A), full(B));
                D = real(sqrt(D)); % To convert them to frequencies (since eigenvalues correspond to squared frequencies in the PWEM)
                K0(:,nbeta) = sort(D);
            end

            % Calculate normalized frequency
            WN = a*K0/(2*pi); % Aligns the frequencies with the physical dimensions of the photonic crystal

            % Find bandgaps
            for i = 1:(size(WN,1)-1)
                % Find the minimum and maximum frequencies of successive bands
                min_upper_band = min(WN(i+1,:));
                max_lower_band = max(WN(i,:));
                
                % Check if the bandgap is valid
                if min_upper_band > max_lower_band
                    % Calculate frequency match score
                    freq_match_score = abs(min_upper_band - desired_freq_max) + abs(max_lower_band - desired_freq_min);
                    
                    % Debugging Information
                    fprintf('Debugging Information:\n');
                    fprintf('a = %.4f, r = %.4f, er = %.4f\n', a, r, er);
                    disp('Computed Frequencies:');
                    disp([max_lower_band, min_upper_band]);
                    disp('Frequency Match Score:');
                    disp(freq_match_score);

                    % Update best parameters based on frequency match
                    if freq_match_score < best_freq_match
                        best_freq_match = freq_match_score;
                        best_params = [a, r, er, max_lower_band, min_upper_band, min_upper_band - max_lower_band];
                    end
                end
            end
        end
    end
end

% Display the best frequency match
if ~isempty(best_params)
    disp('Best Frequency Match:');
    disp('      a        r        er       min_freq   max_freq   bandgap');
    disp(best_params);
else
    disp('No suitable frequency match found.');
end

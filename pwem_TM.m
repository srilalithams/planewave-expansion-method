%INITIALIZE MATLAB 
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL
a = 1;
r = 0.35 * a;
er = 9.0;

% HIGH RESOLUTION GRID 
Nx = 512; %mapping out accuracy of the space that your unit cell will be in
Ny = Nx;

% PWEM Parameters % number of spatial harmonics that we're using
P = 3; % needs to be odd
Q = P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESHGRID
dx = a / Nx;
dy = a / Ny;
xa = [1:Nx] * dx; 
xa = xa - mean(xa);
ya = [1:Nx] * dx; 
ya = ya - mean(ya);
[Y, X] = meshgrid(ya, xa);

% BUILD UNIT CELL
UR = ones(Nx, Ny);
ER = (X.^2 + Y.^2) < r^2; % currently a bool because all the stuff outside of the circle is 0 and inside is 1
ER = 1.0 + (er - 1.0) * ER; % now we get an array of dielectric cylinders

% VISUALIZE UNIT CELL
%figure;
imagesc(xa, ya, ER'); 
colorbar;
axis equal tight;

% COMPUTE CONVOLUTION MATRICES
URC = convmat(UR, P, Q); %convolution of the uniform relative permeability UR over the given number of spatial harmonics P and Q
ERC = convmat(ER, P, Q); %constructing a matrix that, when multiplied with the spatial harmonic components of the wave vector, simulates the effect of the uniform relative permeability on the electromagnetic fields within the unit cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE LIST OF BLOCH WAVE VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RECIPROCAL LATTICE VECTORS
T1 = (2*pi/a) * [1;0];
T2 = (2*pi/a) * [0;1];

% KEY POINT OF SYMMETRY
G = [0;0];
X = 0.5*T1;
M = 0.5*T1 + 0.5*T2;

% GENERATE LIST
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

KP = [ 1, N1, N1+N2-1, N1+N2+N3-2];
KL = {'G' 'X' 'M' 'G'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PWEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE SPATIAL HARMONIC INDICES
p = -floor(P/2):floor(P/2);
q = -floor(Q/2):floor(Q/2);

% INITIALIZE BAND DATA
NBETA = length(BETA(1,:));
K0 = zeros(P*Q, NBETA); %store the sorted eigenvalues (band energies) for each beta value

figure;

%
% MAIN LOOP -- ITERATE OVER BLOCH WAVE VECTOR
%
for nbeta = 1:NBETA

    % GET NEXT BETA
    bx = BETA(1, nbeta);
    by = BETA(2, nbeta);
    %x and y components of the current Bloch wave vector

    % FORM K MATRICES
    KX = bx - 2*pi*p/a; 
    KY = by - 2*pi*q/a;
    %2*pi/a scales the reciprocal lattice vectors to match the units of bx and by
    %p/a and q/a adjust the reciprocal lattice vectors by the spatial harmonics
    %KX and KY are matrices that encode the adjusted wave vector components in reciprocal space
    [KY, KX] = meshgrid(KY, KX);
    KX = diag(KX(:)); %creates a diagonal matrix where the elements of the vector are placed along the diagonal
    KY = diag(KY(:));

    % COMPUTE A MATRIX
    A = (KX / URC) * KX + (KY / URC) * KY; %represents the Hamiltonian matrix in the eigenvalue problem, where KX and KY contribute to the wave vector components adjusted by the dielectric properties (URC)
    B = ERC;
    D = eig(full(A),full(B));
    D = real(sqrt(D)); %to convert them to frequencies (since eigenvalues correspond to squared frequencies in the PWEM)
    K0(:,nbeta) = sort(D);

    %SHOW BANDS
    % subplot(1,3,2:3);
    % plot([1:nbeta],K0(:, 1:nbeta), '.b');
    % xlim([1 NBETA]);
    % ylim([0 5]);
    % drawnow;
end

%CALCULATE NORMALIZED FREQUENCY
WN = a*K0/(2*pi); %aligns the frequencies with the physical dimensions of the photonic crystal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CLEAR FIGURE WINDOW
clf;

%PLOT BANDS
plot([1:NBETA], WN, '-b');
hold on;

%HIGHLIGHT BAND GAPS
for i = 1:(size(WN,1)-1)
    % Find the minimum and maximum frequencies of successive bands
    min_upper_band = min(WN(i+1,:));
    max_lower_band = max(WN(i,:));
    
    % If there is a gap, highlight it
    if min_upper_band > max_lower_band
        fill([1 NBETA NBETA 1], [max_lower_band max_lower_band min_upper_band min_upper_band], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end

%SET AXIS LIMITS
xlim([1 NBETA]);
ylim([0 1]);

%SET TICK MARKS
set(gca, 'XTick', KP, 'XTickLabel',KL);
for n = 1 : length(KP)
    line(KP(n)*[1 1],[0 1],'Color','k','LineStyle',':');
end

%LABEL AXES
xlabel('Block Wave Vector, $\vec{\beta}$', 'Interpreter', 'LaTex');
ylabel('Normalized Frequency, $\omega a/2pi c_0$', 'Interpreter', 'LaTex');

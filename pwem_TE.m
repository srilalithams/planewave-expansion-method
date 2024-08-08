% Clear workspace and command window
clear all; 
clc;

% Define constants
a = 1e-6;       % Period of the photonic crystal (meters)
r = 0.4 * a;    % Radius of the elements (meters)
eps1 = 1;       % Relative permittivity of the background material
eps2 = 12;       % Relative permittivity of the elements
precis = 5;     % Number of k-vector points between high symmetry points
nG = 4;         % Number of plane waves (half-grid size)
precisStruct = 30; % Number of nodes in the discretization mesh

% Initialize unit cell
nx = 1;
for countX = -a/2 : a/precisStruct : a/2
    ny = 1;
    for countY = -a/2 : a/precisStruct : a/2
        if sqrt(countX^2 + countY^2) < r
            struct(nx, ny) = 1/eps2;
        else
            struct(nx, ny) = 1/eps1;
        end
        xSet(nx) = countX;
        ySet(ny) = countY;
        ny = ny + 1;
    end
    nx = nx + 1;
end

% Compute mesh cell area
dS = (a / precisStruct)^2;

% Create 2D arrays of mesh coordinates
[xMesh, yMesh] = meshgrid(xSet(1:end-1), ySet(1:end-1));

% Transform values of inverse dielectric function
structMesh = struct(1:end-1, 1:end-1) * dS / (max(xSet) - min(xSet))^2;

% Define k-path in the Brillouin zone
kx = [0 : pi/a/precis : pi/a, pi/a * ones(1, precis), ...
      pi/a : -pi/a/precis : 0];
ky = [zeros(1, precis+1), pi/a/precis : pi/a/precis : pi/a, ...
      pi/a : -pi/a/precis : 0];

% Define reciprocal lattice vectors
numG = 1;
for Gx = -nG*2*pi/a : 2*pi/a : nG*2*pi/a
    for Gy = -nG*2*pi/a : 2*pi/a : nG*2*pi/a
        G(numG, 1) = Gx;
        G(numG, 2) = Gy;
        numG = numG + 1;
    end
end

% Compute Fourier expansion coefficients
for countG = 1 : numG - 1
    for countG1 = 1 : numG - 1
        CN2D_N(countG, countG1) = sum(sum(structMesh .* ...
            exp(1i * ((G(countG, 1) - G(countG1, 1)) * xMesh + ...
                      (G(countG, 2) - G(countG1, 2)) * yMesh))));
    end
end

% Compute matrix differential operator for TE mode
for countG = 1 : numG - 1
    for countG1 = 1 : numG - 1
        for countK = 1 : length(kx)
            M(countK, countG, countG1) = ...
                CN2D_N(countG, countG1) * ((kx(countK) + G(countG, 1)) * ...
                (kx(countK) + G(countG1, 1)) + ...
                (ky(countK) + G(countG, 2)) * (ky(countK) + G(countG1, 2)));
        end
    end
end

% Compute eigen-states for each wave vector
for countK = 1 : length(kx)
    MM(:,:) = M(countK, :, :);
    [D, V] = eig(MM);
    dispe(:, countK) = sqrt(diag(V)) * a / (2 * pi);
end

% Plot the band structure
figure;
hold on;
for u = 1 : 8
    plot(abs(dispe(u, :)), 'r', 'LineWidth', 2);
    if u < 8 && min(dispe(u + 1, :)) > max(dispe(u, :))
        rectangle('Position', [1, max(dispe(u, :)), ...
            length(kx) - 1, min(dispe(u + 1, :)) - max(dispe(u, :))], ...
            'FaceColor', 'b', 'EdgeColor', 'b');
    end
end
% Labeling the axes
set(gca, 'xtick', [1 precis+1 2*precis+1 3*precis+1]);
set(gca, 'xticklabel', {'Γ', 'X', 'M', 'Γ'});
ylabel('Frequency $\omega a/2\pi c$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('Wavevector', 'FontSize', 16);
xlim([1, length(kx)]);
grid on;
hold off;

% Find the first bandgap frequencies
first_bandgap_found = false;
for u = 1 : 7
    if min(dispe(u + 1, :)) > max(dispe(u, :))
        first_bandgap_min = max(dispe(u, :));
        first_bandgap_max = min(dispe(u + 1, :));
        first_bandgap_found = true;
        break;
    end
end

% Display the first bandgap frequencies
if first_bandgap_found
    fprintf('First bandgap frequencies:\n');
    fprintf('Lower frequency: %.4f\n', first_bandgap_min);
    fprintf('Upper frequency: %.4f\n', first_bandgap_max);
else
    disp('No bandgap found.');
end

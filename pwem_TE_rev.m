% Clear workspace and command window
clear all; 
clc;

% Define constants
eps1 = 1;       % Relative permittivity of the background material
eps2 = 12;      % Relative permittivity of the elements
precis = 5;     % Number of k-vector points between high symmetry points
nG = 4;         % Number of plane waves (half-grid size)
precisStruct = 30; % Number of nodes in the discretization mesh
tolerance = 1e-3; % Tolerance for detecting a bandgap

% Define target bandgap range (example values)
target_bandgap_min = 0.5;
target_bandgap_max = 0.7;

% Define range for 'a' and r_ratio
a_range = 0.1e-6 : 0.1e-6 : 1e-6;
r_ratio_range = 0.1 : 0.1 : 0.5; % r is defined as a fraction of a

% Initialize variables to store optimal 'a' and 'r' values
optimal_a = 0;
optimal_r = 0;
optimal_bandgap_min = 0;
optimal_bandgap_max = 0;
found = false;

% Loop over a and r_ratio values
for a = a_range
    for r_ratio = r_ratio_range
        r = r_ratio * a;

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

        % Find the first bandgap frequencies with tolerance
        first_bandgap_found = false;
        for u = 1 : 7
            if min(dispe(u + 1, :)) > max(dispe(u, :) + tolerance)
                first_bandgap_min = max(dispe(u, :));
                first_bandgap_max = min(dispe(u + 1, :));
                first_bandgap_found = true;
                break;
            end
        end

        % Print the bandgap found for each iteration
        if first_bandgap_found
            fprintf('a = %.4e m, r = %.4e m, Bandgap: %.4f - %.4f\n', ...
                a, r, first_bandgap_min, first_bandgap_max);
        else
            fprintf('a = %.4e m, r = %.4e m, No bandgap found\n', a, r);
        end

        % Check if the first bandgap falls within the target range
        if first_bandgap_found && ...
           first_bandgap_min >= target_bandgap_min && ...
           first_bandgap_max <= target_bandgap_max
            optimal_a = a;
            optimal_r = r;
            optimal_bandgap_min = first_bandgap_min;
            optimal_bandgap_max = first_bandgap_max;
            found = true;
            break;
        end
    end
    if found
        break;
    end
end

% Display the optimal 'a' and 'r' values and the corresponding bandgap
if found
    fprintf('Optimal values found:\n');
    fprintf('a = %.4e m\n', optimal_a);
    fprintf('r = %.4e m\n', optimal_r);
    fprintf('First bandgap frequencies:\n');
    fprintf('Lower frequency: %.4f\n', optimal_bandgap_min);
    fprintf('Upper frequency: %.4f\n', optimal_bandgap_max);
else
    disp('No suitable values found within the specified range.');
end

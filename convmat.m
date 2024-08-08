function C = convmat(A, P, Q)
    % Get the size of the input matrix A
    [Nx, Ny] = size(A);

    % Perform 2D Fourier transform of the input matrix A
    A_fft = fftshift(fft2(A)) / (Nx * Ny);

    % Determine the center index
    cx = floor(Nx / 2) + 1;
    cy = floor(Ny / 2) + 1;

    % Initialize the convolution matrix
    C = zeros(P * Q, P * Q);

    % Loop over the spatial harmonics
    for p = 1:P
        for q = 1:Q
            for m = 1:P
                for n = 1:Q
                    % Compute the indices for the Fourier components
                    i = cx + (p - ceil(P / 2)) - (m - ceil(P / 2));
                    j = cy + (q - ceil(Q / 2)) - (n - ceil(Q / 2));
                    
                    % Handle boundary conditions
                    if (i > 0 && i <= Nx && j > 0 && j <= Ny)
                        % Compute the indices for the convolution matrix
                        row = (p - 1) * Q + q;
                        col = (m - 1) * Q + n;
                        
                        % Assign the Fourier component to the convolution matrix
                        C(row, col) = A_fft(i, j);
                    end
                end
            end
        end
    end
end

function A = calib2mat(C, k)
    % calibration matrix = calib2mat(calibration region, kernel size)
    %
    % Function to create the calibration matrix for ESPIRiT.
    %
    % Inputs:
    %	C - Calibration region [kx, ky, coils]
    %	k - Kernel size [k]
    %
    % Outputs:
    %	A - Calibration matrix.

    kernel = [k, k];
    tmp = im2row(C, kernel); [tsx,tsy,tsz] = size(tmp);
    A = reshape(tmp, tsx, tsy * tsz);

end

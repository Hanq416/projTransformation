# find your camera's parameters
# [1] sensor's dimension x*y in mm
# [2] focal length in mm
# Lens field of vision, 180 degree

# Calibrate camera's vignetting cureve in your lab
# use polynomial fit, degree-3 or degree-4
# Change it in the Functions's Library

#####
function [vcf] = canonVC(angle) % vc correction for Canon t2i + sigma f2.8
# change this function:

vcf = 1./(-4.3909e-09.*angle.^(4) - 3.9024e-07.*angle.^(3) +...
    3.3680e-05.*angle.^(2)-0.0018.*angle + 1.0018);
end

#####

%% Function to covert COE to R,V
function [r,v]=coe2rv(coe,MU)
rad_multiplier=pi/180;
a    = coe(1);
e    = coe(2);
incl = coe(3)*rad_multiplier;
RA   = coe(4)*rad_multiplier;
w    = coe(5)*rad_multiplier;
TA   = coe(6)*rad_multiplier;
h = sqrt(MU*a*(1-e^2));
rp = (h^2/MU) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (MU/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);
R3_W = [ cos(RA)  sin(RA)  0
        -sin(RA)  cos(RA)  0
            0        0     1];
R1_i = [1       0          0
        0   cos(incl)  sin(incl)
        0  -sin(incl)  cos(incl)];
R3_w = [ cos(w)  sin(w)  0 
        -sin(w)  cos(w)  0
           0       0     1];
ROT_MATRIX = (R3_w*R1_i*R3_W)';
r = ROT_MATRIX*rp;
v = ROT_MATRIX*vp;
end
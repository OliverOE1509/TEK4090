A_r = 1;
L_r = 1;
E_r = 1;
rao_r = 6;
M1 = 10;
M2 = 10;
alpha = 0.1;

M_lump = [10 0; 0 10];
M_rod = (A_r * L_r * rao_r) / 6 * [2 1; 1 2];
K_rod = (A_r * E_r) / L_r * [1 -1; -1 1];
C_rod = alpha * K_rod;


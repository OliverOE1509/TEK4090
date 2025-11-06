A = [0 0 1 0;
     0 0 0 1;
     -13/143 13/143 -1.3/143 1.3/143;
     13/143 -13/143 1.3/143 -1.3/143];
B = [0; 0; 12/143; -1/143];
C = [0 1 0 0];
D = 0;

sys = ss(A, B, C, D);
[num, den] = tfdata(sys, 'v');
num = round(num, 6);  den = round(den, 6);
disp('Clean transfer function:');
sys_tf_clean = tf(num, den);
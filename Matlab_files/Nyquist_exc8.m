num = [1,0,1];
den = [1,2,2,1];
sys = tf(num, den);
nyquist(sys)
nC = [1];
dC = [1,0];

nP= [1,2];
dP = [1,0,2];

numProd = conv(nP, nC);
denProd = conv(dC, dP);

ml = max(length(denProd), length(numProd));
denProd_padded = [zeros(1, ml - length(denProd)), denProd];
numProd_padded = [zeros(1, ml - length(numProd)), numProd];

char_poly = denProd_padded + numProd_padded;

r = roots(char_poly);
disp(r)
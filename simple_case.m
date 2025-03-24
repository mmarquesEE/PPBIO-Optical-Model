% Inclua o pacote de simbolos
syms d1 d2 d3 d4 f_c f_f n_p

% Defina as matrizes simbólicas
MLc1 = [1, d1; -1/f_c, 1 - d1/f_c];
MLf = [1, d2; -1/f_f, 1 - d2/f_f];
MS1 = [1, d3; 0, 1/n_p];
Msa = [1, d4; 0, 1];

% Calcule o produto das matrizes
M1 = Msa * MS1 * MLf * MLc1;

% Simplifique o resultado, se necessário
M1 = simplify(M1);

A = M1(1,1); B = M1(1,2); C = M1(2,1); D = M1(2,2);

latex(A)
latex(B)
latex(C)
latex(D)
% Define time vector
t = linspace(0, 10, 1000);

% Constant concentration
c = 1;

% Parameter Set 1 (p1)
% Region A
kon1_A = 1; koff1_A = 2; smax1_A = 3;
alpha1_A = kon1_A * smax1_A;

% Region B
kon1_B = 4; koff1_B = 5; smax1_B = 6;
alpha1_B = kon1_B * smax1_B;

% Compute s for p1
sA_p1 = (alpha1_A * c) / (kon1_A * c + koff1_A) .* (1 - exp(-(kon1_A * c + koff1_A) * t));
sB_p1 = (alpha1_B * c) / (kon1_B * c + koff1_B) .* (1 - exp(-(kon1_B * c + koff1_B) * t));
sobs_p1 = sA_p1 + sB_p1;

% Parameter Set 2 (p2: swapped regions)
kon2_A = kon1_B; koff2_A = koff1_B; smax2_A = smax1_B;
alpha2_A = kon2_A * smax2_A;

kon2_B = kon1_A; koff2_B = koff1_A; smax2_B = smax1_A;
alpha2_B = kon2_B * smax2_B;

% Compute s for p2
sA_p2 = (alpha2_A * c) / (kon2_A * c + koff2_A) .* (1 - exp(-(kon2_A * c + koff2_A) * t));
sB_p2 = (alpha2_B * c) / (kon2_B * c + koff2_B) .* (1 - exp(-(kon2_B * c + koff2_B) * t));
sobs_p2 = sA_p2 + sB_p2;

% Check equivalence between p1 and p2
difference_p1_p2 = max(abs(sobs_p1 - sobs_p2));
disp(['Max difference (p1 vs p2): ', num2str(difference_p1_p2)]);

% Parameter Set 3 (p3: non-equivalent)
kon3_A = kon1_A; koff3_A = koff1_A + 0.1; smax3_A = smax1_A;
alpha3_A = kon3_A * smax3_A;

kon3_B = kon1_B; koff3_B = koff1_B; smax3_B = smax1_B;
alpha3_B = kon3_B * smax3_B;

% Compute s for p3
sA_p3 = (alpha3_A * c) / (kon3_A * c + koff3_A) .* (1 - exp(-(kon3_A * c + koff3_A) * t));
sB_p3 = (alpha3_B * c) / (kon3_B * c + koff3_B) .* (1 - exp(-(kon3_B * c + koff3_B) * t));
sobs_p3 = sA_p3 + sB_p3;

% Check difference between p1 and p3
difference_p1_p3 = max(abs(sobs_p1 - sobs_p3));
disp(['Max difference (p1 vs p3): ', num2str(difference_p1_p3)]);

% Plot results
figure;
subplot(1,2,1);
plot(t, sobs_p1, 'b', t, sobs_p2, 'r--');
legend('p1', 'p2 (swapped)');
title('Equivalent Parameter Sets');
xlabel('Time'); ylabel('s_{obs}');

subplot(1,2,2);
plot(t, sobs_p1, 'b', t, sobs_p3, 'g--');
legend('p1', 'p3 (non-equiv)');
title('Non-Equivalent Parameter Sets');
xlabel('Time'); ylabel('s_{obs}');
% Compute K for each region and parameter set
% For p1
K1_A = (c ./ (kon1_A * c + koff1_A)) .* (1 - exp(-(kon1_A * c + koff1_A) * t));
K1_B = (c ./ (kon1_B * c + koff1_B)) .* (1 - exp(-(kon1_B * c + koff1_B) * t));
sum_alphaK_p1 = alpha1_A * K1_A + alpha1_B * K1_B;

% For p2
K2_A = (c ./ (kon2_A * c + koff2_A)) .* (1 - exp(-(kon2_A * c + koff2_A) * t));
K2_B = (c ./ (kon2_B * c + koff2_B)) .* (1 - exp(-(kon2_B * c + koff2_B) * t));
sum_alphaK_p2 = alpha2_A * K2_A + alpha2_B * K2_B;

figure;
subplot(1,3,1);
plot(t, alpha1_A*K1_A, 'b', t, alpha1_B*K1_B, 'r', t, sum_alphaK_p1, 'k--');
legend('α₁K₁ (A)', 'α₁K₁ (B)', 'Sum');
title('p1: αK Contributions');

subplot(1,3,2);
plot(t, alpha2_A*K2_A, 'b', t, alpha2_B*K2_B, 'r', t, sum_alphaK_p2, 'k--');
legend('α₂K₂ (A)', 'α₂K₂ (B)', 'Sum');
title('p2: αK Contributions');

subplot(1,3,3);
plot(t, sum_alphaK_p1, 'b', t, sum_alphaK_p2, 'r--');
legend('Sum p1', 'Sum p2');
title('Equivalence of Sums');
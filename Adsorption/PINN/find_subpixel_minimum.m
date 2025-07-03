% --- NOVA FUNÇÃO AUXILIAR para encontrar o mínimo com precisão sub-pixel ---
function sub_pixel_angle = find_subpixel_minimum(rp_curve, angle_range)
% Encontra o mínimo de uma curva SPR com precisão sub-pixel usando interpolação parabólica.

    % 1. Encontre o índice do mínimo discreto (como antes)
    [~, idx_min] = min(rp_curve);

    % 2. Verificação de segurança: Se o mínimo estiver nas bordas, não podemos interpolar.
    % Nesse caso, apenas retorne o mínimo discreto.
    if idx_min == 1 || idx_min == length(rp_curve)
        sub_pixel_angle = angle_range(idx_min);
        return;
    end
    sub_pixel_angle = angle_range(idx_min)
%     % 3. Pegue os 3 pontos ao redor do mínimo para formar a parábola
%     x_points = angle_range(idx_min-1 : idx_min+1);
%     y_points = rp_curve(idx_min-1 : idx_min+1);
% 
%     % 4. Ajuste um polinômio de 2º grau (uma parábola) a esses 3 pontos.
%     % p = [a, b, c] para a parábola y = ax^2 + bx + c
%     p = polyfit(x_points, y_points, 2);
%     a = p(1);
%     b = p(2);
% 
%     % 5. O vértice (mínimo) de uma parábola é encontrado em x = -b / (2a)
%     sub_pixel_angle = -b / (2 * a);
end
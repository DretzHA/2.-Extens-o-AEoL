function aoi = analytic_iir(k, delta)
% ANALYTIC_IIR Calcula o AoI teórico para estratégia IIR (Zero-Wait)
% Baseado na Equação (7) de Yates et al. (2017).
%
% Entradas:
%   k     : Número de símbolos necessários para decodificar
%   delta : Probabilidade de apagamento do símbolo (Erasure Probability)
%
% Saída:
%   aoi   : Age of Information médio

    % Eq (7): Delta = (k / (1-delta)) * (3/2 + delta/k)
    term1 = k / (1 - delta);
    term2 = 1.5 + (delta / k);
    
    aoi = term1 * term2;
end
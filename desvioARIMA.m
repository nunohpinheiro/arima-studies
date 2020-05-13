function desvioPadrao = desvioARIMA(B, p, sigmahat, numPrev)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    % Cálculos baseados nas págs.252-256 do livro "Forecasting with UBJ Models", Alan Pankratz
    
    ARcoefs = B(1:p,1);
    MAcoefs = B(p+1:end,1);
    
    desvioPadrao = zeros(1,numPrev);                                        % vector-linha para armazenar um valor de desvio associado a cada previsão
    
    %% CÁLCULO DE ROOT-MEAN-SQUARED ERROR (RMSE)
    
%     n = length(resid);
%     m = length([ARcoefs; MAcoefs]);
%     RMSE = ( 1/(n-m) )* sum( resid.^2 );
    
    %% CÁLCULO DOS COEFS. PSI()
    
    PSIcoefs = zeros(1,numPrev);
    PSIcoefs(1,1) = 1;
    
    % MA() PURO
    
    if (isempty(ARcoefs)==1)
        PSIcoefs = [ 1, -1.*(MAcoefs'), zeros(1,numPrev-length(MAcoefs)) ];     % SEM CERTEZAS
        
        
    % AR() COM/SEM MA()
    
    elseif (isempty(ARcoefs)~=1) && (numPrev>1)
        
        PSIcoefs(1,2) = ARcoefs(1,1) - sum(MAcoefs);        % SEM CERTEZAS
        
        for a=3:1:numPrev
            for b=1:1: min( length(ARcoefs),a-1 )
                PSIcoefs(1,a) = PSIcoefs(1,a) + ARcoefs(b)*PSIcoefs(a-b);
            end
        end
    end
    
    %% CÁLCULO FINAL DO DESVIO-PADRÃO DE CADA FORECAST
    
    for a=1:1:numPrev
        desvioPadrao(1,a) = sigmahat*sqrt( sum( PSIcoefs(1,1:a).^2 ));    % NOTA: sqrt(RMSE) = std(resid)
    end
end


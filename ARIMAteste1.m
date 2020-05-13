function [ y_past, y_predict, resid, lowerConf95, upperConf95, ARcoefs, MAcoefs ] = ARIMAteste1( p, d, q, y_antigo, numPrev )
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    if p > length(y_antigo) || q > length(y_antigo) || p<0 || q<0 || d<0 || nargin ~= 5
        disp('Impossível de calcular, erro de input!')
        y_past = y_antigo; y_predict = []; ARcoefs = []; MAcoefs = [];
        
    else
        %% COLOCAR y_antigo COMO VECTOR-COLUNA
        if size(y_antigo,2) ~= 1
            y_antigo = y_antigo';
        end
        
        %%
        % ISOLAR DATASET INICIAL
        y_past = y_antigo;
        
        % COLOCAR DATASET COM MÉDIA NULA, PARA TRABALHAR
        y = y_antigo - mean(y_antigo);                                                       % ISTO PODE PRECISAR DE MAIS TRANSFORMAÇÕES - PROCESSAMENTO!!! - ATENÇÃO, PODE SER NECESSÁRIO METER MÉDIA NULA!!!
        
        %%
        % CÁLCULO DA AUTOCORRELAÇÃO PARCIAL (PACF) E SIMPLES (ACF)
        
        PACFtotal = parcorr(y_past);                                        % vector-coluna com todos os valores de PACF dos dados de treino, relativamente ao 1º valor do dataset
        ACFtotal = autocorr(y_past);                                        % vector-coluna com todos os valores de ACF dos dados de treino, relativamente ao 1º valor do dataset
        
        
        %% CALCULAR COEFICIENTES AR() E MA()
        
        results = arma_mle( y, p, q, 1 );                                      % MEDIA NULA
        ARcoefs = results.ar;
        MAcoefs = results.ma;
        sigmahat = results.sigma;
        
        resid = results.residuos;
        resid = resid(1,2:end);
        
        B = [ARcoefs; MAcoefs];
        
        % CALCULAR RESÍDUOS
%         resid = ValuesResidues(B, p, d, q, y);               % MEDIA NULA
        
        %% FORMA FINAL - FORECAST
        
        [y_predict, lowerConf95, upperConf95] = ARIMAforecast( B, p, d, q, y, resid, sigmahat, numPrev );  % MEDIA NULA
        
        y_predict = y_predict + mean(y_antigo);
        lowerConf95 = lowerConf95 + mean(y_antigo);
        upperConf95 = upperConf95 + mean(y_antigo);
    end
end

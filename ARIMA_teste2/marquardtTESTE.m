function [Bout, resid, SSR] = marquardtTESTE(B, p, d, q, y_antigo)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    h = ones(length(B),1);                                                  % vector-coluna com PASSO DE CORRECÇÃO a somar em cada optimização dos pesos/coeficientes - ATENÇÃO!!! O h NÃO TEM VALORES DIFERENTES?
    
    SSR = [];                                                               % vector-coluna que vai guardar todos as Sum of Squared Residuals calculadas
    
    tolerancia = 0.001;                                                     % valor máximo de 'h' que já se aceita como convergência
    maxIter = 100;                                                          % nº máximo de iterações aceites para o algoritmo
    iteracao = 1;
    
    while h > tolerancia && iteracao < maxIter
        
        %% CÁLCULO DOS RESÍDUOS PARA OS VALORES DE TREINO, FACE AOS COEFICIENTES PREDEFINIDOS
        
        resid = ValuesResidues(B, p, d, q, y_antigo);                       % vector-coluna com os resíduos entre os valores de treino reais e o modelo ajustado com os coeficientes existentes
        
        %% CÁLCULO DE SSR (SUM OF SQUARE RESIDUALS)
        
        SQRresid = resid.^2;                                                % vector-coluna com cada valor de resíduo elevado ao quadrado
        SSR(iteracao,1) = sum(SQRresid);                                    % soma de todos os quadrados dos resíduos, com valor guardado no vector-coluna SSR
        
        %% CÁLCULO DO PASSO 'h'
        
        B_interno = B + 0.01;                                               % cria-se um vector-coluna 'B_interno' com um incremento nos coefs proposto na pág.217 do "Forecasting with UBJ Models", Alan Pankratz
        
        resid_interno = ValuesResidues(B_interno, p, d, q, y_antigo);       % vector-coluna com os resíduos calculados para esse vector 'B_interno'
        
        derivada = resid - resid_interno;                                   % vector-coluna com a diferença entre os dois vectores-coluna de resíduos, que corresponde a uma derivação parcial que entrará na equação de obtenção de 'h', funcionando como 'variável independente'
        
        % O cálculo de 'h' segue a fórmula:
        % resid = -1 . derivada . h + a , com a='resíduos da fórmula que são minimizados', na pág.216 do "Forecasting with UBJ Models", Alan Pankratz
        % Assim, segue a forma:
        % b = A . x
        %
        % 'h' deverá ser calculado através do método de Least Linear Squares (LLS), com a expressão:
        
        A = diag(-1 * derivada);                                            % o vector 'derivada' tem de ser expresso em forma de matriz diagonal
        h = inv((A')*A) * (A') * resid;                                     % equação normal do método LLS, para melhor estimativa de 'h', na pág.256 do livro "Time Series Analysis: Forecasting and Control", Box, Jenkins
        
        % 'h' (calculado em cima) será vector-coluna com tantos elementos,
        % POIS, HÁ AQUI UM PROBLEMA!!!
        
        % VER ARTIGO MARQUARDT E CENA DE corrcoef DO MATLAB
        
        
        
        iteracao = iteracao + 1;
    end
    
end


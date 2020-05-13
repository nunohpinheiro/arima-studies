function [ y_past, y_predict, ARcoefs, MAcoefs, C ] = ARIMAteste3( p, d, q, y_antigo, numPrev )
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
        
        %%
        % DETERMINAÇÃO DOS VALORES INICIAIS PARA OS COEFICIENTES AR()
        
        if p ~= 0
            
            PACFar = PACFtotal(1:p+1);                                      % vector-coluna com os valores de PACF que importam para o cálculo dos coeficientes AR (depende do valor 'p' escolhido)
            
            A = eye(p,p);                                                   % matriz A das equações de Yule-Walker, vai armazenar valores de PACF na ordem pretendida
            
            for l = 1:1:p
                for c = 1:1:p
                    A(l,c) = PACFar( abs(l-c)+1 );
                end
            end
            
            % CÁLCULO DOS COEFICIENTES AR (FI)
            
            ARcoefs = (inv(A)) * PACFar(2:end);                             % vector-linha coefs AR, pelas equações de Yule-Walker. Não se conta com PACF(1) pq é correlação de y(0) consigo próprio (ou seja, 1)
            
            % CÁLCULO DA VARIÂNCIA ASSOCIADA AOS COEFS AR - ainda não percebi se serve para algo, sem ser diagnóstico...
            
            variancia=0;
            for i=1:1:length(ARcoefs)
                variancia = variancia + ARcoefs(i)*PACFar(i);
            end
            
            % CÁLCULO DA ORDENADA ASSOCIADA À PARTE AR (MIU)
            
            miu = mean(y_past) * (1 - sum(ARcoefs));                        % ATENÇÃO A ISTO! PODE SER PRECISO MUDAR (METER DADOS COM MÉDIA NULA, P.EX.)
            
            % Atenção: não sei se está bem calculado... A fórmula é:
            % miu= MEAN *( 1 - sum(ARcoefs) )
            % MEAN= 'mean of differenced series'
            
        elseif p == 0
            ARcoefs = [];
            miu = mean(y_past);                                             % ATENÇÃO! VER BEM ISTO DO MIU!
        end
        
        %%
        % DETERMINAÇÃO DOS VALORES INICIAIS PARA OS COEFICIENTES MA()
        
        if q ~= 0
            
            % CÁLCULO DOS COEFICIENTES MA (TETA)
            
            if q==1 || q==2
                MAcoefs = -1* ACFtotal(2:q+1);                              % Aplicando o sugerido na pág.307 do "Forecasting with UBJ Models", Alan Pankratz. Exclui-se ACF(1) pq é correlação de y(0) consigo próprio
                
            else
                MAcoefs = 0.1* ones(q,1);
            end
            
        elseif q == 0
            MAcoefs = [];
        end
        
        %%
        % OPTIMIZAÇÃO DOS COEFICIENTES AR() E/OU MA()
        
        if p~=0 || q~=0                                                     % estes coeficientes têm de ser optimizados, para se adequarem aos dados de treino fornecidos. O método utilizado (referido em baixo) está sugerido nas págs.200/209 do "Forecasting with UBJ Models", Alan Pankratz
            
            B = [ARcoefs; MAcoefs];
            
            %
            %             [Bout, resid, SSR] = marquardtTESTE(B,p,d,q,y_antigo);          % Função para aplicar o Marquadt's Compromise (ou Levenberg-Marquardt). É uma junção das vantagens do método do gradiente descendente, com as vantagens do método de Newton
            
            
            resid = ValuesResidues(B, p, d, q, y_antigo);
            
            t = 1:1:length(y_antigo)-2;
            
            res = @(B) -1*y_antigo(t+1) + B(1)*y_antigo(t) ;  %   anonym. funct. for residua
            
            [C,ssq,cnt] = LMFnlsq(res,B,'Display',-1);% without displ. lambda
            
            
            
            
        elseif p==0 && q==0
            % NÃO ACONTECE NADA
        end
        
        %% FORMA FINAL - FORECAST
        
        y_predict = ARIMAforecast( ARcoefs, MAcoefs, p, d, q, y_antigo, resid, numPrev );
    end
end

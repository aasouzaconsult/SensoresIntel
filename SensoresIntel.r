options(max.print=8.5E5) #Numero de linhas - retorno
options(error=recover)   #Debugar

############################################################################################################
# Sao 14400 epocas, cada epoca possue a medicao de temperatura de 52 sensores distribuidos no laboratorio. #
# http://db.lcs.mit.edu/labdata/labdata.html                                                               #
# http://www.stat.wisc.edu/~mchung/teaching/MIA/reading/diffusion.gaussian.kernel.pdf.pdf                  #
############################################################################################################

############################################################################################################
############################################################################################################
# Este trabalho tem como objetivo programar uma função de interpolação com kernel gaussiano
# para estimar a superfície de temperatura de uma época qualquer dos dados. Calcular o valor
# da temperatura em 50 pontos de teste localizados aleatoriamente na área sensoriada. Mostrar
# em um gráfico a superfície de interpolação, localizando os pontos dos sensores e os 50 pontos
# de teste. Para aprender o parâmetro abertura do kernel gaussiano, utilizar as primeiras 300
# épocas dos dados. Os valores de temperatura estimados por essa função de interpolação serão
# considerados os valores de referência (ground truth) para os pontos de teste. Saída: códigos
# dos procedimentos, a função de interpolação e os gráficos do campo sensor para algumas épocas
# escolhidas.

############
## ITEM A ##
############

# Funcao para converter em numerico
converte = function(dados) {
  c<-NULL
  for (i in 1:ncol(dados))
    c<-cbind(c,as.numeric(dados[,i]))
  
  return (c)
}

# Importando os dados - ja convertendo
# Usaremos os dados pré-processados, disponíveis em http://www.ulb.ac.be/di/labo/code/PCAgExpe.zip.

dados = converte(read.table("C:/temp/subsfin.txt"))	
locs  = converte(read.table("C:/temp/mote_locs.txt")[,2:3]) # Posicionamento do sensor
locs  = locs[c(-5,-15),]

# Função para calcular o parâmetro de abertura com base nas 300 primeiras épocas
estima_sigma = function(dados){
  sigma = c() # Vetor Sigma
  g_t = dados[1:300,]
  for(i in 1:ncol(dados)){
    sigma[i] = sd(g_t[,i]) # Desvio Padrão para cada um dos 52 sensores (sd(dados[,1]))
  }
  return (sigma)
}

# Gerar pontos (na area) para validacao
gera_pontos = function(n){
  # Cria matriz de 0 com N linhas e 2 colunas
  pontos = matrix(0, nrow = n, ncol=2)
  
  # Gera n pontos entre o mínimo e o máximo de x(locs[,1]) e de y (locs[,2])
  pontos[,1] = sample(min(locs[,1]):max(locs[,1]), n, replace = TRUE) 
  pontos[,2] = sample(min(locs[,2]):max(locs[,2]), n, replace = TRUE) 	
  return (pontos)
}

# Funcao Kernel (Fórmula 2)
kernel_gauss = function(locs , x, sigma){ 
  
  m_x = t(replicate(52, x)) # x é cada um dos pontos gerados 
  # m_x = (t(replicate(52, pontos[1,])))
  # transposto da repetição de cada um dos pontos gerados (Matriz - 52X2)
  
  # Distancia Euclidiana
  dst = (m_x - locs)^2 # Distancia de cada ponto gerado para com os demais (originais)
  dst = dst[,1]+dst[,2]
  
  pt1 = 1/(2*pi*(sigma^2))
  pt2 = exp(-(dst/(2*(sigma^2))))
  result = pt1*pt2
  
  return(result)
}

############
# EXECUCAO #
############
epoca = 301 # Época Testada

sigma = estima_sigma(dados) # Desvio padrão de cada um dos 52 sensores
dados_epoca = dados[epoca,] # 52 temperaturas da Época Testada
pontos = gera_pontos(50)    # Gera 50 pontos na area

result = c() # cria um vetor de resultados

# Chamada da função para estimar as temperaturas nos 50 pontos gerados
for(i in 1: nrow(pontos)){
  # (Fórmula 2 - Kernel Gaussiano)
  # Entrada: 52 Pontos originais(x.y), Pontos gerados (cada um - repetido 52x), sigma (desvio padrão dos 52 pontos)  
  k = kernel_gauss(locs, pontos[i,], sigma) # k são 52 números (exemplo: kernel_gauss(locs, pontos[1,], sigma))
  
  # (Fórmula 1 - Nadaraya-Watson)
  result[i] = (k%*%dados_epoca)/sum(k) # (Os 52 k´s gerados * as 52 temperaturas da Época Testada) / pelo somatório dos 52 k´s
}

# Temperatura em cada um dos 50 pontos escolhidos aleatoriamente 
result

##############################################
# Concatenando as localizacoes e Temperatura #
##############################################
# Vetores de pontos (x,y) e de temperaturas (tn) 
x  = c(locs[,1], pontos[,1]) # x (102 pontos - 52 originais e 50 gerados)
y  = c(locs[,2], pontos[,2]) # y (102 pontos - 52 originais e 50 gerados)
tn = c(dados_epoca, result) # Temperatura (102 temperaturas - 52 originais e 50 gerados)

# Matriz (102x4) concatenando todos os dados
# Coluna 1: Valores de x (52 originais e 50 gerados)
# Coluna 2: Valores de y (52 originais e 50 gerados)
# Coluna 3: Temperaturas (52 originais e 50 gerados)
# Coluna 4: (1 para dados originais e 2 para dados gerados)
coord = matrix(0, ncol = 4, nrow=102) 

coord[1:52,1]   = locs[,1]    # x dos dados originais
coord[53:102,1] = pontos[,1]  # x dos 50 pontos - Previstos
coord[1:52,2]   = locs[,2]    # y dos dados originais
coord[53:102,2] = pontos[,2]  # y dos 50 pontos - Previstos
coord[1:52,3]   = dados_epoca # temperaturas originais
coord[53:102,3] = result      # temperaturas previstas
coord[1:52,4]   = 1           # Originais 
coord[53:102,4] = 2           # 50 pontos - Previstos


# Visão fixa dos pontos
library(scatterplot3d)
scatterplot3d(coord[,1:3],
              main="3D Sensores - Disposição",
              xlab = "x",
              ylab = "y",
              zlab = "z")

#Gerando a superfície (Malha)
library(akima)
library(rgl)
shape = interp(locs[,1], locs[,2], dados_epoca,
               xo=seq(min(locs[,1]), max(locs[,1]), length=600), 
               yo=seq(min(locs[,2]), max(locs[,2]), length=600))
# Esse interp, ele recebe os locs, e as temperaturas, daí aplica pra uma 
#sequencia de tamanho 600 entre o o max e min de cada coordenada (pra ficar 
#dentro do laboratorio), aí gera os valores (600 triplas) com base nos que passei 
#pra construir a malha ou superficie.

# Visão 3D
# Dados Originais (invertido y por z)
plot3d(coord[1:52,1], coord[1:52,3], coord[1:52,2], ylim=c(0,40),
       type = "s",
       col = "black", 
       size = 1,
       xlab = "x",
       ylab = "z",
       zlab = "y")

# Desenhando os pontos gerados (invertido y por z)
plot3d(coord[53:102,1], coord[53:102,3], coord[53:102,2], ylim=c(0,40), 
       type = "s", 
       col = "red", 
       size = 1,
       xlab = "x",
       ylab = "z",
       zlab = "y")

# Plotando a superficie
rgl.surface(shape$x,shape$y,shape$z, color = "orange", alpha=c(0.5))

############################################################################################################
############################################################################################################
# A estratégia (decisão descentralizada) para economizar energia consiste em cada nó sensor
# decidir autonomamente transmitir ou não transmitir sua medição de cada época com probabilidade
# p. Nesse caso, em cada época o nó central conhecerá os valores reais da temperatura
# apenas nos locais dos nós que transmitiram. Para estimar a temperatura nos pontos de teste
# ele utiliza os últimos valores disponíveis ou estimados nos pontos dos nós sensores. Todos os
# nós sensores são programados para transmitirem as primeiras 5 épocas. Assim, os dados dessas
# épocas estarão disponíveis no nó central.
# (b) Estimar os valores de temperatura nos pontos de teste utilizando regressão linear com
# regularização. Ajustar experimentalmente o parâmetro de regularização. Calcular, sobre todo o
# dataset: o NRMSE (raiz do erro quadrático médio normalizado) de épocas, os NRMSEs máximos
# e mínimo, e reter o id e os valores das 10 épocas com NRMSEs máximos e 10 épocas com NRMSEs
# mínimos. Fazer isso para p= 0.1:0.1:0.9.
# Saída: Gráfico e tabela com NRMSEs, médio, máximo e mínimo em função de p.

############
## ITEM B ##
############

# Função para calcular a probabilidade (probabilidade de quem vai ou não transmitir)
# Quanto maior a probabilidade, maior o numero de sensores (o "melhor")

trans_p = function(p){
  #aqui colocar 52 ou 102 no numero de colunas
  colunas = 52
  # Matriz 14400 x n colunas(informado na linha anterior)
  result = matrix(0, nrow = nrow(dados), ncol = colunas)
  
  for(i in 1:nrow(dados)){
    i_tr = c() # vetor de probabilidades usadas
    count = 0
    
    for (j in 1:ncol(dados)) {
      if(runif(1)<p){ # Numero aleatório é menor que a probabilidade de entrada?
        count = count + 1
        i_tr[count] = j
      }
      if(length(i_tr)<=1){ # numero baixo de sensores (probabilidades mais baixas)
        i_tr = sample(1:52, 2)
      }
    }
    #opcional para incluir pontos de teste
    #c = gera_pontos(50)		
    #c = rbind(locs, c)
    #caso não...
    c = locs
    # Calcula a regressão Linear com base nos sensores escolhidos
    result[i,] = reg_lin(i_tr, dados[i,], c)
  }
  return(result)
}

## Testando a função acima (trans_p)
# colunas = 52
# Alex = matrix(0, nrow = nrow(dados), ncol = colunas)
# 
# for(i in 1:nrow(dados)){
#   i_tr = c() # vetor de probabilidades usadas
#   count = 0
#   
#   for (j in 1:ncol(dados)) {
#     if(runif(1)<0.8){ # Numero aleatório é menor que a probabilidade de entrada?
#       count = count + 1
#       i_tr[count] = j
#     }
#     if(length(i_tr)<=1){ # numero baixo de sensores
#       i_tr = sample(1:52, 2)
#     }
#   }
#   c = locs
#   Alex[i,] = reg_lin(i_tr, dados[i,], c)
# }
# summary(Alex)

# Função da Regressão Linear
reg_lin = function(i_tr, dados_epoca, nlocs){
  result = array(0, dim = nrow(nlocs))        # array de nrow(nlocs) posições
  X = matrix(1, ncol = 3, nrow = nrow(nlocs)) # Matriz de 1´s (nrow(nlocs) linhas e 3 colunas)
 #X[,1] é o bias
  X[,2] = nlocs[,1]
  X[,3] = nlocs[,2]   
  #...
  # Faz a regressão apenas se teve mudanças (sensores a menos)
  if(nrow(nlocs) != length(dados_epoca) || length(i_tr) != length(dados_epoca)){
    x = X[i_tr,] # pontos que foram escolhidos
    lambda = 0.1	
    I = diag(ncol(x)) # 3x3
    
    # Calculo do B (Beta -> w)? v-> w
    v = solve(t(x)%*%x + lambda*I)%*%t(x)%*%dados_epoca[i_tr] # dados_epoca[i_tr] - Temperaturas na época Testada # 3x1
    
    # length(result), esse for vai preencher o vetor resultado...verificando se foi transmitido e copiando pro result, 
    #e se não, estimando.
    for(i in 1:length(result)){
      if(i %in% i_tr){ # Se tiver entre os "trasmitidos" na probabilidade (pega os dados da época)
        result[i] = dados_epoca[i]
      }
      else{ # Senão, estima.
        # y = XB + e?
        result[i] = (v[1,1] + v[2,1]*nlocs[i,1] + v[3,1]*nlocs[i,2])
      }
    }
  }
  else{
    return (dados_epoca)
  }
  return (result) # Retorna as temperaturas estimadas (pontos escolhidos - probabilidade) e as reais dos pontos não escolhidos
}

# Testando a regressão...
#  Alex = array(0, dim = nrow(c)) # array de nrow(nlocs) posições
#  X = matrix(1, ncol = 3, nrow = nrow(c)) # Matriz de 1´s (nrow(nlocs) linhas e 3 colunas)
#  X[,2] = nlocs[,1]
#  X[,3] = nlocs[,2]    
#  #...
#  # Faz a regressão apenas se teve mudanças (sensores a menos)
#  if(nrow(c) != length(dados_epoca) || length(i_tr) != length(dados_epoca)){
#    x = X[i_tr,]	# pontos que foram escolhidos
#    lambda = 0.1	
#    I = diag(ncol(x)) # 3x3
#    v = solve(t(x)%*%x + lambda*I)%*%t(x)%*%dados_epoca[i_tr] # 3x1
#    for(i in 1:length(result)){
#      if(i %in% i_tr){
#        result[i] = dados_epoca[i]
#      }
#      else{		
#        result[i] = (v[1,1] + v[2,1]*nlocs[i,1] + v[3,1]*nlocs[i,2])
#      }
#    }
#  }
##  else{
##    return (dados_epoca)
##  }
##  return (result)
##}

# NRMSE - Erro máximo
erromax = function(result, dados){
  re2_max_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    re2_max_ep[i] = max(sqrt(a))
  }
  return (re2_max_ep)
}

# NRMSE - Erro médio
erromedio = function(result, dados){
  rmse_epoca = c()
  re2_min_ep = c()
  re2_max_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    rmse_epoca[i] = sqrt(sum(a))/ncol(result)
  }
  return (rmse_epoca)
}

# NRMSE - Erro mínimo
erromin = function(result, dados){
  rmse_epoca = c()
  re2_min_ep = c()
  re2_max_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    re2_min_ep[i] = min(sqrt(a))
  }
  
  return (re2_min_ep)
}

# Função para reter os 10 NRMSEs máximos
max10 = function(re2_max_ep){
  res2 = re2_max_ep
  
  ind = c()
  for(i in 1:10){
    a = which.max(res2) # Índice do maior valor
    ind[i] = a
    res2[a] = 0
  }
  return(ind)
}
# Função para reter os 10 NRMSEs mínimos
min10 = function(re2_min_ep){
  
  res2 = re2_min_ep
  
  ind = c()
  for(i in 1:10){
    #a = which(mm == min(re2_min_ep), arr.ind = TRUE)
    a = which.min(res2)
    ind[i] = a
    res2[a] = 1000
  }
  return(ind)
}

##############
# Executando #
##############

# Gerando 50 pontos
a = gera_pontos(50)

# Monta a regressão com probabilidade de 0.1 a 0.9
rs1 = trans_p(0.1)	
rs2 = trans_p(0.2)
rs3 = trans_p(0.3)
rs4 = trans_p(0.4)
rs5 = trans_p(0.5)
rs6 = trans_p(0.6)
rs7 = trans_p(0.7)
rs8 = trans_p(0.8)
rs9 = trans_p(0.9)

# Plotando os erros máximos
 plot(erromax(rs1, dados), type="l", col="blue", ylim=c(0,400), xlab = "", ylab = "" )
lines(erromax(rs2, dados), type="l", pch=22, lty=2, col="yellow")
lines(erromax(rs3, dados), type="l", pch=22, lty=2, col="black")
lines(erromax(rs4, dados), type="l", pch=22, lty=2, col="red")
lines(erromax(rs5, dados), type="l", pch=22, lty=2, col="green")
lines(erromax(rs6, dados), type="l", pch=22, lty=2, col="pink")
lines(erromax(rs7, dados), type="l", pch=22, lty=2, col="orange")
lines(erromax(rs8, dados), type="l", pch=22, lty=2, col="gray")
lines(erromax(rs9, dados), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros máximos")

# Plotando os erros médios
 plot(erromedio(rs1, dados), type="l", col="blue", ylim=c(0,30), xlab = "", ylab = "" )
lines(erromedio(rs2, dados), type="l", pch=22, lty=2, col="yellow")
lines(erromedio(rs3, dados), type="l", pch=22, lty=2, col="black")
lines(erromedio(rs4, dados), type="l", pch=22, lty=2, col="red")
lines(erromedio(rs5, dados), type="l", pch=22, lty=2, col="green")
lines(erromedio(rs6, dados), type="l", pch=22, lty=2, col="pink")
lines(erromedio(rs7, dados), type="l", pch=22, lty=2, col="orange")
lines(erromedio(rs8, dados), type="l", pch=22, lty=2, col="gray")
lines(erromedio(rs9, dados), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros médios")

# Plotando os erros mínimos
 plot(erromin(rs1, dados), type="l", col="blue", ylim=c(0,0.1), xlab = "", ylab = "" )
lines(erromin(rs2, dados), type="l", pch=22, lty=2, col="yellow")
lines(erromin(rs3, dados), type="l", pch=22, lty=2, col="black")
lines(erromin(rs4, dados), type="l", pch=22, lty=2, col="red")
lines(erromin(rs5, dados), type="l", pch=22, lty=2, col="green")
lines(erromin(rs6, dados), type="l", pch=22, lty=2, col="pink")
lines(erromin(rs7, dados), type="l", pch=22, lty=2, col="orange")
lines(erromin(rs8, dados), type="l", pch=22, lty=2, col="gray")
lines(erromin(rs9, dados), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros mínimos")

# Gerando dados para a Tabela Máximo
max_matrix     = matrix(0, nrow = 9, ncol=10)
max_matrix[1,] = min10(erromax(rs1, dados))
max_matrix[2,] = min10(erromax(rs2, dados))
max_matrix[3,] = min10(erromax(rs3, dados))
max_matrix[4,] = min10(erromax(rs4, dados))
max_matrix[5,] = min10(erromax(rs5, dados))
max_matrix[6,] = min10(erromax(rs6, dados))
max_matrix[7,] = min10(erromax(rs7, dados))
max_matrix[8,] = min10(erromax(rs8, dados))
max_matrix[9,] = min10(erromax(rs9, dados))
View(max_matrix)

max_matrix     = matrix(0, nrow = 9, ncol=10)
max_matrix[1,] = max10(erromax(rs1, dados))
max_matrix[2,] = max10(erromax(rs2, dados))
max_matrix[3,] = max10(erromax(rs3, dados))
max_matrix[4,] = max10(erromax(rs4, dados))
max_matrix[5,] = max10(erromax(rs5, dados))
max_matrix[6,] = max10(erromax(rs6, dados))
max_matrix[7,] = max10(erromax(rs7, dados))
max_matrix[8,] = max10(erromax(rs8, dados))
max_matrix[9,] = max10(erromax(rs9, dados))
View(max_matrix)

# Gerando dados para a Tabela Médio
med_matrix     = matrix(0, nrow = 9, ncol=10)
med_matrix[1,] = min10(erromedio(rs1, dados))
med_matrix[2,] = min10(erromedio(rs2, dados))
med_matrix[3,] = min10(erromedio(rs3, dados))
med_matrix[4,] = min10(erromedio(rs4, dados))
med_matrix[5,] = min10(erromedio(rs5, dados))
med_matrix[6,] = min10(erromedio(rs6, dados))
med_matrix[7,] = min10(erromedio(rs7, dados))
med_matrix[8,] = min10(erromedio(rs8, dados))
med_matrix[9,] = min10(erromedio(rs9, dados))
View(med_matrix)

med_matrix     = matrix(0, nrow = 9, ncol=10)
med_matrix[1,] = max10(erromedio(rs1, dados))
med_matrix[2,] = max10(erromedio(rs2, dados))
med_matrix[3,] = max10(erromedio(rs3, dados))
med_matrix[4,] = max10(erromedio(rs4, dados))
med_matrix[5,] = max10(erromedio(rs5, dados))
med_matrix[6,] = max10(erromedio(rs6, dados))
med_matrix[7,] = max10(erromedio(rs7, dados))
med_matrix[8,] = max10(erromedio(rs8, dados))
med_matrix[9,] = max10(erromedio(rs9, dados))
View(med_matrix)

# Gerando dados para a Tabela Mínimo
min_matrix     = matrix(0, nrow = 9, ncol=10)
min_matrix[1,] = min10(erromin(rs1, dados))
min_matrix[2,] = min10(erromin(rs2, dados))
min_matrix[3,] = min10(erromin(rs3, dados))
min_matrix[4,] = min10(erromin(rs4, dados))
min_matrix[5,] = min10(erromin(rs5, dados))
min_matrix[6,] = min10(erromin(rs6, dados))
min_matrix[7,] = min10(erromin(rs7, dados))
min_matrix[8,] = min10(erromin(rs8, dados))
min_matrix[9,] = min10(erromin(rs9, dados))
View(min_matrix)

min_matrix     = matrix(0, nrow = 9, ncol=10)
min_matrix[1,] = max10(erromin(rs1, dados))
min_matrix[2,] = max10(erromin(rs2, dados))
min_matrix[3,] = max10(erromin(rs3, dados))
min_matrix[4,] = max10(erromin(rs4, dados))
min_matrix[5,] = max10(erromin(rs5, dados))
min_matrix[6,] = max10(erromin(rs6, dados))
min_matrix[7,] = max10(erromin(rs7, dados))
min_matrix[8,] = max10(erromin(rs8, dados))
min_matrix[9,] = max10(erromin(rs9, dados))
View(min_matrix)

# Testes extras #
# Concatena os 50 pontos dos sensores originais com os 50 pontos gerados
# c = rbind(locs, a)
# Calculando a Regressão Linear em uma época
# reg_lin(seq(1,52), dados[1,], c)

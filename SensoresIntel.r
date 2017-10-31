#options(max.print=8.5E5) #Numero de linhas - retorno
#options(error=recover)   #Debugar

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

  # transposto da repetição de cada um dos pontos gerados (Matriz - 52X2)
  m_x = t(replicate(52, x)) # x é cada um dos pontos gerados 
  
  # Distancia Euclidiana (eleva ao quadrado e depois tira a raiz - tirar os negativos)
  dst = (m_x - locs)^2         # Distancia de cada ponto gerado para os originais (52)
  dst = sqrt(dst[,1]+dst[,2])  # Raiz quadrada
  
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

# Chamada da função para estimar as temperaturas nos 50 pontos gerados (ground truth)
result = c() # cria um vetor de resultados
for(i in 1: nrow(pontos)){
  # (Fórmula 2 - Kernel Gaussiano)
  # Entrada: 52 Pontos originais(x.y), Pontos gerados (cada um - repetido 52x), sigma (desvio padrão dos 52 sensores)  
  k = kernel_gauss(locs, pontos[i,], sigma) # k são 52 números (exemplo: kernel_gauss(locs, pontos[1,], sigma))
  
  # (Fórmula 1 - Nadaraya-Watson)
  result[i] = sum(k%*%dados_epoca)/sum(k) # (Os 52 k´s gerados * as 52 temperaturas da Época Testada) / pelo somatório dos 52 k´s
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
# Função da Regressão Linear
reg_lin = function(i_tr, dados_epoca, nlocs){
  result = array(0, dim = nrow(nlocs))       # array de 0 de nrow(nlocs) + 50 posições
  X = matrix(1, ncol = 3, nrow = nrow(locs)) # Matriz de 1´s (nrow(nlocs) linhas e 3 colunas)
  #X[,1] é o bias  
  X[,2] = locs[,1]
  X[,3] = locs[,2]    
  #...
  # Faz a regressão apenas se teve mudanças (sensores a menos)
  if(nrow(nlocs) != length(dados_epoca) || length(i_tr) != length(dados_epoca)){
    x = X[i_tr,]	# pontos que foram escolhidos
    lambda = 0.1	
    I = diag(ncol(x)) # 3x3
    
    # Fórmula 3 (modelo linear com regularização)
    v = solve(t(x)%*%x + lambda*I)%*%t(x)%*%dados_epoca[i_tr]
    
    # length(result), esse for vai preencher o vetor resultado...verificando se foi transmitido e copiando pro result, 
    #e se não, estimando.    
    for(i in 1:length(result)){
      if(i %in% i_tr){  # Se tiver entre os "trasmitidos" na probabilidade (pega os dados da época)
        result[i] = dados_epoca[i]
      }
      else{	# Senão, estima
        # y = XB + e - Fórmula 4 (cálculo da estimação)
        result[i] = (v[1,1] + v[2,1]*nlocs[i,1] + v[3,1]*nlocs[i,2])
      }
    }
  }
  else{
    return (dados_epoca)
  }
  
  return (result) # Retorna as temperaturas estimadas (ponto não escolhido - probabilidade) e as reais dos pontos escolhidos
}

# Função para calcular a probabilidade (probabilidade de quem vai ou não transmitir)
# Quanto maior a probabilidade, maior o numero de sensores (o "melhor")

trans_p = function(p, pontos){
  #matriz 14400 x 50 colunas(pontos gerados)
  result = matrix(0, nrow = nrow(dados), ncol = nrow(pontos))
  
  for(i in 1:nrow(dados)){
    i_tr = c() # vetor de probabilidades usadas
    count = 0
    
    for (j in 1:ncol(dados)) {
      if(runif(1)<p){ # Numero aleatório é menor que a probabilidade de entrada
        count = count + 1
        i_tr[count] = j
      }
      if(length(i_tr)<=1){ # numero baixo de sensores (probabilidades mais baixas)
        i_tr = sample(1:52, 2)
      }
    }
    
    #print(length(i_tr))
    
    # Calcula a regressão Linear com base nos sensores escolhidos (50)
    # Para cada Epoca i, retorna um vetor de 50 temperaturas (sensores escolhidos ou estimadas)
    result[i,] = reg_lin(i_tr, dados[i,], pontos) 
  }
  return(result)
}
result

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
   #a = sqrt((dados[i,] - result[i,])^2)
   #rmse_epoca[i] = sum(a)/ncol(result)
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
    a = which.max(res2) # índice do maior valor
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
sigma = estima_sigma(dados)
pontos = gera_pontos(50)

# Chamada da função para estimar as temperaturas de cada época (14400) nos 50 pontos gerados 
# (Kernel Gaussiano)  (ground truth)
resultT = matrix(0, nrow = nrow(dados), ncol = nrow(pontos))
for(i in 1:nrow(resultT)){
  for(j in 1:ncol(resultT)){
    k = kernel_gauss(locs, pontos[j,], sigma) # Distancia Euclidiana
    resultT[i,j] = (k%*%dados[i,])/sum(k)
  }
}
summary(resultT)

# Monta a regressão linear com probabilidade de 0.1 a 0.9 para as 50 temperaturas para cada época
rs1 = trans_p(0.1,pontos)	
rs2 = trans_p(0.2,pontos)
rs3 = trans_p(0.3,pontos)
rs4 = trans_p(0.4,pontos)
rs5 = trans_p(0.5,pontos)
rs6 = trans_p(0.6,pontos)
rs7 = trans_p(0.7,pontos)
rs8 = trans_p(0.8,pontos)
rs9 = trans_p(0.9,pontos)

# Numero de sensores (exemplo)
# 0.1 - 3
# 0.2 - 7
# 0.3 - 19
# 0.4 - 22
# 0.5 - 26
# 0.6 - 35
# 0.7 - 36
# 0.8 - 41
# 0.9 - 46

#############
## TABELAS ##
#############
# Gerando dados para a Tabela Máximo
max_matrix     = matrix(0, nrow = 9, ncol=10)
max_matrix[1,] = min10(erromax(rs1, resultT))
max_matrix[2,] = min10(erromax(rs2, resultT))
max_matrix[3,] = min10(erromax(rs3, resultT))
max_matrix[4,] = min10(erromax(rs4, resultT))
max_matrix[5,] = min10(erromax(rs5, resultT))
max_matrix[6,] = min10(erromax(rs6, resultT))
max_matrix[7,] = min10(erromax(rs7, resultT))
max_matrix[8,] = min10(erromax(rs8, resultT))
max_matrix[9,] = min10(erromax(rs9, resultT))
View(max_matrix)

max_matrix     = matrix(0, nrow = 9, ncol=10)
max_matrix[1,] = max10(erromax(rs1, resultT))
max_matrix[2,] = max10(erromax(rs2, resultT))
max_matrix[3,] = max10(erromax(rs3, resultT))
max_matrix[4,] = max10(erromax(rs4, resultT))
max_matrix[5,] = max10(erromax(rs5, resultT))
max_matrix[6,] = max10(erromax(rs6, resultT))
max_matrix[7,] = max10(erromax(rs7, resultT))
max_matrix[8,] = max10(erromax(rs8, resultT))
max_matrix[9,] = max10(erromax(rs9, resultT))
View(max_matrix)

# Gerando dados para a Tabela Médio
med_matrix     = matrix(0, nrow = 9, ncol=10)
med_matrix[1,] = min10(erromedio(rs1, resultT))
med_matrix[2,] = min10(erromedio(rs2, resultT))
med_matrix[3,] = min10(erromedio(rs3, resultT))
med_matrix[4,] = min10(erromedio(rs4, resultT))
med_matrix[5,] = min10(erromedio(rs5, resultT))
med_matrix[6,] = min10(erromedio(rs6, resultT))
med_matrix[7,] = min10(erromedio(rs7, resultT))
med_matrix[8,] = min10(erromedio(rs8, resultT))
med_matrix[9,] = min10(erromedio(rs9, resultT))
View(med_matrix)

med_matrix     = matrix(0, nrow = 9, ncol=10)
med_matrix[1,] = max10(erromedio(rs1, resultT))
med_matrix[2,] = max10(erromedio(rs2, resultT))
med_matrix[3,] = max10(erromedio(rs3, resultT))
med_matrix[4,] = max10(erromedio(rs4, resultT))
med_matrix[5,] = max10(erromedio(rs5, resultT))
med_matrix[6,] = max10(erromedio(rs6, resultT))
med_matrix[7,] = max10(erromedio(rs7, resultT))
med_matrix[8,] = max10(erromedio(rs8, resultT))
med_matrix[9,] = max10(erromedio(rs9, resultT))
View(med_matrix)

# Gerando dados para a Tabela Mínimo
min_matrix     = matrix(0, nrow = 9, ncol=10)
min_matrix[1,] = min10(erromin(rs1, resultT))
min_matrix[2,] = min10(erromin(rs2, resultT))
min_matrix[3,] = min10(erromin(rs3, resultT))
min_matrix[4,] = min10(erromin(rs4, resultT))
min_matrix[5,] = min10(erromin(rs5, resultT))
min_matrix[6,] = min10(erromin(rs6, resultT))
min_matrix[7,] = min10(erromin(rs7, resultT))
min_matrix[8,] = min10(erromin(rs8, resultT))
min_matrix[9,] = min10(erromin(rs9, resultT))
View(min_matrix)

min_matrix     = matrix(0, nrow = 9, ncol=10)
min_matrix[1,] = max10(erromin(rs1, resultT))
min_matrix[2,] = max10(erromin(rs2, resultT))
min_matrix[3,] = max10(erromin(rs3, resultT))
min_matrix[4,] = max10(erromin(rs4, resultT))
min_matrix[5,] = max10(erromin(rs5, resultT))
min_matrix[6,] = max10(erromin(rs6, resultT))
min_matrix[7,] = max10(erromin(rs7, resultT))
min_matrix[8,] = max10(erromin(rs8, resultT))
min_matrix[9,] = max10(erromin(rs9, resultT))
View(min_matrix)

####################################
## PLOTANDO ERROS TODAS AS EPOCAS ##
####################################
# Plotando os erros máximos
 plot(erromax(rs1, resultT), type="l", col="blue", ylim=c(0,400), xlab = "", ylab = "" )
lines(erromax(rs2, resultT), type="l", pch=22, lty=2, col="yellow")
lines(erromax(rs3, resultT), type="l", pch=22, lty=2, col="black")
lines(erromax(rs4, resultT), type="l", pch=22, lty=2, col="red")
lines(erromax(rs5, resultT), type="l", pch=22, lty=2, col="green")
lines(erromax(rs6, resultT), type="l", pch=22, lty=2, col="pink")
lines(erromax(rs7, resultT), type="l", pch=22, lty=2, col="orange")
lines(erromax(rs8, resultT), type="l", pch=22, lty=2, col="gray")
lines(erromax(rs9, resultT), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros máximos")

# Plotando os erros médios
 plot(erromedio(rs1, resultT), type="l", col="blue", ylim=c(0,220), xlab = "", ylab = "" )
lines(erromedio(rs2, resultT), type="l", pch=22, lty=2, col="yellow")
lines(erromedio(rs3, resultT), type="l", pch=22, lty=2, col="black")
lines(erromedio(rs4, resultT), type="l", pch=22, lty=2, col="red")
lines(erromedio(rs5, resultT), type="l", pch=22, lty=2, col="green")
lines(erromedio(rs6, resultT), type="l", pch=22, lty=2, col="pink")
lines(erromedio(rs7, resultT), type="l", pch=22, lty=2, col="orange")
lines(erromedio(rs8, resultT), type="l", pch=22, lty=2, col="gray")
lines(erromedio(rs9, resultT), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros médios")

# Plotando os erros mínimos
 plot(erromin(rs1, resultT), type="l", col="blue", ylim=c(0,3.5), xlab = "", ylab = "" )
lines(erromin(rs2, resultT), type="l", pch=22, lty=2, col="yellow")
lines(erromin(rs3, resultT), type="l", pch=22, lty=2, col="black")
lines(erromin(rs4, resultT), type="l", pch=22, lty=2, col="red")
lines(erromin(rs5, resultT), type="l", pch=22, lty=2, col="green")
lines(erromin(rs6, resultT), type="l", pch=22, lty=2, col="pink")
lines(erromin(rs7, resultT), type="l", pch=22, lty=2, col="orange")
lines(erromin(rs8, resultT), type="l", pch=22, lty=2, col="gray")
lines(erromin(rs9, resultT), type="l", pch=22, lty=2, col="brown")
title("Plotando os erros mínimos")

########################################################
## GRAFICOS TOTAIS - MAX, MEAN e MIN dos Erros Médios ##
########################################################
# Máximo dos erros médios
mx = c()
mx[1] = max(erromedio(rs1, resultT))
mx[2] = max(erromedio(rs2, resultT))
mx[3] = max(erromedio(rs3, resultT))
mx[4] = max(erromedio(rs4, resultT))
mx[5] = max(erromedio(rs5, resultT))
mx[6] = max(erromedio(rs6, resultT))
mx[7] = max(erromedio(rs7, resultT))
mx[8] = max(erromedio(rs8, resultT))
mx[9] = max(erromedio(rs9, resultT))
plot(mx, type="l", col="blue", ylim=c(0,180), xlab = "", ylab = "" )
title("Plotando os erros máximos")

# Média dos erros médios
m = c()
m[1] = mean(erromedio(rs1, resultT))
m[2] = mean(erromedio(rs2, resultT))
m[3] = mean(erromedio(rs3, resultT))
m[4] = mean(erromedio(rs4, resultT))
m[5] = mean(erromedio(rs5, resultT))
m[6] = mean(erromedio(rs6, resultT))
m[7] = mean(erromedio(rs7, resultT))
m[8] = mean(erromedio(rs8, resultT))
m[9] = mean(erromedio(rs9, resultT))
plot(m, type="l", col="blue", ylim=c(0,5), xlab = "", ylab = "" )
title("Plotando os erros médios")

# Minimo dos erros médios
n = c()
n[1] = min(erromedio(rs1, resultT))
n[2] = min(erromedio(rs2, resultT))
n[3] = min(erromedio(rs3, resultT))
n[4] = min(erromedio(rs4, resultT))
n[5] = min(erromedio(rs5, resultT))
n[6] = min(erromedio(rs6, resultT))
n[7] = min(erromedio(rs7, resultT))
n[8] = min(erromedio(rs8, resultT))
n[9] = min(erromedio(rs9, resultT))
plot(n, type="l", col="blue", ylim=c(0,0.5), xlab = "", ylab = "" )
title("Plotando os erros mínimos")

#####################################################################
plot(mx, type="l", col="red"  , ylim=c(0,180), xlab = "", ylab = "" )
lines(m, type="l", col="blue" , ylim=c(0,5)  , xlab = "", ylab = "" )
lines(n, type="l", col="green", ylim=c(0,0.5), xlab = "", ylab = "" )
title("Plotando os erros")


################################################################
## GRAFICOS TOTAIS - MEAN dos Erros máximos, médios e mínimos ##
################################################################
# Máximo dos erros médios
mx1 = c()
mx1[1] = mean(erromax(rs1, resultT))
mx1[2] = mean(erromax(rs2, resultT))
mx1[3] = mean(erromax(rs3, resultT))
mx1[4] = mean(erromax(rs4, resultT))
mx1[5] = mean(erromax(rs5, resultT))
mx1[6] = mean(erromax(rs6, resultT))
mx1[7] = mean(erromax(rs7, resultT))
mx1[8] = mean(erromax(rs8, resultT))
mx1[9] = mean(erromax(rs9, resultT))
plot(mx1, type="l", col="blue", ylim=c(0,30), xlab = "", ylab = "" )
title("Plotando os erros max dos máximos")

# Média dos erros médios
m1 = c()
m1[1] = mean(erromedio(rs1, resultT))
m1[2] = mean(erromedio(rs2, resultT))
m1[3] = mean(erromedio(rs3, resultT))
m1[4] = mean(erromedio(rs4, resultT))
m1[5] = mean(erromedio(rs5, resultT))
m1[6] = mean(erromedio(rs6, resultT))
m1[7] = mean(erromedio(rs7, resultT))
m1[8] = mean(erromedio(rs8, resultT))
m1[9] = mean(erromedio(rs9, resultT))
plot(m1, type="l", col="blue", ylim=c(0,2.3), xlab = "", ylab = "" )
title("Plotando os max dos erros médios")

# Minimo dos erros médios
n1 = c()
n1[1] = mean(erromin(rs1, resultT))
n1[2] = mean(erromin(rs2, resultT))
n1[3] = mean(erromin(rs3, resultT))
n1[4] = mean(erromin(rs4, resultT))
n1[5] = mean(erromin(rs5, resultT))
n1[6] = mean(erromin(rs6, resultT))
n1[7] = mean(erromin(rs7, resultT))
n1[8] = mean(erromin(rs8, resultT))
n1[9] = mean(erromin(rs9, resultT))
plot(n1, type="l", col="blue", ylim=c(0,1), xlab = "", ylab = "" )
title("Plotando os max dos erros mínimos")

# Variancia dos erros médios
v1 = c()
v1[1] = var(erromedio(rs1, resultT))
v1[2] = var(erromedio(rs2, resultT))
v1[3] = var(erromedio(rs3, resultT))
v1[4] = var(erromedio(rs4, resultT))
v1[5] = var(erromedio(rs5, resultT))
v1[6] = var(erromedio(rs6, resultT))
v1[7] = var(erromedio(rs7, resultT))
v1[8] = var(erromedio(rs8, resultT))
v1[9] = var(erromedio(rs9, resultT))
plot(v1, type="l", col="blue", ylim=c(0,2.3), xlab = "", ylab = "" )
title("Plotando os max dos erros médios")

#####################################################################
 plot(mx1  , type="l", col="red"  , ylim=c(0,15) , xlab = "", ylab = "" )
lines(m1   , type="l", col="blue" , ylim=c(0,5)  , xlab = "", ylab = "" )
lines(n1   , type="l", col="green", ylim=c(0,0.5), xlab = "", ylab = "" )
lines(m1+v1, type="l", col="black", ylim=c(0,0.5), xlab = "", ylab = "", lty=2)
lines(m1-v1, type="l", col="black", ylim=c(0,0.5), xlab = "", ylab = "", lty=2)
title("Média dos Erros máximos, médios e mínimos (item b)")


###############################################################################
# Repetir o item (b) utilizando regressão kernel gaussiano (Nadaraya-Watson). #
###############################################################################

############
## ITEM C ##
############

# Funcao Kernel (Fórmula 2) - para o Nadaraya Watson
kernel_gauss_nw = function(i_tr , dados_epoca, sigma, ponto){
  
 #Ao inves de 52x - Agora replica pelo numero de sensores escolhidos
  m_x = t(replicate(length(i_tr), ponto))
  dst = (m_x - locs[i_tr,])^2
  dst = sqrt(dst[,1]+dst[,2])
  
  pt1 = 1/(2*pi*(sigma[i_tr]^2))
  pt2 = exp(-(dst/(2*(sigma[i_tr]^2)))) # sigma dos pontos escolhidos
  result = pt1*pt2
  
  return(result)
}

# Função - Nadaraya Watson (Formula 1)
nadaraya = function(i_tr, dados_epoca, sigma, ponto){	
  # Mudou a funcao, para pegar informações dos sensores escolhidos
  k = kernel_gauss_nw(i_tr, dados_epoca[i_tr], sigma, ponto) # saída = qtd num. de sensores escolhidos
  return( (k%*%dados_epoca[i_tr])/sum(k) ) # saída = 1 temperatura ("ponto a ponto")
}

# Função para calcular a probabilidade (probabilidade de quem vai ou não transmitir)
# Quanto maior a probabilidade, maior o numero de sensores (o "melhor")
trans_p_nw = function(p, pontos){
  sigma = estima_sigma(dados)
  result = matrix(0, nrow = nrow(dados), ncol = nrow(pontos))
  for(i in 1:nrow(dados)){
    i_tr = c()
    count = 0
    for(j in 1:nrow(pontos)) { # De cada época a probabilidade de cada ponto 
      if(runif(1)<p){
        count = count + 1
        i_tr[count] = j
      }
    }
    if(length(i_tr)<=1){
      i_tr = sample(1:52, 2)
    }
    #print(i)
    for(d in 1:nrow(pontos)){
      result[i,d] = nadaraya(i_tr, dados[i,], sigma, pontos[d,]) # saida = 1 temperatura
    }
  }
  return(result)
}

################
## Executando ##
################
sigma  = estima_sigma(dados)
pontos = gera_pontos(50) 

# Chamada da função para estimar as temperaturas de cada época nos 50 pontos gerados s(ground truth)
resultT2 = matrix(0, nrow = nrow(dados), ncol = nrow(pontos))
for(i in 1:nrow(resultT2)){
  for(j in 1:ncol(resultT2)){
    k = kernel_gauss(locs, pontos[j,], sigma)
    resultT2[i,j] = (k%*%dados[i,])/sum(k)
  }
  #print(i)
}

# Monta a regressão kernel gaussiano com probabilidade de 0.1 a 0.9 para as 50 temperaturas para cada época
# Nadaraya Watson
rs1 = trans_p_nw(0.1,pontos)	
rs2 = trans_p_nw(0.2,pontos)
rs3 = trans_p_nw(0.3,pontos)
rs4 = trans_p_nw(0.4,pontos)
rs5 = trans_p_nw(0.5,pontos)
rs6 = trans_p_nw(0.6,pontos)
rs7 = trans_p_nw(0.7,pontos)
rs8 = trans_p_nw(0.8,pontos)
rs9 = trans_p_nw(0.9,pontos)

#############################################################
## GRAFICOS TOTAIS - MEAN, MIN, MAX e VAR dos Erros Médios ##
#############################################################
# Vetor de media dos erros médios
m = c()
m[1] = mean(erromedio(rs1, resultT2))
m[2] = mean(erromedio(rs2, resultT2))
m[3] = mean(erromedio(rs3, resultT2))
m[4] = mean(erromedio(rs4, resultT2))
m[5] = mean(erromedio(rs5, resultT2))
m[6] = mean(erromedio(rs6, resultT2))
m[7] = mean(erromedio(rs7, resultT2))
m[8] = mean(erromedio(rs8, resultT2))
m[9] = mean(erromedio(rs9, resultT2))

# Vetor de minimos dos erros médios
n = c()
n[1] = min(erromedio(rs1, resultT2))
n[2] = min(erromedio(rs2, resultT2))
n[3] = min(erromedio(rs3, resultT2))
n[4] = min(erromedio(rs4, resultT2))
n[5] = min(erromedio(rs5, resultT2))
n[6] = min(erromedio(rs6, resultT2))
n[7] = min(erromedio(rs7, resultT2))
n[8] = min(erromedio(rs8, resultT2))
n[9] = min(erromedio(rs9, resultT2))

# Vetor de maximos dos erros médios
o = c()
o[1] = max(erromedio(rs1, resultT2))
o[2] = max(erromedio(rs2, resultT2))
o[3] = max(erromedio(rs3, resultT2))
o[4] = max(erromedio(rs4, resultT2))
o[5] = max(erromedio(rs5, resultT2))
o[6] = max(erromedio(rs6, resultT2))
o[7] = max(erromedio(rs7, resultT2))
o[8] = max(erromedio(rs8, resultT2))
o[9] = max(erromedio(rs9, resultT2))

# Vetor de variancias dos erros médios
v = c()
v[1] = var(erromedio(rs1, resultT2))
v[2] = var(erromedio(rs2, resultT2))
v[3] = var(erromedio(rs3, resultT2))
v[4] = var(erromedio(rs4, resultT2))
v[5] = var(erromedio(rs5, resultT2))
v[6] = var(erromedio(rs6, resultT2))
v[7] = var(erromedio(rs7, resultT2))
v[8] = var(erromedio(rs8, resultT2))
v[9] = var(erromedio(rs9, resultT2))

##########################
# Plotando os resultados #
##########################

 plot(m  , type="l", col="blue", ylim=c(0,2.1), xlab = "Probabilidade", ylab = "Erros" )
lines(n  , type="l", pch=10, col="Orange", xlab = "", ylab = "" )
lines(o  , type="l", pch=22, col="red"   , xlab = "", ylab = "" )
lines(m+v, type="l", pch=22, lty=2, col="black" , xlab = "", ylab = "" )
lines(m-v, type="l", pch=22, lty=2, col="black" , xlab = "", ylab = "" )
title("MEAN, MIN, MAX e VAR dos Erros Médios (Item C)")

####################################################################################
## GRAFICOS TOTAIS - MEAN dos Erros Maximo, Medios e Minimos + Variancia do Medio ##
####################################################################################
# Vetor max dos erros maximos
m2 = c()
m2[1] = mean(erromax(rs1, resultT2))
m2[2] = mean(erromax(rs2, resultT2))
m2[3] = mean(erromax(rs3, resultT2))
m2[4] = mean(erromax(rs4, resultT2))
m2[5] = mean(erromax(rs5, resultT2))
m2[6] = mean(erromax(rs6, resultT2))
m2[7] = mean(erromax(rs7, resultT2))
m2[8] = mean(erromax(rs8, resultT2))
m2[9] = mean(erromax(rs9, resultT2))

# Vetor max dos erros médios
n2 = c()
n2[1] = mean(erromedio(rs1, resultT2))
n2[2] = mean(erromedio(rs2, resultT2))
n2[3] = mean(erromedio(rs3, resultT2))
n2[4] = mean(erromedio(rs4, resultT2))
n2[5] = mean(erromedio(rs5, resultT2))
n2[6] = mean(erromedio(rs6, resultT2))
n2[7] = mean(erromedio(rs7, resultT2))
n2[8] = mean(erromedio(rs8, resultT2))
n2[9] = mean(erromedio(rs9, resultT2))

# Vetor max dos erros minimos
o2 = c()
o2[1] = mean(erromin(rs1, resultT2))
o2[2] = mean(erromin(rs2, resultT2))
o2[3] = mean(erromin(rs3, resultT2))
o2[4] = mean(erromin(rs4, resultT2))
o2[5] = mean(erromin(rs5, resultT2))
o2[6] = mean(erromin(rs6, resultT2))
o2[7] = mean(erromin(rs7, resultT2))
o2[8] = mean(erromin(rs8, resultT2))
o2[9] = mean(erromin(rs9, resultT2))

# Vetor de variancias dos erros médios
v2 = c()
v2[1] = var(erromedio(rs1, resultT2))
v2[2] = var(erromedio(rs2, resultT2))
v2[3] = var(erromedio(rs3, resultT2))
v2[4] = var(erromedio(rs4, resultT2))
v2[5] = var(erromedio(rs5, resultT2))
v2[6] = var(erromedio(rs6, resultT2))
v2[7] = var(erromedio(rs7, resultT2))
v2[8] = var(erromedio(rs8, resultT2))
v2[9] = var(erromedio(rs9, resultT2))

##########################
# Plotando os resultados #
##########################
 plot(m2  , type="l" , pch=22, col="blue", ylim=c(0,5)  , xlab = "Probabilidade", ylab = "Erros" )
lines(n2  , type="l" , pch=22, col="Orange", xlab = "", ylab = "" )
lines(o2  , type="l" , pch=22, col="red"   , xlab = "", ylab = "" )
lines(n2+v2, type="l", pch=22, col="black" , xlab = "", ylab = "", lty=2)
lines(n2-v2, type="l", pch=22, col="black" , xlab = "", ylab = "", lty=2)
title("MEAN dos Erros Máximos, Médios e Mínimos + Variância dos Médios (Item C)")

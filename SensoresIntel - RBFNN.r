#options(max.print=8.5E5) #Numero de linhas - retorno
#options(error=recover)   #Debugar

############################################################################################################
# Sao 14400 epocas, cada epoca possue a medicao de temperatura de 52 sensores distribuidos no laboratorio. #
# http://db.lcs.mit.edu/labdata/labdata.html                                                               #
# http://www.stat.wisc.edu/~mchung/teaching/MIA/reading/diffusion.gaussian.kernel.pdf.pdf                  #
############################################################################################################
# Sensores Intel - RBFNN

# Função para converter em numérico
converte = function(dados) {
  c<-NULL
  for (i in 1:ncol(dados))
    c<-cbind(c,as.numeric(dados[,i]))
  
  return (c)
}

# Importando os dados - já convertendo
# Usaremos os dados pré-processados, disponíveis em http://www.ulb.ac.be/di/labo/code/PCAgExpe.zip.
dados = converte(read.table("C:/temp/subsfin.txt"))	
locs  = converte(read.table("C:/temp/mote_locs.txt")[,2:3]) # Posicionamento do sensor
locs  = locs[c(-5,-15),]

d_tr=dados[1,]

# Apenas 1 Épocas e 10 Sensores
# Y = dados[1,1:10]
# X  = locs[1:10,]

# NRMSE - Erro máximo
erromax = function(result, dados){
  re2_max_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    re2_max_ep[i] = sqrt(max(a))
  }
  return (re2_max_ep)
}

# NRMSE - Erro médio
erromed = function(result, dados){
  re2_med_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    re2_med_ep[i] = sqrt(mean(a))
  }
  return (re2_med_ep)
}

# NRMSE - Erro mínimo
erromin = function(result, dados){
  re2_min_ep = c()
  for(i in 1:nrow(result)){
    a = (dados[i,] - result[i,])^2
    re2_min_ep[i] = sqrt(min(a))
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

# Função para calcular a probabilidade (probabilidade de quem vai ou não transmitir)
# Quanto maior a probabilidade, maior o numero de sensores que irão transmitir ("melhor")

trans_p = function(p, pontos){
  d_tr = dados[1,] # Primeira linha das temperaturas (base)
  
  result = matrix(0, nrow = nrow(dados), ncol = nrow(pontos))
  #matriz 14400 x 50 colunas(pontos gerados)
  
  # Treinando o y
  ytreino = dados * 0
  ytreino[1,] = d_tr
  
  for(i in 2:nrow(dados)){
    for (j in 1:nrow(locs)) {
      if(runif(1)<p){                # Se o ponto for escolhido (será transmitido)
        d_tr[j] = dados[i,j]         # Atualiza os dados base (d_tr)
      }
    }
    ytreino[i,]	= d_tr
  }
  
  k = kmeans(locs, 10)               # K-Means, com k = 10
  
  w = treino_rbf(ytreino, 10, k)     # Treinando
  
  a = teste_rbf(pontos, k, 10, w)    # Testando com base nos pontos gerados
  
  return(a)
}

# Teste da RBF
teste_rbf = function(pontos, k, n_neuronios, w){
  print("AH")
  H        = matrix(1, nrow=nrow(pontos), ncol = n_neuronios+1) # 50X11 - 1 Bias
  means    = k$centers
  clusters = k$cluster
  v = array(0, c(2,2, n_neuronios))
  for(i in 1:n_neuronios){
    group = which(clusters %in% i) # Unir os do mesmo grupo
    v[,,i] = var(locs[group,])	   # variancia do grupo
  }
  
  for(i in 1:nrow(pontos)){
    for(j in 1:n_neuronios){
      sigma = diag(1,2)
      diag(sigma) = diag(v[,,j])
      #H[i,j+1] = exp(-(rowSums((locs[i,]-means)^2)/var%*%var))
      H[i,j+1] = exp(-(t(pontos[i,]-means[j,]) %*% solve(sigma) %*% (pontos[i,]-means[j,])))
    }
  }
  # yT
  return(t(H%*%w))
}

# Treinamento da RBF
treino_rbf = function(y, n_neuronios, k){	
  H        = matrix(1, nrow=nrow(locs), ncol = n_neuronios+1) # Matriz de 1 - 52x11 (1 - Bias)
  means    = k$centers
  clusters = k$cluster

 # Variancia dos grupos  
  v = array(0, c(2,2, n_neuronios))
  for(i in 1:n_neuronios){
    group = which(clusters %in% i)
    v[,,i] = var(locs[group,])	
  }
  
  # H
  for(i in 1:nrow(locs)){
    for(j in 1:n_neuronios){
      #H[i,j+1] = exp(-(rowSums((locs[i,]-means)^2)/var%*%var))
      sigma = diag(1,2)
      diag(sigma) = diag(v[,,j])
      #print(sigma)
      H[i,j+1] = exp(-(t(locs[i,]-means[j,]) %*% solve(sigma) %*% (locs[i,]-means[j,])))
    }
  }
  
  # W
  return(MASS::ginv(t(H)%*%H)%*%t(H)%*%t(y)) # 
}

# Gerar pontos (na area) para validacao
gera_pontos = function(n){
  pontos = matrix(0, nrow = n, ncol=2) # Cria matriz de 0 com N linhas e 2 colunas
  # Gera n pontos entre o mínimo e o máximo de x(locs[,1]) e de y (locs[,2])
  pontos[,1] = sample(min(locs[,1]):max(locs[,1]), n, replace = TRUE) 
  pontos[,2] = sample(min(locs[,2]):max(locs[,2]), n, replace = TRUE) 	
  return (pontos)
}

# Função para calcular o parâmetro de abertura com base nas 300 primeiras épocas
estima_sigma = function(dados){
  sigma = c()
  g_t = dados[1:300,]
  for(i in 1:ncol(dados)){
    sigma[i] = sd(g_t[,i]) # Desvio Padrão para cada um dos 52 sensores (sd(dados[,1]))
  }
  return (sigma)
}

############
# EXECUCAO #
############
sigma  = estima_sigma(dados)
pontos = gera_pontos(50)

# Funcao Kernel (Fórmula 2)
kernel_gauss_gt = function(locs , x, sigma){
  # transposto da repetição de cada um dos pontos gerados (Matriz - 52X2)
  m_x = t(replicate(52, x))   # x é cada um dos pontos gerados
  
  # Distancia Euclidiana (eleva ao quadrado e depois tira a raiz - tirar os negativos)
  dst = (m_x - locs)^2        # Distancia de cada ponto gerado para os originais (52)
  dst = sqrt(dst[,1]+dst[,2]) # Raiz quadrada
  
  pt1 = 1/(2*pi*(sigma^2))
  pt2 = exp(-(dst/(2*(sigma^2))))
  result = pt1*pt2
  
  return(result)
}

# Chamada da função Kernel - para estimar as temperaturas nos 50 pontos gerados (ground truth)
resultT = matrix(0, nrow = nrow(dados), ncol = nrow(pontos)) # 14400 X 50
for(i in 1:nrow(resultT)){
  for(j in 1:ncol(resultT)){
    k = kernel_gauss_gt(locs, pontos[j,], sigma)  # Distancia Euclidiana
    resultT[i,j] = (k%*%dados[i,])/sum(k)
  }
  print(i)
}

# Monta a RBF com probabilidade de 0.1 a 0.9 para as 50 temperaturas para cada época
rs1 = trans_p(0.1,pontos)	
rs2 = trans_p(0.2,pontos)
rs3 = trans_p(0.3,pontos)
rs4 = trans_p(0.4,pontos)
rs5 = trans_p(0.5,pontos)
rs6 = trans_p(0.6,pontos)
rs7 = trans_p(0.7,pontos)
rs8 = trans_p(0.8,pontos)
rs9 = trans_p(0.9,pontos)

# Plotando os erros
#############################################################
## GRAFICOS TOTAIS - MEAN, MIN, MAX e VAR dos Erros Médios ##
#############################################################
# Vetor de media dos erros médios
m = c()
m[1] = mean(erromed(rs1, resultT))
m[2] = mean(erromed(rs2, resultT))
m[3] = mean(erromed(rs3, resultT))
m[4] = mean(erromed(rs4, resultT))
m[5] = mean(erromed(rs5, resultT))
m[6] = mean(erromed(rs6, resultT))
m[7] = mean(erromed(rs7, resultT))
m[8] = mean(erromed(rs8, resultT))
m[9] = mean(erromed(rs9, resultT))

# Vetor de minimos dos erros minimos
n = c()
n[1] = min(erromin(rs1, resultT))
n[2] = min(erromin(rs2, resultT))
n[3] = min(erromin(rs3, resultT))
n[4] = min(erromin(rs4, resultT))
n[5] = min(erromin(rs5, resultT))
n[6] = min(erromin(rs6, resultT))
n[7] = min(erromin(rs7, resultT))
n[8] = min(erromin(rs8, resultT))
n[9] = min(erromin(rs9, resultT))

# Vetor de maximos dos erros médios
o = c()
o[1] = max(erromax(rs1, resultT))
o[2] = max(erromax(rs2, resultT))
o[3] = max(erromax(rs3, resultT))
o[4] = max(erromax(rs4, resultT))
o[5] = max(erromax(rs5, resultT))
o[6] = max(erromax(rs6, resultT))
o[7] = max(erromax(rs7, resultT))
o[8] = max(erromax(rs8, resultT))
o[9] = max(erromax(rs9, resultT))

# Vetor de variancias dos erros médios
v = c()
v[1] = var(erromed(rs1, resultT))
v[2] = var(erromed(rs2, resultT))
v[3] = var(erromed(rs3, resultT))
v[4] = var(erromed(rs4, resultT))
v[5] = var(erromed(rs5, resultT))
v[6] = var(erromed(rs6, resultT))
v[7] = var(erromed(rs7, resultT))
v[8] = var(erromed(rs8, resultT))
v[9] = var(erromed(rs9, resultT))

##########################
# Plotando os resultados #
##########################
plot (m  , type="l", col="blue", ylim=c(0,29), xlab = "Probabilidade", ylab = "Erros" )
lines(n  , type="l", pch=10, col="Orange", xlab = "", ylab = "" ) # Erro Minimo
lines(o  , type="l", pch=22, col="red"   , xlab = "", ylab = "" ) # Erro Maximo
lines(m+v, type="l", pch=22, lty=2, col="black" , xlab = "", ylab = "" ) # Erro Medio
lines(m-v, type="l", pch=22, lty=2, col="black" , xlab = "", ylab = "" ) # Erro Medio
title("MEAN, MIN, MAX e VAR dos Erros Médios (Item D - RBF)")

############
## Testes ##
############

d_tr = dados[1,] # Primeira linha das temperaturas (base)

# Treinando o y
ytreino = dados * 0
ytreino[1,] = d_tr

for(i in 2:nrow(dados)){
  for (j in 1:nrow(locs)) {
    if(runif(1)<0.1){                # Se o ponto for escolhido (será transmitido)
      d_tr[j] = dados[i,j]         # Atualiza os dados base (d_tr)
    }
  }
  ytreino[i,]	= d_tr
}

# Plotando - Analises
plot(dados[,1],ytreino[,1])
plot(dados[,14]) # Outliers
plot(dados[,7])
plot(dados[,39])

# Outliers
i = 1
j = 1

for (j in 1:ncol(dados)){
  print('Epoca: ')
  print(j)
  
  for (i in 1:nrow(dados)){
    med = mean(dados[,j])
    var =  var(dados[,j])
    std = med - var
    
    if (dados[i,j] < std){
      #Outliers[i,j] = dados[i,j]
      print(dados[i, j])
    } 
  }
}

# Observações
# - Colocar o K pra fora do trans_p
# - Outliers - uma solução: Suavização exponencial

############################################################################################################
# Sao 14400 epocas, cada epoca possue a medicao de temperatura de 52 sensores distribuidos no laboratorio. #
# http://db.lcs.mit.edu/labdata/labdata.html                                                               #
# http://www.stat.wisc.edu/~mchung/teaching/MIA/reading/diffusion.gaussian.kernel.pdf.pdf                  #
############################################################################################################

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
locs = converte(read.table("C:/temp/mote_locs.txt")[,2:3]) # Posicionamento do sensor
locs=locs[c(-5,-15),]

# Função para calcular o parametro de abertura
estima_sigma = function(dados){
  sigma = c()
  g_t = dados[1:300,]
  for(i in 1:ncol(dados)){
    sigma[i] = sd(g_t[,i]) # desvio padrão
  }
  return (sigma)
}

# Gerar pontos (na area) para validacao
gera_pontos = function(n){
  pontos = matrix(0, nrow = n, ncol=2)
  pontos[,1] = sample(min(locs[,1]):max(locs[,1]), n, replace = TRUE) 
  pontos[,2] = sample(min(locs[,2]):max(locs[,2]), n, replace = TRUE) 	
  return (pontos)
}

# Funcao Kernel
kernel_gauss = function(locs , x, sigma){
  
  m_x = t(replicate(52, x)) # x é o ponto
  
  dst = (m_x - locs)^2
  dst = dst[,1]+dst[,2] # Distancia Euclidiana
  
  pt1 = 1/(2*pi*(sigma^2))
  pt2 = exp(-(dst/(2*(sigma^2))))
  result = pt1*pt2
  
  return(result)
}

############
# EXECUCAO #
############
epoca = 301 # Época 301 - Testada

sigma = estima_sigma(dados) # para cada sensor
dados_epoca = dados[epoca,]
pontos = gera_pontos(50)

result = c()
for(i in 1: nrow(pontos)){
  k = kernel_gauss(locs, pontos[i,], sigma)
  result[i] = (k%*%dados_epoca)/sum(k)
}

# Temperatura em cada ponto escolhido aleatoriamente
result

##############################################
# Concatenando as localizacoes e Temperatura #
##############################################

x = c(locs[,1], pontos[,1]) # x
y = c(locs[,2], pontos[,2]) # y
tn = c(dados_epoca, result) # Temperatura

# Matriz com todos os dados
coord = matrix(0, ncol = 4, nrow=102) # Matriz de Zeros 102X4
coord[1:52,4]   = 1           # Originais 
coord[53:102,4] = 2           # 50 pontos - Previstos
coord[1:52,1]   = locs[,1]    # x dos dados originais
coord[1:52,2]   = locs[,2]    # y dos dados originais
coord[1:52,3]   = dados_epoca # temperaturas originais
coord[53:102,1] = pontos[,1]  # x dos 50 pontos - Previstos
coord[53:102,2] = pontos[,2]  # y dos 50 pontos - Previstos
coord[53:102,3] = result      # temperaturas previstas

# Visão fixa dos pontos
library(scatterplot3d)
scatterplot3d(coord[,1:3],
              main="3D Sensores - Disposição",
              xlab = "x",
              ylab = "y",
              zlab = "z")

library(akima)
library(rgl)
#Gerando a superfície
shape = interp(locs[,1], locs[,2], dados_epoca,
               xo=seq(min(locs[,1]), max(locs[,1]), length=600), 
               yo=seq(min(locs[,2]), max(locs[,2]), length=600))

# Visão 3D
# Dados Originais
plot3d(coord[1:52,1], coord[1:52,3], coord[1:52,2], ylim=c(0,40),
       type = "s",
       col = "black", 
       size = 1,
       xlab = "x",
       ylab = "z",
       zlab = "y")

# Desenhando os pontos gerados
plot3d(coord[53:102,1], coord[53:102,3], coord[53:102,2], ylim=c(0,40), 
       type = "s", 
       col = "red", 
       size = 1,
       xlab = "x",
       ylab = "z",
       zlab = "y")

# Plotando a superficie
rgl.surface(shape$x,shape$y,shape$z, color = "orange", alpha=c(0.5))

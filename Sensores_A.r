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

locsTot <- matrix(nrow=102,ncol=4)

# Localizacoes originais (locs)
for(i in 1:nrow(locs)) { 
  locsTot[i,] = matrix(c(locs[i,],dados[epoca,i],1)) # 1 Originais
}

#localizacoes geradas (pontos)
j = 1
for(i in 53:nrow(locsTot)) {
  locsTot[i,] = matrix(c(pontos[j,],result[j],2)) # 2 Resultados
  j = j+1
}

locsTot # Localizações concatenadas + Temperatura
locsTot[,1:3] # Sem identificacao

################################
# Concatenando as temperaturas #
################################

TempTot <- matrix(nrow=1,ncol=102)

dados[epoca,52]

# Temperaturas originais (dados)
for(i in 1:ncol(dados)) { 
  TempTot[i] = matrix(c(dados[epoca,i]))
}

# Temperaturas geradas (result)
j = 1
i = 53
for(i in 53:ncol(TempTot)) {
  TempTot[i] = matrix(c(result[j]))
  j = j+1
}

TempTot # Temperaturas concatenadas


#########################
### PLOTAR SUPERFICIE ###
#########################

# Gerando x para Matrix
# x <- matrix(nrow=1,ncol=102)
#for(i in 1:ncol(x)) { 
#  x[1,i] = matrix(c(locsTot[i,1]))
#}
#x = locsTot[,1]

# Gerando y para Matrix
#y <- matrix(nrow=1,ncol=102)
#for(i in 1:ncol(y)) { 
#  y[1,i] = matrix(c(locsTot[i,2]))
#}
#y = locsTot[,2]

#z = TempTot #z

scatterplot3d(locsTot[,1:3],
              main="3D Sensores - Disposição",
              xlab = "x",
              ylab = "y",
              zlab = "z")

#library(plotly)
#kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
#p <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
#Desenho <- plot_ly(x = locsTot[,1], y = locsTot[,2], z = TempTot) %>% add_surface(z)

# Required for using persp3D() function below.
#library(plot3D)
## We call persp3D function for same Gaussian kernal data generated above.
#persp3D(locsTot[,1:3],theta=30, phi=50, axes=TRUE,scale=2, box=TRUE, nticks=5, 
#        ticktype="detailed",xlab="X-value", ylab="Y-value", zlab="Z-value", 
#        main="Gaussian Kernal with persp3D()")

#kd <- with(MASS::geyser, MASS::kde2d(locsTot[,1], locsTot[,2], n = 102))

# SURFACE
#p <- plot_ly(x = locsTot[,2], y = locsTot[,1], z = TempTot) %>% add_surface(z*2)


################
# Novos testes #
################

# Originais
n = 52
x = locs[,1]
y = locs[,2]
z = dados[epoca,]
z1 = dados[epoca:352,] # 52 após o 301

locais = matrix(nrow=52,ncol=52) # t tem que ser algo 52X52
for(i in 1:ncol(locais)) { 
  locais[i,] = z
}

# Gerados
x1 = pontos[,1]
y1 = pontos[,2]
z1 = result

# Desenhando os pontos originais
ponts = plot3d(x, y, z, type = "s", col = "red", size = 1, forceClipregion = TRUE, xlim=c(0,40),
               ylim=c(0,32),
               zlim=c(0,50))

# Desenhando os pontos gerados
#pontsG = plot3d(x1, y1, z1, type = "s", col = "black", size = 1, forceClipregion = TRUE, xlim=c(0,40),
#                ylim=c(0,32),
#                zlim=c(0,50))

surface3d(x, y, locais, back = 'line', front = 'line', col = 'black', lwd = 1.0, alpha = 0.4)
#axes3d()

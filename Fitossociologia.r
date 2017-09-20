#----------------------------------------------------------------------------------------------
# Fun��o para c�lculo dos descritores fitossociol�gicos e similaridade entre s�tios
# 
# Reposit�rio: https://github.com/ricds/fitoR
#
# Contribui��es: Alexandre Gabriel Christo, Pedro Higuchi, Ricardo Dal'Agnol, Arthur Vinicius Rodrigues
#
#
# Observa��es:
# a) O arquivo de entrada precisa das seguintes colunas com os nomes exatos em minusculo:
#	- parc (identifica��o das parcelas, numerico ou letra, tanto faz)
#	- spp (nome vulgar ou cient�fico do indiv�duo, ou utilizacao de rotulos - sp1, sp2, etc)
#	- dap (di�metro a altura do peito)
#	ou
#	- cap (circunfer�ncia a altura do peito)
#	e quando tiver algum indiv�duo com mais de um fuste/rebrota, criar colunas com nome:
#	- cap ou dap1, dap2, dap3, dap4, etc.
#	OBS: o script buscar� pelas palavras-chave 'cap' e 'dap' nas colunas para localiza-los.
#	Portanto, os dados n�o podem possuir nenhuma coluna extra com 'cap' e 'dap' contidos em
#	seu nome, ou causar� diferen�a no c�lculo, por ex.: uma coluna chamada 'dap_medio'.
#
# b) O script cont�m duas fun��es:
#	1) fitoR(variavel_input, area_de_cada_parcela_em_m2, 'nome_do_arquivo')
#	2) similaridade(variavel_input1, variavel_input2)
#	Para utiliza��o, carregue seus dados, execute as fun��es do IN�CIO at� o FIM das fun��es
#	marcado abaixo com v�rios ###. Para tal selecione o texto e aperte Ctrl+R. Em seguida,
#	utilize a fun��o desejada.
#	Na fun��o fitoR(), o nome_do_arquivo � o nome para o arquivo .csv que cont�m os resultados,
#	colocar entre ' ', ou pode n�o funcionar direito
# 
# c) Mais algumas informa��es da sintaxe do uso no fim do arquivo
#
#-----------------------------------------------------------------------------------------------


# # # # # # # # # # # # # # # # # # # # # # # # 
# IN�CIO DAS FUN��ES - CARREGAR A PARTIR DAQUI
# # # # # # # # # # # # # # # # # # # # # # # # 

# Fun��o fitoR
 
fitoR<-function(x,area,filename)
{
matriz<-table(x$spp,x$parc)

#numero de parcelas
nparc<-length(levels(as.factor(x$parc)))

#area total amostrada
area.parc=(area*nparc)

#densidade
dta=length(x$spp)/(area.parc/10000)

#desvio da densidade entre parcelas
dtadesv = 0
dtai = 1
vetor=1
while(dtai <= nparc)
{
length(vetor) <- nparc
vetor[dtai] = sum(matriz[,dtai])
dtadesv = sd(vetor)/(area/10000)
dtai = dtai + 1
}

#calcula o numero de ind amostrados
N<-apply(matriz,1,sum)

#calcula densidades
DA<-apply(matriz,1,sum)/(area.parc/10000)
DR<-DA/sum(DA)*100

#calcula frequencias 
freq<-(if (length(dim(matriz)) > 1)
 {apply(matriz > 0,1,sum)} else sum(matriz > 0))
FA<-(freq/nparc)*100
FR<-(FA/sum(FA))*100

#checa por NAs nos dados e transforma em zeros
x[is.na(x)] <- 0

#determina se existe "caps" ou "daps" e quais colunas est�o
cols = grep('cap', colnames(x))
ncols = length(cols)
if (ncols>0) param="cap"

cols2 = grep('dap', colnames(x))
ncols2 = length(cols2)
if (ncols2>0) param="dap"

if (param=="dap") cols=cols2
if (param=="dap") ncols=ncols2

#calcula a area da se��o transversal para cada cap/dap e faz a soma por individuo
i=1
x$areasec=0
while (i<=ncols)
{
if (param=="cap") x$areasec<-x$areasec+((pi*(x[,cols[i]]/pi)^2)/40000)
if (param=="dap") x$areasec<-x$areasec+((pi*x[,cols[i]]^2)/40000)
i=i+1
}

#calcula as dominancias
DoA<-tapply(x$areasec, x$spp, sum)/(area.parc/10000)
DoR<-DoA/sum(DoA) *100
 
# area basal por esp�cie
AB<-tapply(x$areasec, x$spp, sum)

#area basal
abta=sum(DoA)
 
#desvio da area basal entre parcelas
somag<-tapply(x$areasec, x$parc, sum)/(area/10000)
abdesv = sd(somag)

#calcula o indice de valor de importancia
VI<-(DR+DoR+FR)/3

#monta a tabela
fito=data.frame(N=N,AB=AB,DA=DA,DR=DR,DoA=DoA,DoR=DoR,FA=FA,FR=FR,VI=VI)
fito$AB<-round(fito$AB,digits=2)
fito$DR<-round(fito$DR,digits=2)
fito$DA<-round(fito$DA,digits=2)
fito$FR<-round(fito$FR,digits=2)
fito$FA<-round(fito$FA,digits=2)
fito$DoR<-round(fito$DoR,digits=2)
fito$DoA<-round(fito$DoA,digits=2)
fito$VI<-round(fito$VI,digits=2)
fito <- fito[order(VI, decreasing = TRUE),]
print(fito)

#calcula os indices de diversidade
Pi<-N/sum(N)
Pi<-Pi*log(Pi)
SW=-sum(Pi)
S=nrow(fito)
J=SW/log(S)

cat("Densidade total por �rea = ",round(dta,digits=2),"\u00B1",round(dtadesv,digits=2),"ind/ha", fill=TRUE)
cat("�rea basal total por �rea = ",round(abta,digits=2),"\u00B1",round(abdesv,digits=2),"m2/ha", fill=TRUE)
cat("Riqueza = ",S,"esp.", fill=TRUE)
cat("�ndice de Shannon-Wiener (H') = ",SW, fill=TRUE)
cat("Equabilidade de Pielou (J) = ",J, fill=TRUE)


if(!missing(filename)) filename = paste(filename, ".csv", sep="") else filename = 'fito.csv'
write.table(fito, file = filename, row.names = TRUE, dec=",", sep=";", quote=FALSE, col.names=NA)
write.table(' ', file = filename, sep=";", quote=TRUE, append=TRUE, row.names=FALSE, col.names=FALSE)
cat("Densidade total por �rea = ",round(dta,digits=2),"\u00B1",round(dtadesv,digits=2),"ind/ha", fill=TRUE, file=filename, append=TRUE)
cat("�rea basal total por �rea = ",round(abta,digits=2),"\u00B1",round(abdesv,digits=2),"m2/ha", fill=TRUE, file=filename, append=TRUE)
cat("Riqueza = ",S,"esp.", fill=TRUE, file=filename, append=TRUE)
cat("�ndice de Shannon-Wiener (H') = ",SW, fill=TRUE, file=filename, append=TRUE)
cat("Equabilidade de Pielou (J) = ",J, fill=TRUE, file=filename, append=TRUE)


}

# Fun��o similaridade

similaridade<-function(x,y)
{
x=levels(x$spp)
y=levels(y$spp)
sim=data.frame(site=1, spp=x)
sim2=data.frame(site=2, spp=y)
sim=rbind(sim,sim2)
table=table(sim$site, sim$spp)
dif=table[1,]-table[2,]

i=1
xspp=0
xyspp=0
yspp=0
while(i<=length(dif))
{
if (dif[i] == 1) xspp = xspp + 1
if (dif[i] == 0) xyspp = xyspp + 1
if (dif[i] == -1) yspp = yspp + 1
i = i + 1
}

jac=(xyspp/(xyspp+xspp+yspp))
sor=((2*xyspp)/((2*xyspp)+xspp+yspp))

cat("A similaridade de Jaccard entre os dois s�tios � de",jac,fill=TRUE)
cat("A similaridade de Sorensen entre os dois s�tios � de",sor,fill=TRUE)
cat("As duas �reas compartilham ",xyspp,"esp., sendo que a primeira tem",xspp,"esp. exclusivas e a segunda tem",yspp,fill=TRUE)

}


# # # # # # # # # # # # # # # # # # # #
# FIM DAS FUN��ES - CARREGAR AT� AQUI 
# # # # # # # # # # # # # # # # # # # #










# APLICANDO AS FUN��ES
 
# Fun��o fitoR()
# Calcula os principais descritores ecol�gicos, �ndices de diversidade, etc.
# Escreve na tela o resultado e em arquivo .csv

# Sintaxe do uso da fun��o fitoR():
# fitoR(variavel_input, area_de_cada_parcela_em_m2, 'nome_do_arquivo')

# Exemplo do uso com os dados teste: area da parcela de 2500 m2 e vai gerar um arquivo "Fito_Exemplo.csv" como output
dados=read.table("Dados_exemplo.csv", header=TRUE, sep=";", dec=",")
fitoR(dados, 2500, 'Fito_Exemplo')


# Fun��o similaridade()
# Calcula a similaridade flor�stica entre os sitios, comparando as esp�cies que os sitios compartilham ou n�o
# �ndices de Jaccard e Sorensen: Baseados na matriz de presen�a-aus�ncia de esp�cies

# Sintaxe do uso da fun��o similaridade()
# similaridade(variavel_input1, variavel_input2)

# Exemplo do uso com o dado teste, obviamente vai resultar em similaridade igual a 1 pois � o mesmo input em ambas as vari�veis
similaridade(dados, dados)

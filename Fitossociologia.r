#----------------------------------------------------------------------------------------------
# Função para cálculo dos descritores fitossociológicos e similaridade entre sítios
#
# Autores originais:  Alexandre Gabriel Christo e Pedro Higuchi - 25/03/2012
#
# Modificações: Ricardo Dal'Agnol da Silva
#
# Data: 06/06/2013
# *Add: Riqueza, Índices Shannon-Wiener e J de Pielou
# *Add: função similaridade(area1,area2)
#
# Data: 10/06/2013
# *Corrigido: equação da área basal
# *Add: desvios da densidade e área basal
# *Add: escreve no .csv a densidade, ab, riq, shannon e pielou
# *Add: testa se existe a variavel cap ou dap, pode-se usar tanto um quanto o outro
#
# Data: 12/06/2013
# *Corrigido: equação de sorensen
# *Add: possibilidade de colocar um nome para o arquivo.csv dos resultados
#
# Data: 17/06/2013
# *Add: contabiliza CAP/DAP de múltiplos fustes no DoA e área basal, max=10
#
# Data: 27/06/2013
# *Add: multiplos fustes não tem mais limite de número (antes era 10)
#
# Observações:
# a) O arquivo de entrada precisa das seguintes colunas com os nomes exatos em minusculo:
#	- parc (identificação das parcelas, numerico ou letra, tanto faz)
#	- spp (nome vulgar ou científico do indivíduo, ou utilizacao de rotulos - sp1, sp2, etc)
#	- dap (diâmetro a altura do peito)
#	ou
#	- cap (circunferência a altura do peito)
#	e quando tiver algum indivíduo com mais de um fuste/rebrota, criar colunas com nome:
#	- cap ou dap1, dap2, dap3, dap4, etc.
#	OBS: o script buscará pelas palavras-chave 'cap' e 'dap' nas colunas para localiza-los.
#	Portanto, os dados não podem possuir nenhuma coluna extra com 'cap' e 'dap' contidos em
#	seu nome, ou causará diferença no cálculo, por ex.: uma coluna chamada 'dap_medio'.
#
# b) O script contém duas funções:
#	1) fitoR(variavel_input, area_de_cada_parcela_em_m2, 'nome_do_arquivo')
#	2) similaridade(variavel_input1, variavel_input2)
#	Para utilização, carregue seus dados, execute as funções do INÍCIO até o FIM das funções
#	marcado abaixo com vários ###. Para tal selecione o texto e aperte Ctrl+R. Em seguida,
#	utilize a função desejada.
#	Na função fitoR(), o nome_do_arquivo é o nome para o arquivo .csv que contém os resultados,
#	colocar entre ' ', ou pode não funcionar direito
# 
# c) Mais algumas informações da sintaxe do uso no fim do arquivo
#
#
# Tem dúvida ou sugestão? Quer contribuir? Me envie: ricds@hotmail.com
#
#-----------------------------------------------------------------------------------------------


# # # # # # # # # # # # # # # # # # # # # # # # 
# INÍCIO DAS FUNÇÕES - CARREGAR A PARTIR DAQUI
# # # # # # # # # # # # # # # # # # # # # # # # 

# Função fitoR
 
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

#determina se existe "caps" ou "daps" e quais colunas estão
cols = grep('cap', colnames(x))
ncols = length(cols)
if (ncols>0) param="cap"

cols2 = grep('dap', colnames(x))
ncols2 = length(cols2)
if (ncols2>0) param="dap"

if (param=="dap") cols=cols2
if (param=="dap") ncols=ncols2

#calcula a area da seção transversal para cada cap/dap e faz a soma por individuo
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
 
#area basal
abta=sum(DoA)
 
#desvio da area basal entre parcelas
somag<-tapply(x$areasec, x$parc, sum)/(area/10000)
abdesv = sd(somag)

#calcula o indice de valor de importancia
VI<-(DR+DoR+FR)/3

#monta a tabela
fito=data.frame(N=N,DA=DA,DR=DR,DoA=DoA,DoR=DoR,FA=FA,FR=FR,VI=VI)
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

cat("Densidade total por área = ",round(dta,digits=2),"±",round(dtadesv,digits=2),"ind/ha", fill=TRUE)
cat("Área basal total por área = ",round(abta,digits=2),"±",round(abdesv,digits=2),"m2/ha", fill=TRUE)
cat("Riqueza = ",S,"esp.", fill=TRUE)
cat("Índice de Shannon-Wiener (H') = ",SW, fill=TRUE)
cat("Equabilidade de Pielou (J) = ",J, fill=TRUE)


if(!missing(filename)) filename = paste(filename, ".csv", sep="") else filename = 'fito.csv'
write.table(fito, file = filename, row.names = TRUE, dec=",", sep=";", quote=FALSE)
write.table(' ', file = filename, sep=";", quote=TRUE, append=TRUE, row.names=FALSE, col.names=FALSE)
cat("Densidade total por área = ",round(dta,digits=2),"±",round(dtadesv,digits=2),"ind/ha", fill=TRUE, file=filename, append=TRUE)
cat("Área basal total por área = ",round(abta,digits=2),"±",round(abdesv,digits=2),"m2/ha", fill=TRUE, file=filename, append=TRUE)
cat("Riqueza = ",S,"esp.", fill=TRUE, file=filename, append=TRUE)
cat("Índice de Shannon-Wiener (H') = ",SW, fill=TRUE, file=filename, append=TRUE)
cat("Equabilidade de Pielou (J) = ",J, fill=TRUE, file=filename, append=TRUE)


}

# Função similaridade

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

cat("A similaridade de Jaccard entre os dois sítios é de",jac,fill=TRUE)
cat("A similaridade de Sorensen entre os dois sítios é de",sor,fill=TRUE)
cat("As duas áreas compartilham ",xyspp,"esp., sendo que a primeira tem",xspp,"esp. exclusivas e a segunda tem",yspp,fill=TRUE)

}


# # # # # # # # # # # # # # # # # # # #
# FIM DAS FUNÇÕES - CARREGAR ATÉ AQUI 
# # # # # # # # # # # # # # # # # # # #










# APLICANDO AS FUNÇÕES
 

# ABRE CONJUNTO DE DADOS

# Conjunto de dados do http://labdendro.com/blog
dados=read.table(file="http://dl.dropbox.com/u/36531552/fito1.txt",header=TRUE)
dados2=read.table(file="http://dl.dropbox.com/u/36531552/fito1.txt",header=TRUE)

# Dica: Para carregar um .csv salvo diretamente do office (separadores são ponto-e-virgula e decimais são virgula)
# dados=read.table("Nome_do_Arquivo.csv", header=TRUE, sep=";", dec=",")


# Função fitoR()
# Calcula os principais descritores ecológicos, índices de diversidade, etc.
# Escreve na tela o resultado e em arquivo .csv

# Sintaxe do uso da função fitoR():
# fitoR(variavel_input, area_de_cada_parcela_em_m2, 'nome_do_arquivo')

# Exemplos:

# Sem nome do arquivo (gera um fito.csv)
fitoR(dados, 250)

# Com nome do arquivo (gera um Fito_DadosLabDendro.csv)
fitoR(dados, 250, 'Fito_DadosLabDendro')


# Função similaridade()
# Calcula a similaridade florística entre os sitios, comparando as espécies que os sitios compartilham ou não
# Índices de Jaccard e Sorensen: Baseados na matriz de presença-ausência de espécies

# Sintaxe do uso da função similaridade()
# similaridade(variavel_input1, variavel_input2)

# Exemplo:
similaridade(dados, dados2)







###

dados=read.table("Dados_exemplo.csv", header=TRUE, sep=";", dec=",")
fitoR(dados, 2500, 'Fito_Exemplo')


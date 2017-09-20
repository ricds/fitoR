# fitossociologia.R

Função para cálculo dos descritores fitossociológicos e similaridade entre sítios

Autores originais:  Alexandre Gabriel Christo e Pedro Higuchi - 25/03/2012

Modificações: Ricardo Dal'Agnol da Silva


Changelog:

Data: 06/06/2013

*Add: Riqueza, Índices Shannon-Wiener e J de Pielou

*Add: função similaridade(area1,area2)


Data: 10/06/2013

*Corrigido: equação da área basal

*Add: desvios da densidade e área basal

*Add: escreve no .csv a densidade, ab, riq, shannon e pielou

*Add: testa se existe a variavel cap ou dap, pode-se usar tanto um quanto o outro

Data: 12/06/2013

*Corrigido: equação de sorensen

*Add: possibilidade de colocar um nome para o arquivo.csv dos resultados

Data: 17/06/2013

*Add: contabiliza CAP/DAP de múltiplos fustes no DoA e área basal, max=10

Data: 27/06/2013

*Add: multiplos fustes não tem mais limite de número (antes era 10)

Observações:
a) O arquivo de entrada precisa das seguintes colunas com os nomes exatos em minusculo:
- parc (identificação das parcelas, numerico ou letra, tanto faz)
- spp (nome vulgar ou científico do indivíduo, ou utilizacao de rotulos - sp1, sp2, etc)
- dap (diâmetro a altura do peito)
ou
- cap (circunferência a altura do peito)
e quando tiver algum indivíduo com mais de um fuste/rebrota, criar colunas com nome:
- cap ou dap1, dap2, dap3, dap4, etc.
OBS: o script buscará pelas palavras-chave 'cap' e 'dap' nas colunas para localiza-los.
Portanto, os dados não podem possuir nenhuma coluna extra com 'cap' e 'dap' contidos em
seu nome, ou causará diferença no cálculo, por ex.: uma coluna chamada 'dap_medio'.

b) O script contém duas funções:
1) fitoR(variavel_input, area_de_cada_parcela_em_m2, 'nome_do_arquivo')
2) similaridade(variavel_input1, variavel_input2)
Para utilização, carregue seus dados, execute as funções do INÍCIO até o FIM das funções
marcado abaixo com vários ###. Para tal selecione o texto e aperte Ctrl+R. Em seguida,
utilize a função desejada.
Na função fitoR(), o nome_do_arquivo é o nome para o arquivo .csv que contém os resultados,
colocar entre ' ', ou pode não funcionar direito

c) Mais algumas informações da sintaxe do uso no fim do arquivo

Tem dúvida ou sugestão? Quer contribuir? Me envie: ricds@hotmail.com ou contribua aqui no github

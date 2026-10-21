## Encontro 29/09/2026
**Responsáveis:** Herqüimedes Glaudys, Murilo Durigon, Pedro Wenceslau

## Objetivo
Desenvolver um script em Python capaz de de predizer, em um cromossomo bacteriano, o ponto mais provável para a origem de replicação.

## Contexto
Cromossomos e plasmídes bacterianos possuem várias utilidades práticas na biotecnologia moderna, e para melhor utilizar essas moléculas como ferramentas é crucia conhecermos suas estruturas e mecanismos. Dentre esses conhecimentos essenciais, está a iniciação do processo de replicação natural da informação genética nesses microorganismos.
Ao longo do genoma bacteriano alguns padrões costumam indicar onde o processo de replicação do cromossomo tem início. Entre essas pistas para localização da origem de replicação(oriC) estão a presença de sítios específicos de reconhecimento da proteína iniciadora da replicação(dnaA boxes), a assimetria na distribuição de citosina e guanina e maior riqueza de adenina e timina. No entanto a análise individual desses fatores pode levar a enganos, sendo necessária uma avaliação conjunta desses indicadores.

## Dados de entrada
Arquivos no formato [fasta](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) com cromossomo bacteriano
Padrão de reconhecimento da proteína dnaA da espécie.
Origem de replicação(posição inicial), se conhecida.  

## Saída esperada
Posição da oriC no genoma bacteriano.  
Erro: distância, em pares de base(bp - base pair), entre a posição real e a estimada.  
Gráfico representando posição estimada e posição real no contexto de desvio GC acumulado.

## Ferramentas
Bibliotecas python:  
- biopython para leitura de arquivos fasta e manipulação de sequências;  
- matplotlib para construção de gráficos

[NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) para obter os arquivos fasta.  
[DoriC](https://tubic.org/doric/home) para obter a posição da oriC já conhecida e a sequência da dnaA box.  

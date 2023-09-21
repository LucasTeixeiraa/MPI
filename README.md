[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/BbTfNgtQ)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-718a45dd9cf7e7f842a935f5ebbe5719a5e09af4491e668f4dbf3b35d5cca122.svg)](https://classroom.github.com/online_ide?assignment_repo_id=11905607&assignment_repo_type=AssignmentRepo)
# Resolvendo a Equação Poisson - Bidimensional 

## Introdução

A Equação de Poisson é dada por:

∇²φ = f

Onde:

φ é uma função desconhecida (usualmente chamada de função potencial) que depende de várias variáveis, como x, y e z em três dimensões ou x e y em duas dimensões.
∇² é o operador laplaciano, que é a soma das segundas derivadas parciais em relação a todas as variáveis independentes. Em coordenadas cartesianas tridimensionais, o operador laplaciano é dado por: ∇² = ∂²/∂x² + ∂²/∂y² + ∂²/∂z².
f é uma função conhecida, também dependente das mesmas variáveis que φ, que é conhecida como fonte ou termo fonte.
Essa equação é muito utilizada em diversas áreas da física, como eletrostática, gravitação, teoria do calor, entre outras. A solução da equação de Poisson depende das condições de contorno (ou condições de fronteira) impostas ao problema, que variam conforme o contexto físico ou matemático em que ela é aplicada.

A resolução numérica da equação de Poisson por diferenças finitas é um método comum e eficiente para obter uma solução aproximada em problemas onde não é possível obter uma solução analítica de forma direta.

Para resolver a equação de Poisson usando o método de diferenças finitas, é necessário discretizar o domínio espacial em uma malha de pontos e substituir as derivadas parciais por diferenças finitas. Vamos considerar o problema em duas dimensões (x e y) para simplificar a explicação.

A equação de Poisson em duas dimensões é dada por:

∂²φ/∂x² + ∂²φ/∂y² = f(x, y)

O passo inicial é discretizar o domínio espacial. Suponha que tenhamos um domínio retangular com uma malha de pontos espaçados por Δx e Δy ao longo dos eixos x e y, respectivamente. Então, os pontos da malha serão denotados como φ(i, j), onde i é o índice para a coordenada x e j é o índice para a coordenada y. Vamos considerar um mesmo espaçamento
tanto da direção x quando na direção y, isto é, Δx  = Δy = h.

Agora, podemos aproximar as derivadas parciais por diferenças finitas. Para a segunda derivada em relação a x, usamos a seguinte aproximação:

∂²φ/∂x² ≈ (φ(i+1, j) - 2φ(i, j) + φ(i-1, j)) / h²

E, para a segunda derivada em relação a y, usamos:

∂²φ/∂y² ≈ (φ(i, j+1) - 2φ(i, j) + φ(i, j-1)) / h²

Substituindo essas aproximações na equação de Poisson, obtemos:

(φ(i+1, j) - 2φ(i, j) + φ(i-1, j)) / h² + (φ(i, j+1) - 2φ(i, j) + φ(i, j-1)) / h² = f(i, j)

Agora, podemos rearranjar a equação para isolar o ponto φ(i, j):

φ(i, j) = (φ(i+1, j) + φ(i-1, j) + φ(i, j+1) + φ(i, j-1) - h²f(i, j))/4.0

Agora, temos uma fórmula para calcular o valor aproximado de φ(i, j) com base nos valores dos pontos vizinhos e na função fonte f(i, j). Esse processo deve ser repetido para todos os pontos da malha, exceto aqueles que estão nas bordas, para os quais as condições de contorno devem ser aplicadas diretamente.

## Estudo de Caso

O código serial implementado em [serial_stencil](serial_stencil.c) resolve o seguinte problema:

∇²φ = f em [0,1] x [0,1]

com φ = 0 no contorno e com o termo fonté definido como
 - f(x,y) =  10.0 se x=0.25 e y=0.25
 - f(x,y) = -10.0 se x=0.75 e y=0.75

### Estrutura do código

No arquivo [appctx.h](appctx.h) está a definição da estrutura que armazena as informações do estudo de caso e as descrições de cada atributo.

```C++
typedef struct 
{
    int global_n;      // número global de pontos em cada direção
    int local_n;       // numero local de pontos em cada direção
    int ndof;          // número total de pontos
    int niters;        // número de iterações
    int energy;        // energia que será ingetada no sistema pelo termo fonte
    int rank;          // rank do processo
    int size;          // numero de processos
    int neighbors[4];  // vizinhos
    int coords[2];     // coordenadas do processo
    double L;          // tamanho do dominio
    double h;          // espaçamento
    double tol;

} AppCtx; 
```

A implementação númerica da resolução do problema está no arquivo [serial_stencil.c](serial_stencil.c). Seus esforços devem 
se concentrar nele. 

Por fim, o arquivo [save_array.c](save_array.c), salva o resultado no formato VTK ImageData (serial) que permite a visualização do resultado no ParaView. 

### Para compilacação e execução

O código disponibilizado é serial e compila em um sistema sem o uso do MPI. A compilação é gerenciado pelo utilitário make.
Para tanto, basta executar o seguinte comando:
```bash 
make 
```
Será gerado o executável *serial_stencil*. Para executar as opções são:
```bash 
./serial_stencil
```

É possivel alterar o tamanho do grid, tolerancia do solver e o número máximo de iterações usando
comandos por linha de argumento. Vejas os opções disponíveis fazendo:

```bash 
./serial_stencil -h
```


## O que deve ser feito?

1. Crie um arquivo mpi_stencil.c e implemente uma versão paralela.
2. Faça um estudo de desempenho, analisando *speedup* e eficiência de sua implementação.

**Observação:**
 - Descomente as linhas do Makefile referente a compilação do mpi_stencil. 
 - Altere a flag USE_MPI para 1.





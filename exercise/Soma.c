#include <stdio.h>
#include <stdlib.h>

#define VALOR 0.6f
#define NUM_ELEMENTOS 10000

float somaSequencial(float *dados, unsigned int tam){
    float soma = 0.0;
    while (tam--){
        soma += dados[tam];
    }

    return soma;
}

float somaPar(float *dados, unsigned int tam){
    if (tam == 2)
        return dados[0] + dados[1];
    if (tam == 1)
        return dados[1];

    unsigned int div = tam / 2;
    //dados+div move o ponteiro para div posições a frente
    return somaPar(dados, div) + somaPar(dados+div, tam-div); 
}

float KahanSoma(float *dados, unsigned int tam){
    float soma = 0.0;
    // compensação tenta recuperar os valores pequenos da mantissa perdidos no arredondamento
    float compensacao = 0.0;
    
    while (tam--){
        // subtrai o erro acumulado da compensação da iteração anterior
        float y = dados[tam] - compensacao;
        // realiza a soma normalmente, só que com o valor adequado
        float t = soma + y;
        // (t - soma) é a parte que realmente foi adicionada a soma
        // a parte que foi adicionada - erro acumulado, resultado é o valor da próxima compensação 
        compensacao = (t - soma) - y;
        soma = t;

    }
    // há um valor perdido na ultima iteração de compensação, ele é irrelevante?
    // podemos adequar o valor final para ser ainda mais exato?
    return soma;
}

void main(){

    float *dados = (float*) malloc(NUM_ELEMENTOS * sizeof(float));

    for (unsigned int i = 0; i < NUM_ELEMENTOS; i++){
        dados[i] = VALOR;
    }

    float soma1 = somaSequencial(dados, NUM_ELEMENTOS);
    printf("Soma sequencia: %1.15f\n", soma1);
    

    float soma2 = somaPar(dados, NUM_ELEMENTOS);
    printf("Soma par: %1.15f\n", soma2);

    printf("Diferença: %1.15f\n", soma2 - soma1);

    float soma3 = KahanSoma(dados, NUM_ELEMENTOS);
    printf("Soma Kahan: %1.15f\n", soma3);

    free(dados);

    return;
}

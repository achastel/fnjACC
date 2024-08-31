#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENACC
#include <openacc.h>
#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX_THREADS_PER_BLOCK 1024

struct p_node {
	long int d;
	int i,j;
};

struct otuName {
   char *name;
   int ordem;
};

typedef struct p_node p_node_type;



void heap_fica_seq(p_node_type vet[], int i, int qtde){
   int g, min=i;
   p_node_type aux;
   while( (2*min) <= qtde ){
      g = 2*min;
      if ((g<qtde) && (vet[g].d > vet[g+1].d))
         g = g+1;
      if (vet[min].d <= vet[g].d)
         min = qtde;
      else {
         aux = vet[min];
         vet[min] = vet[g];
         vet[g] = aux;
         min = g;
      }
   }
   return;
}


//#pragma acc heap_fica
void heap_fica(p_node_type vet[], int i, int qtde){
   int g, min=i;
   //printf("heap_fica i=%d, qtde=%d\n",i,qtde);
   p_node_type aux;
   while( (2*min) <= qtde ){
      g = 2*min;
      if ((g<qtde) && (vet[g].d > vet[g+1].d))
         g = g+1;
      if (vet[min].d <= vet[g].d)
         min = qtde;
      else {
         aux = vet[min];
         vet[min] = vet[g];
         vet[g] = aux;
         min = g;
      }
   }
   return;
}

//**************************************
// Driver Code
//**************************************
int main(int argc, char *argv[])
{
   FILE* ptr = fopen(argv[1], "r");	//arquivo contendo matriz de distancia triangular
   double f = strtod(argv[2],NULL);	//f=p/q percentual de pares a serem selecionados de uma vez
   clock_t t_start, t_end, t_1,t_2;
   double dt;

   if(ptr == NULL) {
     printf("no such file.\n");
     return 0;
   }

   int  size;						// No. OTU's
   fscanf(ptr,"%d", &size);
   printf("acc_hfica SIZE=%d\n",size);

   float s = 0.0;
   s = (size/2);
   s = (size+1)*s;
   int t = s-size;	//tamanho da matriz triangular
   
   int i,j,k;
   long int **m;			// MATRIZ DE DISTANCIAS
   long int *Q[size];	// MATRIZ Q
   long int D[size];    // SOMA DAS DISTANCIAS DE UMA OTU(linha)
   int indices[size];
   int limite=(size*f)+1;;
   p_node_type L[limite];          // LISTA DE MENORES
   struct otuName otus[size];	  // armazena nome/nós/ramos da árvore
   m=(long int**)malloc(size*sizeof(long int));
   int t_vet=t;
   int menor,maior;
   long int mij, Di, Dj;
   printf ("t_vet=%d size=%d \n",t_vet, size);
   
   char str[12];
   printf("ALOCAÇÃO LINHAS DAS MATRIZES E INICIALIZAÇÃO\n");
   for(int i=0;i<size;i++){
      D[i]=0;
      m[i]= (long int*) malloc((i+1)*sizeof(long int));
      Q[i]= (long int*) malloc(i*sizeof(long int));
      sprintf(str,"%d",i);
      for(int j=0;j<i;j++) m[i][j]=-1;
      otus[i].name = (char *) malloc(strlen(str)*sizeof(char));
      sprintf(otus[i].name,"%d",i);
      otus[i].ordem=i;
      indices[i]=0;
   }
   printf("Alocando Q\n");
   for(int i=1;i<size;i++)
	  Q[i]= (long int *) malloc(i*sizeof(long int));

   printf("Leirura da Matriz Triangular:\n");
   i=0;
   while (i<size){						// LENDO DISTANCIAS DO ARQUIVO...
     j=0;
     while(j<=i){
       fscanf(ptr, "%ld ", &m[i][j]);				// ...inserindo na matriz
       j++;
     }
     i++;
   }
   fclose(ptr);

   limite=(size*f)+1;
   p_node_type pares[limite];						//HEAP COM PARES DE VIZINHOS
   
   printf("Calculando soma das distancias das OTUS Paralelo- D[i]\n");
   #pragma acc data copy(D,m)
   #pragma acc parallel
   {
      #pragma acc loop
      for(i=0;i<size;i++){
         #pragma acc loop
         for(j=0;  j<size;j++){
            for(j=0;  j<size;j++){
               if (j!=i ) D[i]=D[i]+m[i][j];
               //if (j>i)  D[i]=D[i]+m[j][i];
            }
         }
      }
   }
   
   #pragma acc update self(D)
   //printf("Soma das Distacias na liha i:\n");
   //for(i=0;i<size;i++) printf("D[%d] = %ld == %ld \n",i,D[i],Dseq[i]);

   printf("N. OTUs=%d f=%f vet[t_vet]=%ld \n",size,f,sizeof(p_node_type)*t_vet);
   p_node_type *vet;	//vet CONTENDO VALORES DA MATRIZ...PARA ORDENAÇÃO   
   vet = (p_node_type*)malloc(t_vet*sizeof(p_node_type));
//**************************************** LAÇO PRINCIPAL *****************************
   t_start = clock();
   while(size>=2){
      printf("_______________________________________________________________\n");
      limite = (size*f);
      if (limite==0) limite=1;
      k=0;
      if (size>2){
        //printf("********************************** CALCULO DA MATRIZ Q\n");
        t_1 = clock();
        #pragma acc data copy(Q,m,D,vet,indices,k)		//copia dados para gpu
        #pragma acc parallel private(i,j)				   //inicia paralelo
        {
           #pragma acc loop
           for(i=0;i<size;i++){
              indices[i]=0;
              #pragma acc loop
              for(j=0;j<i;j++){
                 Q[i][j]=((long int)(size - 2) * m[i][j]) - (D[i] + D[j]);
		 if(i!=j){
		    vet[k].i=i;
                    vet[k].j=j;
                    vet[k].d = Q[i][j];
                    k++;
		 }
             }
           }
           t_2=clock();
           printf("k=%d MatrizQ %18.4f\n",k,((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
        }
        //#pragma acc serial
        #pragma acc update self(vet,k)
        t_1=clock();
	printf("k=%d\n",k);
	for(int h=k/2; h>0; h--)
           heap_fica(vet, h, k);
	t_2=clock();
	printf("criando heap %8.4f\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
  
        //for(int i=0;i<k;i++) printf("vet[%d]=%ld (%d,%d)\n",i,vet[i].d,vet[i].i,vet[i].j);

        // obter a lista de nohs pares[u] = True se o noh u foi escolhido; False caso contrario.
        // printf("//selecionando os menores %d do heap z\n",limite);
	t_vet=k;
        if(limite == 1) pares[0]=vet[1];
        else{
           t_1=clock();
	   j=0;
	   while((j<limite)&&(t_vet>0)){		//seleciona os menores o topo do heap   
              if(indices[vet[1].i]==0 && indices[vet[1].j]==0){
                pares[j++]=vet[1];
                indices[vet[1].i]=1; 
                indices[vet[1].j]=1; }
             /*remover elemento do heap*/
             vet[1]=vet[t_vet];
             t_vet--;
             heap_fica(vet, 1, t_vet);
           }//while()
           t_2=clock();
           printf("Selecionar pares %18.4f\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
        }
        //printf("\n PARES\n");
        //for(i=0;i<limite;i++) {
        //    printf("[%d] (%d,%d) %ld \n",i, pares[i].i, pares[i].j, pares[i].d);
        //}
        //printf("SIZE=%d\n",size);

        // printf("Juntando %d Pares de OTU's\n",limite);
        //************************************************************
        //juntando pares - paralelo mas sequencial ????
        //************************************************************
        //#pragma acc loop seq
        t_1=clock();
        for(int c=0; c<limite; c++) {
           if(size >2){
              double d_i_novo,d_j_novo=0.0;
	      double t1, t2;
              //printf("// determina menor coordenada (%d,%d)- nomenclaturas %d \n",pares[c].i,pares[c].j,c);
              if (pares[c].i > pares[c].j) { // determina menor coordenada - nomenclaturas
                 mij=m[pares[c].i][pares[c].j]; Di=D[pares[c].j]; Dj=D[pares[c].i];
                 menor=pares[c].j;
                 maior=pares[c].i; }
              else{
                 mij = m[pares[c].j][pares[c].i]; Di = D[pares[c].i]; Dj = D[pares[c].j];
                 menor=pares[c].i;
                 maior=pares[c].j;}
              
              // printf("//nomeando e determinando a ordem da nova OTU\n");
 // printf("otus[menor=%d].name=%d otus[maior=%d].name=%d\n",menor,strlen(otus[menor].name),maior,strlen(otus[maior].name));
	      char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
               
              if (otus[menor].ordem < otus[maior].ordem)
	         sprintf(aux, "(%s,%s)",otus[menor].name,otus[maior].name);
	      else{
	          sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);
	          otus[menor].ordem=otus[maior].ordem;
	      }
              // printf("sizeof(otus[menor].name)=%d sizeof(aux)=%d aux=%s\n",sizeof(otus[menor].name), sizeof(aux),aux);
	      if(sizeof(otus[menor].name)<sizeof(aux)) {
		 free(otus[menor].name);
		 otus[menor].name=(char*)malloc(sizeof(aux)*sizeof(char));
	      }
              strcpy(otus[menor].name,aux);
	      otus[maior].name = otus[size-1].name;
	      otus[maior].ordem=otus[size-1].ordem;
              // printf("Corrigindo coordenadas.\n");
              //for(int x1=0;x1<size-1;x1++)
              //   printf("[%d]=(%d,%d)\n",x1,pares[x1].i,pares[x1].j);

              //
              //#pragma acc data copy(pares)		//copia dados para gpu
              //#pragma acc loop
              for(int e=c+1; e<limite; e++){
                 if(pares[e].i==size-1)
                    if(pares[e].j > maior) {
                       pares[e].i=pares[e].j;
                       pares[e].j=maior;
	            }
	            else pares[e].i=maior;
	            if (pares[e].j==size-1)
	               pares[e].j=maior;
	      }
              //printf("CORRIGIDO\n");
              //for(int x1=0;x1<size-1;x1++)
              //   printf("[%d]=(%d,%d)\n",x1,pares[x1].i,pares[x1].j);
	           D[maior]=0;
	           //#pragma acc data copy(m,D)		//copia dados para gpu
	           //#pragma acc loop
              for(int k=0;k<size-1;k++) 
                 D[maior]=D[maior]+m[size-1][k];
                     
              // #pragma acc update self(m,D)
              //if (size >2)
              size--;
           }//if(size>2)
        }//for(int c=0; c<limite; c++)
        t_2=clock();
        //printf("demais tarefas %18.4f\n\n",((double)(t_2-t_1))/(double)CLOCKS_PER_SEC);
     }else size--;
   }//while(size>=2)
   t_end = clock();
   dt = ((double)(t_end-t_start)) / (double)CLOCKS_PER_SEC;
   printf("\nTime =%18.4f sec\n",dt);
   char aux[(strlen(otus[0].name)+strlen(otus[1].name)+4+1)];
   if(otus[0].ordem  < otus[1].ordem)
     sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
   else
     sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);
   // printf("%s\n",aux);
   return 0;
}

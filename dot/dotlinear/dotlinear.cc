#include <cmath>
#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <float.h>
#include <fstream>
#include <iomanip>
#include "giv.h"

const double D = 1.0;
const double LAMB = 5.0;
const double ED = 0.001;
const int DIM = 40;
const double N = 1000;
const double DELTA = 2.*D/N;

using namespace std;

int ind(int i, int j)
{
  int m;

  if (i >= j)
    {
      m = i*(i-1)/int(2) + j;
    }

  else
    {
      m = j*(j-1)/int(2) + i;
      //Podemos trocar a ordem dos índices porque a matriz manipulada é Hermitiana.
    }

  return m;
}

double sinal(double x)
{
  return x/sqrt(x*x);
}

double Em(double m)
  {
    return m*DELTA;
    //Nesta função, são considerados os casos de energia positivas E de negativas
  }

double V(double E, double alfa)
{
  double gamma = 0.005*(1. + alfa*cos(4.*M_PI*E/D));
    return sqrt(gamma/M_PI);
}

double delta(double ergaval, double z)
{
  //Em nossa proposta, os coeficientes de Anderson da diagonal são da forma D Lambda ^ - (m + z).
  //Isso significa que um autovalor deve ser reescrito na forma D Lambda ^ - ( m + z + gamma).
  //Podemos, portanto, descobrir em qual intervalo (m, m+1) o autovalor está e, a partir disso, determinar o gamma correspondente.
  //A convenção usual, contudo, é escolher gamma variando entre -0.5 a 0.5, por isso é melhor trabalhar com intervalos da forma (m - 0.5, m + 0.5).
  
  //Módulo de um autovalor pode ser escrito na forma D Lambda ^ -(s + z)
  //Então s é dado por:
  
  double s;

  s = -ergaval/DELTA - z;
  
  //s = m + gamma, com m ainda indeterminado e gamma também indeterminado
  // Como m é sempre um inteiro, basta descobrir qual é o objeto fracionário correspondente.
  //Por exemplo, se s = 5.4, então m = 5 e gamma = 0.4.
  //Já se s = 6.8, então m = 7 e gamma = -0.2
  //Em geral, devemos olhar se a parte fracionária de um número é maior ou menor que 0.5. 
  
  double m; //m é o inteiro mais próximo de s

  if(s > 0)
  {
  if( (s - int(s)) >= 0.5) //Se a condição neste if é verdadeira, estamos mais próximos do inteiro superior, e gamma é negativo.
  {
    m = int(s) + 1;
  }
  else //Do contrário, então estamos mais próximos do inteiro inferior.
  {
    m = int(s);
  }
  }

if(s < 0)
  {
  if( (int(s) - s) >= 0.5) //Se a condição neste if é verdadeira, estamos mais próximos do inteiro superior -desconsiderando sinal-, e gamma é negativo.
  {
    m = int(s) - 1;
  }
  else //Do contrário, então estamos mais próximos do inteiro inferior, desconsiderando o sinal.
  {
    m = int(s);
  }
  }

  //Determinado m, gamma está determinado. gamma = s - m.
  double gamma;
  gamma = s - m;
  
  //gamma e delta, na malha logarítmica, estão relacionados por:  delta = pi gamma * sinal(ergaval)
  
  double deltinha;
  deltinha = -gamma*M_PI;
  
  return deltinha;
}

int main(int argc, char* argv[])
{
  double alfa = 0.0;
  switch(argc){
  case 2:
    alfa = atof(argv[1]);
    break;

  default:
    cerr << "***No. argumentos incorreto: " << argv[0] << " alfa" << endl;
    exit(1);
  }
  const int MAXCHAR = 100;
  char nome[MAXCHAR];
  int ialf = int(1000. * alfa);
  cout << ialf << endl;
  sprintf(nome, "deltas.%04d.plt", ialf);
  ofstream delts(nome);
  
  ofstream avetd("avetd.plt");
  
  ofstream spectrad("spectrad.plt");
  

  int dim = N+1; //N níveis da banda de condução e um devido ao quantum-dot.
  
  //O termo da ponta da direita da matriz é da forma H[dim][dim]. Na forma de um vetor, isso corresponde a H[ind(dim, dim)].
  
  double *ham = new double[ind(dim, dim) + 1];
  double *aval = new double[int(dim + 1)];

  double *avalaux = new double [int(dim + 1)]; //essa matriz será necessária para calcular derivada dos autovalores com relação a z.

  double **avet = new double*[dim+1];
  for(int lin = 1; lin <=dim; lin++){
    avet[lin] = new double[dim+1];
  }
 
  int i = 0;
  int j = 0; //índices a serem usados para percorrer os vários elementos de matrizes
    
//Esse for produz a matriz nula.
    for(i = 1; i <= ind(dim, dim); i++)
      {
	ham[i] = 0.0;
      }
  
    ham[ind(1,1)] = ED;
  
    //O Hamiltoniano possui elementos na diagonal da forma Em- e Em+;  e possui elementos nas primeiras linha e coluna da forma Vm- e Vm+
    //Para manter as energias em ordem crescente, começamos com energias negativas.
    //Como escolhemos a dimensão como ímpar, há um número par disponível de posições para esses coeficientes.
    //Metade para energias negativas. Metade para energias positivas.
      
    i = 2; //i guarda a primeira linha em que aparecem os elementos negativos; nesse contexto, é a segunda.
    for(double m = 0; m <= maxm; m++)
      {
	ham[ind(i, i)] = Em(m+z, -1);
	ham[ind(i, 1)] = Vred(m+z, -1, alfa)/sqrt(2.0*D);
	i++; //para mudar de linha.
      }
  
    //Isso deu conta dos termos de energia negativa. Se tudo foi feito corretamente, i sofreu maxm + 1 incrementos (1 para cada m possível)
    //Seu valor inicial era 2. Então i está em (dim -1)/2 + 2.
    //No exemplo de dim = 5, isso é 4.
  
    for(double m = maxm; m >= 0; m--)
      {
	ham[ind(i, i)] = Em(m+z, 1);
	ham[ind(i, 1)] = Vred(m+z, 1, alfa)/sqrt(2.0*D); //isso é Vm(Em+)
	i++; //para mudar de linha.
      }
  
    //Se i estava em 2 + maxm + 1, agora recebeu outro incremento de maxm + 1, totalizando 4 + 2 maxm = 4 - 2 + dim - 1 =  dim + 1
    //Isso é esperado. Se tudo foi feito corretamente, o índice i deve ser tal que estamos na primeira linha/coluna fora da matriz
  
    // cout << i - dim << endl;
 
 
 
    //Para checar se a matriz foi corretamente escrita, podemos fazer com que o programa a escreva na tela.
/*     for(i = 1; i <= dim; i++)
       {
 	for(j = 1; j <= dim; j++)
 	  {
 	    cout << ham[ind(i, j)] << " " ;
 	  }
 	cout << endl;
       }
*/
    
    givens(dim, dim, dim, ham, aval, avet);

    for(int ell = 1; ell <= dim; ell++){
      double del = delta(aval[ell], z);
      double delTeor = atan(M_PI * pow(V(aval[ell], alfa), 2.0)/2.0 / (ED - aval[ell]));
      delts << setw(11) << aval[ell] << "\t" << del << "\t" << delTeor << "\t" << del/delTeor << endl;
      
      double dotcomponent = avet[ell][1];
      avetd << aval[ell] << "\t" << dotcomponent << endl;
      }

//Antes de encerrarmos este for, temos que guardar os autovalores em uma matriz. Quando a próxima rodada começar, a matriz aval será perdida.
//Para isso, usaremos o vetor avalaux, sobre o qual os autovalores serão guardados.
//Além disso, para o cálculo da densidade espectral de c_d, precisaremos guardar a componente do dot.

   if((z - w) < 0.5e-4)  //Esse if diz que estamos no primeiro turno.
   {
    for(int ell = 1; ell <= dim; ell++)
    {
      avalaux[ell] = aval[ell];
    }
   }
   
//Já se estivermos na segunda rodada, é hora de fazer a comparação e calcular a densidade espectral
  if((z - w) > 0.5e-4)
  {
    for(int ell = 1; ell <= dim; ell++)
    {
      double diffze = (z - w)/(aval[ell] - avalaux[ell]);
      
      //cout << diffze << endl;
      
      double spectral = avet[ell][1]*avet[ell][1]*diffze*sinal(-aval[ell]);
      
      double delTeor = atan(M_PI * pow(V(aval[ell], alfa), 2.0)/2.0 / (ED - aval[ell]));
      
      double spectraldteor = 2.0*D*sin(delTeor)*sin(delTeor)/(M_PI*M_PI*V(aval[ell], alfa)*V(aval[ell], alfa));
      spectrad << aval[ell] << "\t" << spectral << "\t" << spectraldteor << "\t" << spectral/spectraldteor << endl;
    }
  }
   
}
}
}
      


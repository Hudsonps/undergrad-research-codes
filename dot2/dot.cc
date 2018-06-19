#include <cmath>
#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <float.h>
#include <fstream>
#include <iomanip>
#include "giv.h"

const double D = 1.0; //largura da banda
const double LAMB = 10; //parâmetro de discretização logarítmica
const double ED = 0.002; //energia associada ao sítio do ponto quântico
const int DIM = 40; //dimensão da matriz que será diagonalizada

using namespace std;

int ind(int i, int j)  //Função que dá os índices de uma matriz em uma forma Hermitiana
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

double sinal(double x)   //Essa função retorna o sinal de x.
{
  return x/sqrt(x*x);
}


double valordes(double ergaval, double dim, double t)
{
  double s;
  
  if(fabs(ergaval) < 2.*t)
  {
    s =  dim*acos(ergaval/(2.*t))/M_PI;
  }
  else
  {
    s = 1000;
  }
  
  return s;
}
  

double defasagem(double ergaval, double z, double dim, double t) //dim = número de sítios CONTANDO o dot. 2t = meia largura da banda
{
  //Em nossa proposta, os coeficientes de Anderson da diagonal são da forma -2 t cos (m PI/ (N+1)), em que N é o número de níveis da banda de condução.
  //Isso significa que um autovalor deve ser reescrito na forma -2 t cos ( s PI / (N+1) ).
  //Podemos, portanto, descobrir em qual intervalo (m, m+1) o autovalor está e, a partir disso, determinar o gamma correspondente.
  // s = m + gamma + z.
  //A convenção usual, contudo, é escolhr gamma variando entre -0.5 a 0.5, por isso é melhor trabalhar com intervalos da forma (m - 0.5, m + 0.5).
  
  double s;
  
  if(fabs(ergaval) < 2.*t)
  {
    s =  dim*acos(-ergaval/(2.*t))/M_PI - z;
  }
  else
  {
    s = 0;
  }
  
  // Como m é sempre um inteiro, basta descobrir qual é o objeto fracionário correspondente.
  //Por exemplo, se s = 5.4, então m = 5 e gamma = 0.4.
  //Já se s = 6.8, então m = 7 e gamma = -0.2
  //Em geral, devemos olhar se a parte fracionária de um número é maior ou menor que 0.5. 
  
  double m;
  
  if( (s - int(s)) >= 0.5) //Se a condição neste if é verdadeira, estamos mais próximos do inteiro superior, e gamma é negativo.
  {
    m = int(s) + 1;
  }
  else //Do contrário, então estamos mais próximos do inteiro inferior.
  {
    m = int(s);
  }
  
  //Determinado m, gamma está determinado. gamma = s - m.
  double gamma;
  gamma = s - m;
  
  //gamma e delta, na malha logarítmica, estão relacionados por:  delta = pi gamma * sinal(ergaval)
  
  double deltinha;
  deltinha = gamma*M_PI; //O sinal de delta é algo que ainda tenho que verificar qual é o apropriado.
  
  return deltinha;
}

int main()
{
  
  char nome[100];
  char nome2[100];
  
  ofstream avals("avals2.plt");
  
  ofstream favalW0("avalW0.plt");
  
  ofstream avetVW2("avetVW2.plt");

  
//-------------------------------------------------------------------------------------------- 
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
  
  double Ed = 0.1;
  double V = 0.2;
  double t = 1;
  int dim = 4;

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

  //O termo da ponta da direita da matriz é da forma H[dim][dim]. Na forma de um vetor, isso corresponde a H[ind(dim, dim)].
  //A matriz do Hamiltoniano terá dim linhas e dim colunas
  
  double *ham = new double[ind(dim, dim) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval = new double[int(dim + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet = new double*[dim+1];
  for(int lin = 1; lin <= dim; lin++)
  {
    avet[lin] = new double[dim+1];
  }
  
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
  
  /*---- Vamos diagonalizar a matriz no caso em que W = 0 e V = 0 --- */
  
    for(int i = 1; i <= ind(dim, dim); i++)
    {
      ham[i] = 0.0;  //Com esse for, anulamos todos os elementos da matriz.
    }
  
    ham[ind(1,1)] = 1000; //Normalmente aqui estaria Ed, mas não queremos que ele interfira na obtenção das curvas, por isso colocamos 1000.
    ham[ind(1,2)] = 0; //Normalmene aqui estaria V. Mas tomamos V = 0.
    
    for(int i = 3; i <= (dim-1); i++)
    {
      ham[ind(i, i-1)] = t;
      ham[ind(i, i+1)] = t;
    }
    
    ham[ind(dim, dim-1)] = t;
    ham[ind(dim, dim)] = 0; //aqui seria W, mas tomamos W = 0.
    
    for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
      {
	cout << ham[ind(i, j)] << " ";
      }
      cout << endl;
    }
    cout << endl;
   
    
    givens(dim, dim, dim, ham, aval, avet);
    
    ofstream avalV0W0("avalV0W0.plt");
    ofstream avetV0W0("avetV0W0.plt");
    for(int i = 1; i<= dim; i++)
    {
      if( (dim-1) != 0)
      avalV0W0 << i << "\t" << aval[dim-i] << endl; //Queremos os autovalores em ordem decrescente. Por isso o argumento de aval é dim-1.
                                                    //O aval[dim] não é contado porque, neste, está Ed = 1000.
						    //aval[0] também não é um objeto pertinente, por isso o if.
    }
    
    int chosen = 3; //chosen denota o autovetor que escolhemos para escrever no arquivo.
    for(int i = 2; i <= dim; i++)
    {
      avetV0W0 << i-1 << "\t" << avet[dim-chosen][i] << endl; //O índice é i-1 em vez de i para desconsiderar o dot, que está no sítio 1.
							      //Assim, o sítio do dot passa a ser 0 e os outros são contados a partir de 1, como deve ser.
    }
    cout << "O autovetor escrito no arquivo avetV0W0.plt corresponde ao autovalor " << aval[dim-chosen] << endl; //escreve na tela o autovalor correspondente ao autovetor escolhido.
    
    
    //Aqui acaba o caso em em que W = 0 e V = 0
    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------    
  
/*  
  double Wall = 5.;
  
  for(int i = 1; i <= ind(dim, dim); i++)
  {
    ham[i] = 0.0;  //Com esse for, anulamos todos os elementos da matriz.
  }
  
    ham[ind(1,1)] = Ed; //Na primeira posição, fica a energia do dot.
    ham[ind(1,2)] = V;
    
    for(int i = 3; i <= (dim-1); i++)
    {
      ham[ind(i, i-1)] = t;
      ham[ind(i, i+1)] = t;
    }
    
    ham[ind(dim, dim-1)] = t;
    ham[ind(dim, dim)] = Wall;  //aqui entra o valor de W.
    
 */   
    
    /*
    for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
      {
	cout << ham[ind(i, j)] << " ";
      }
      cout << endl;
    }
    cout << endl;
    */
 
 /*   
    givens(dim, dim, dim, ham, aval, avet);
    
    ofstream avalvsdeltaW0("avalvsdeltaW00.plt");
    
    for(int i = 1; i <= dim; i++)
    {
      avalvsdeltaW0 << aval[i]  << "\t" << defasagem(aval[i], 0., dim, t) << endl; 
    }
*/
  
  /*int auxN = 1;
  for(int i = 1; i <= dim; i++)
  {
    if(fabs(avalW0[i] - Ed) > 1e-13)
    {
      double auxiliar = -2.*t*cos(auxN*M_PI/(dim));
      favalW0 << auxN << "\t" << avalW0[i] <<  "\t" << auxiliar << "\t" << avalW0[i]/auxiliar << endl;
      //favalW0 << auxN << "\t" << avetW0[40][i] << "\t" << endl;
      //cout << avalW0[40] << endl;
      auxN++;
    }
  }*/
  
  
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------

  dim = 4;
  
  //Nessa parte do programa, diagonalizaremos as matrizes no caso em que V = 0 e no caso em que V != 0.
  //Mediremos a defasagem com relação aos autovalores obtidos no caso em que V = 0, excluíndo E_d.
  //Isso porque o nosso objetivo, ao introduzir W, é gerar diversas bandas de condução.
  
  double *ham2 = new double[ind(dim, dim) + 1]; //Essa matriz guardará os elementos do Hamiltoniano quando V != 0
  double *aval2 = new double[int(dim + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano quando V != 0
  double **avet2 = new double*[dim+1]; //Guarda autovetores no caso em que V != 0
  for(int lin = 1; lin <= dim; lin++)
  {
    avet2[lin] = new double[dim+1];
  }
  
  double *hamV0 = new double [ind(dim, dim) + 1]; //Essa matriz guardará os elementos do Hamiltoniano no caso particular em que V = 0.
  double *avalV0 = new double[int(dim + 1)]; //Essa matriz guarda os autovalores no caso particular em que V = 0.
  double **avetV0 = new double*[dim+1]; //Guarda autovetores no caso em que V = 0.
  for(int lin = 1; lin <= dim; lin++)
  {
    avetV0[lin] = new double[dim+1];
  } 
  
  ofstream avalparaV0todos("avalV0todos.plt");
  
  for(double W = -20; W <= 20; W = W + 0.1)
  {
    //O if abaixo é para diferenciar o nome dos arquivos tais que W é positivo e tais que W é negativo
    if(W > 0)
    {
      sprintf(nome, "avalvsdeltaW%3.1f.plt", W);
      sprintf(nome2, "avalV0W%3.1f.plt", W);
    }
    else
    {
      sprintf(nome, "avalvsdeltamW%3.1f.plt", -W);
      sprintf(nome2, "avalV0mW%3.1f.plt", -W);
    }
      
    
    //ofstream avalvsdeltaW(nome);
    //ofstream avalparaV0(nome2);

    for(int i = 1; i <= ind(dim, dim); i++)
    {
      ham2[i] = 0.0;  //Com esse for, anulamos todos os elementos da matriz.
    }
  
    ham2[ind(1,1)] = Ed; //Na primeira posição, fica a energia do dot.
    ham2[ind(1,2)] = V;
    
    for(int i = 3; i <= (dim-1); i++)
    {
      ham2[ind(i, i-1)] = t;
      ham2[ind(i, i+1)] = t;
    }
    
    ham2[ind(dim, dim-1)] = t;
    ham2[ind(dim, dim)] = W;
    
    /*
    for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
      {
	cout << ham2[ind(i, j)] << " ";
      }
      cout << endl;
    }
    cout << endl;
    */

    
    givens(dim, dim, dim, ham2, aval2, avet2);
    
    
    for(int i = 1; i <= ind(dim, dim); i++)
    {
      hamV0[i] = 0.0;  //Com esse for, anulamos todos os elementos da matriz.
    }
  
    hamV0[ind(1,1)] = 1000; //Aqui normalmente seria Ed. Entretanto, Ed não faz parte da banda de condução, por isso não usaremos ele para medir defasagem.
    hamV0[ind(1,2)] = 0.; //Normalmente aqui seria V. Mas tomamos V = 0.
    
    for(int i = 3; i <= (dim-1); i++)
    {
      hamV0[ind(i, i-1)] = t;
      hamV0[ind(i, i+1)] = t;
    }
    
    hamV0[ind(dim, dim-1)] = t;
    hamV0[ind(dim, dim)] = W;
    
    /*
    for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
      {
	cout << hamV0[ind(i, j)] << " ";
      }
      cout << endl;
    }
    
    cout << endl;
    */

    
    givens(dim, dim, dim, hamV0, avalV0, avetV0);
    
    int chosen2 = 3;  //autovetor escolhido para ser escrito em um arquivo
    

    for(int i = 1; i <= dim; i++)
    {
     //avalparaV0 << valordes(avalV0[dim-i], dim, t) << "\t" << avalV0[dim-i] << "\t" << endl;  //Aqui, para desconsiderar o Ed, devemos começar de dim - i
     
     avalparaV0todos << valordes(avalV0[i], dim, t)/dim << "\t" << avalV0[i] << "\t" << endl;    
     //avalparaV0 << i << "\t" << avalV0[dim-i] << "\t" << endl;  //Aqui, para desconsiderar o Ed, devemos começar de dim - i
    
    /*for(int i = 1; i <= dim; i++)
    {
      avalvsdeltaW << aval2[i]  << "\t" << defasagem(aval2[i], 0., dim, t) << endl;
      avalparaV0 << avalV0[i] << "\t" << defasagem(avalV0[i], 0, dim, t) << endl;
    }
    */
    
	/*for(int j = 1; j <= dim; j++)
	{
	avetVW2 << j-1 << "\t" << avetV0[dim-1-chosen2][j] << endl;
	}
	
	cout << "O autovetor escrito no arquivo avetVW2.plt corresponde ao autovalor " << aval[dim-1-chosen] << endl;
	cout << "Isso corresponde a V = 0 e a W = " << W << endl;
	*/
    }

  }
   
}
      


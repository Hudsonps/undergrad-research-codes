# -*- coding: utf-8 -*-
from math  import *
import numpy as np
from matplotlib.pyplot as ppl
from os import system as sys
from scipy import interpolate as ipl
from scipy.stats import linregress
from numpy import *

D = 1.0 # largura da banda
LAMB = 10; #parâmetro de discretização logarítmica
ED = 0.002 #energia associada ao sítio do ponto quântico
DIM = 40 #dimensão da matriz que será diagonalizada

def ind(i, j):
  if (i >= j):
      m = i*(i-1)/2 + j
      
  else
      m = j*(j-1)/2 + i;
      #Podemos trocar a ordem dos índices porque a matriz manipulada é Hermitiana.
  return m
  


def sinal(double x):
  return x/sqrt(x*x)   #Essa função retorna o sinal de x.


def defasagem(ergaval, z, dim, t) /#dim = número de sítios CONTANDO o dot. 2t = meia largura da banda
{
  #Em nossa proposta, os coeficientes de Anderson da diagonal são da forma -2 t cos (m PI/ (N+1)), em que N é o número de níveis da banda de condução.
  #Isso significa que um autovalor deve ser reescrito na forma -2 t cos ( s PI / (N+1) ).
  #Podemos, portanto, descobrir em qual intervalo (m, m+1) o autovalor está e, a partir disso, determinar o gamma correspondente.
  #s = m + gamma + z.
  #A convenção usual, contudo, é escolhr gamma variando entre -0.5 a 0.5, por isso é melhor trabalhar com intervalos da forma (m - 0.5, m + 0.5).
  
  #s -> 
  
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
  int dim = 71;

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
      avalV0W0 << i << "\t" << aval[dim-i] << endl;
    }
    
    int chosen = 10;
    for(int i = 2; i<= dim; i++)
    {
      avetV0W0 << i-1 << "\t" << avet[dim-chosen][i] << endl;
    }
    cout << chosen << endl;
    cout << aval[dim-chosen] << endl;
    
    
    //Aqui acaba o caso em em que W = 0 e V = 0
    
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------    
  
  
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
    
    givens(dim, dim, dim, ham, aval, avet);
    
    ofstream avalvsdeltaW0("avalvsdeltaW00.plt");
    
    for(int i = 1; i <= dim; i++)
    {
      avalvsdeltaW0 << aval[i]  << "\t" << defasagem(aval[i], 0., dim, t) << endl; 
    }
  
  
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

  dim = 71;
  
  //Nessa parte do programa, diagonalizaremos as matrizes no caso em que V = 0 e no caso em que V != 0.
  //Mediremos a defasagem com relação aos autovalores obtidos no caso em que V = 0, excluíndo E_d.
  //Isso porque o nosso objetivo, ao introduzir W, é gerar diversas bandas de condução.
  
  double *ham2 = new double[ind(dim, dim) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval2 = new double[int(dim + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet2 = new double*[dim+1];
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
  
  
  for(double W = 27; W <= 27.1; W = W + 1)
  {
    sprintf(nome, "avalvsdeltaW%3.1f.plt", W);
    sprintf(nome2, "avalV0W%3.1f.plt", W);
    
    ofstream avalvsdeltaW(nome);
    ofstream avalparaV0(nome2);

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
    
    int chosen2 = 10;
    

    
    for(int i = 1; i <= dim; i++)
    {
      avalvsdeltaW << aval2[i]  << "\t" << defasagem(aval2[i], 0., dim, t) << endl;
      avalparaV0 << avalV0[i] << "\t" << defasagem(avalV0[i], 0, dim, t) << endl;
    }
    
	for(int j = 1; j <= dim; j++)
	{
	avetVW2 << j-1 << "\t" << avetV0[dim-chosen2][j] << endl;
	}
      
  }
  
  
    
      
}
      


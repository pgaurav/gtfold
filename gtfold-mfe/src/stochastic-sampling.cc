#include "stochastic-sampling.h"

#include "global.h"
#include <assert.h>
#include <stdio.h>
#include <stack>

std::stack<base_pair> g_stack;


double randdouble()
{
  return rand()/(double(RAND_MAX)+1);
}

bool feasible(int i, int j)
{
  return j-i > TURN && canPair(RNA[i],RNA[j]);
}

double U_0(int i, int j)
{
  return 1.0/u[i][j];
}

double U_ij(int i, int j)
{
  return (feasible(i,j) == true)?up[i][j]*exp(-auPenalty(i,j)/RT)/u[i][j]:0;
}

double U_hj(int i, int h, int j)
{
  return (feasible(h,j) == true)?up[h][j]*exp(-(ED5_new(h,j,h-1)+auPenalty(h,j))/RT)/u[i][j]:0;
}

double U_il(int i, int j)
{
  return ud[i][j]/u[i][j];
}

double U_s1h(int i, int h, int j)
{
  return s1[h][j]/u[i][j];
}

double U_ihlj_case1(int i, int h, int l, int j)
{
  return feasible(h,l)?up[h][l]*exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT)* (exp(-ED3_new(h,l,l+1)/RT)*u[l+2][j]) /s1[h][j]:0;
}

double U_ihlj_case2(int i, int h, int l, int j)
{
  return feasible(h,l)?up[h][l]*exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT)* (ud[l+1][j]) /s1[h][j]:0;
}

double U_ihlj_case3(int i, int h, int l, int j)
{
  return feasible(h,l)?up[h][l]*exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT)* (up[l+1][j]*exp(-(auPenalty(l+1,j)/RT))) /s1[h][j]:0;
}

double UD_il_case1(int i, int l, int j)
{
  return feasible(i,l)?(up[i][l]*exp(-1*auPenalty(i,l)/RT)*exp(-1*ED3_new(i,l,l+1)/RT)*u[l+2][j])/ud[i][j]:0;  
}

double UD_il_case2(int i, int l, int j)
{
  return feasible(i,l)?(up[i][l]*exp(-1*auPenalty(i,l)/RT)*ud[l+1][j])/ud[i][j]:0;
}

double UD_il_case3(int i, int l, int j)
{
  return feasible(i,l)?up[i][l]*exp(-1*auPenalty(i,l)/RT)*up[l+1][j]*exp(-1*auPenalty(l+1,j))/ud[i][j]:0;
}

double Q_ijH(int i, int j)
{
  return exp(-eH(i,j)/RT)/up[i][j];
}

double Q_ijS(int i, int j)
{
  return exp(-eS(i,j)/RT)*up[i+1][j-1]/up[i][j];
}

double Q_ijM(int i, int j)
{
  return upm[i][j]/up[i][j];
}

double Q_ijhlBI(int i, int j, int h, int l)
{
  return feasible(h,l)?exp(-1*eL(i,j,h,l)/RT)*up[h][l]/up[i][j]:0;
}

double UPM_ip1l_case1(int i, int l, int j)
{
  return feasible(i+1,l)?(up[i+1][l] * exp((-1)*(Ea+2*Ec+auPenalty(i+1,l))/RT) * exp(-1*(ED3_new(i+1,l,l+1)+Eb)/RT) * u1[l+2][j-1])/upm[i][j]:0;
}

double UPM_ip1l_case2(int i, int l, int j)
{
  return feasible(i+1,l)?(up[i+1][l] * exp((-1)*(Ea+2*Ec+auPenalty(i+1,l))/RT) * u1d[l+1][j-1])/upm[i][j]:0;
}

double UPM_ip2l_case1(int i, int l , int j)
{
  return feasible(i+2,l)?up[i+2][l]*exp((-1)*(Ea+2*Ec+Eb+ED3_new(i,j,i+1)+auPenalty(i+2,l))/RT) * (exp((-1)*(ED3_new(i+2,l,l+1)+Eb)/RT)*u1[l+2][j-1])/upm[i][j]:0;
}

double UPM_ip2l_case2(int i, int l , int j)
{
  return feasible(i+2,l)?(up[i+2][l]*exp((-1)*(Ea+2*Ec+Eb+ED3_new(i,j,i+1)+auPenalty(i+2,l))/RT) * u1d[l+1][j-1])/upm[i][j]:0;
}


double UPM_ijs2h(int i, int h , int j)
{
  return exp((-1)*ED3_new(i,j,i+1)/RT)* (s2[h][j] * exp((-1)*(Ea+2*Ec+(h-i-1)*Eb)/RT))/upm[i][j];
}


double UPM_ijhl_case1(int i, int h, int l, int j)
{
  return feasible(h,l)?up[h][l]*(exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT))* (exp(-ED3_new(h,l,l+1)/RT)*u1[l+2][j-1])/s2[h][j]:0;
}

double UPM_ijhl_case2(int i, int h, int l, int j)
{
  return feasible(h,l)?up[h][l]*(exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT))* (u1d[l+1][j-1])/s2[h][j]:0;
}

// u1d : case 1
double U1D_ij_il_case1(int i, int l, int j)
{
  return feasible(i,l)?(up[i][l]*exp((-1)*(Ec+auPenalty(i,l))/RT) * (f(j+1,i,l)*exp((-1)*(j-l)*Eb/RT)))/u1d[i][j]:0;
}

// u1d : case 2
double U1D_ij_il_case2(int i, int l, int j)
{
  return feasible(i,l)?(up[i][l]*exp((-1)*(Ec+auPenalty(i,l))/RT)*(exp((-1)*(ED3_new(i,l,l+1)+Eb)/RT)*u1[l+2][j]))/u1d[i][j]:0;
}

// u1d : case 3
double U1D_ij_il_case3(int i, int l, int j)
{
  return feasible(i,l)?(up[i][l]*exp((-1)*(Ec+auPenalty(i,l))/RT) * u1d[l+1][j])/u1d[i][j]:0;
}

// u1
double U1_ij(int i, int j)
{
  return u1d[i][j]/u1[i][j];
}

// u1 : sample h
double U1_ij_s3h(int i, int h, int j)
{
  return (s3[h][j] * exp((-1)*(Ec+(h-i)*Eb)/RT))/u1[i][j];
}

// u1 : sample l
double U1_j_hl_case1(int h, int l, int j)
{
  return feasible(h,l)?(up[h][l]*(exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT))* (f(j+1,h,l)*exp(-((j-l)*Eb)/RT))) /s3[h][j]:0;
}

double U1_j_hl_case2(int h, int l, int j)
{
  return feasible(h,l)?(up[h][l]*(exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT))* (exp(-(ED3_new(h,l,l+1)+Eb)/RT)*u1[l+2][j])) /s3[h][j]:0;
}

double U1_j_hl_case3(int h, int l, int j)
{
  return feasible(h,l)?(up[h][l]*(exp(-(ED5_new(h,l,h-1)+auPenalty(h,l))/RT))* (u1d[l+1][j])) /s3[h][j]:0;
}

void rnd_u(int i, int j, int* structure)
{
  double rnd = randdouble();
  double cum_prob = 0.0;

  cum_prob += U_0(i,j);
  if (rnd < cum_prob)
  {
    //set_single_stranded(i,j, structure);
    return;
  }
  
  cum_prob += U_ij(i, j);
  if (rnd < cum_prob)
  {
    base_pair bp(i,j,UP);
    g_stack.push(bp);
    return;
  }

  for (int h = i+1; h < j; ++h)
  {
    cum_prob += U_hj(i,h,j);
    if (rnd < cum_prob)
    {
      base_pair bp(h,j,UP);
      //set_single_stranded(i,h-1,structure);
      g_stack.push(bp);
      return;
    }
  }
  
  cum_prob += U_il(i,j);
  if (rnd < cum_prob)
  {
    base_pair bp1(i,j,UD);
    g_stack.push(bp1);
    return;
  }

  int h1 = -1;
  for (int h = i+1;  h < j-1; ++h)
  {
    cum_prob += U_s1h(i,h,j);
    if (rnd < cum_prob)
    {
      h1 = h;
      break;
    }
  }
  assert (h1 != -1) ;

  rnd = randdouble();
  cum_prob = 0;
  for (int l = h1+1; l < j; ++l)
  {
    cum_prob += U_ihlj_case1(i,h1,l,j);
    if (rnd < cum_prob)
    {
      //set_single_stranded(i,h1-1,structure);
      base_pair bp1(h1,l,UP);
      base_pair bp2(l+2,j,U);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return ;
    }

    cum_prob += U_ihlj_case2(i,h1,l,j);
    if (rnd < cum_prob)
    {
      //set_single_stranded(i,h1-1,structure);
      base_pair bp1(h1,l,UP);
      base_pair bp2(l+1,j,UD);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return ;
    }
    
    cum_prob += U_ihlj_case3(i,h1,l,j);
    if (rnd < cum_prob)
    {
      //set_single_stranded(i,h1-1,structure);
      base_pair bp1(h1,l,UP);
      base_pair bp2(l+1,j,UP);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return ;
    }
  }
  assert(0);
}

void rnd_ud(int i, int j, int* structure)
{
  double rnd = randdouble();
  double cum_prob = 0.0;
  for (int l = i+1; l < j ; ++l)
  {
    cum_prob += UD_il_case1(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      base_pair bp2(l+2,j,U);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }

    cum_prob += UD_il_case2(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      base_pair bp2(l+1,j,UD);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }

    cum_prob += UD_il_case3(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      base_pair bp2(l+1,j,UP);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }
  }
  assert(0);
}

void rnd_up(int i, int j, int* structure)
{
  double rnd = randdouble();
  assert(structure[i] == 0);
  assert(structure[j] == 0);
  
  set_base_pair(i,j,structure);

  rnd -= Q_ijH(i,j);
  if (rnd < 0)
  {
    //set_single_stranded(i+1,j-1,structure);
    return ;
  }

  rnd -= Q_ijS(i,j);
  if (rnd  < 0)
  {
    base_pair bp(i+1,j-1,UP);
    g_stack.push(bp);
    return ;
  }

  rnd -= Q_ijM(i,j);
  if (rnd < 0)
  {
    rnd_upm(i,j,structure);
    return;
  }

  for (int h = i+1; h < j-1; ++h)
    for (int l = h+1; l < j; ++l)
    {
      if (h == i+1 && l == j-1) continue;
      rnd -= Q_ijhlBI(i,j,h,l);
      if (rnd < 0)
      {
        base_pair bp(h,l,UP);
        g_stack.push(bp);
        return;
      }
    }
}

void rnd_u1(int i, int j, int* structure)
{
  double rnd = randdouble();
  double cum_prob = 0;
  
  cum_prob += U1_ij(i,j);
  if (rnd < cum_prob)
  {
    base_pair bp(i,j,U1D);
    g_stack.push(bp);
    return;
  }

  int h1 = -1;
  for (int h = i+1; h < j-1; ++h) 
  {
    cum_prob += U1_ij_s3h(i,h,j);
    if (rnd < cum_prob)
    {
      h1 = h;
      break;
    }
  }
 
  assert(h1 != -1);
  // sample l given h1 
  rnd = randdouble();
  cum_prob = 0;
  for (int l = h1+1; l <= j ; ++l)
  {
    cum_prob += U1_j_hl_case1(h1,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp(h1,l,UP);
      g_stack.push(bp);
      return;
    }

    cum_prob += U1_j_hl_case2(h1,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(h1,l,UP);
      base_pair bp2(l+2,j,U1);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }
    
    cum_prob += U1_j_hl_case3(h1,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(h1,l,UP);
      base_pair bp2(l+1,j,U1D);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }
  }
}

void rnd_u1d(int i, int j, int* structure)
{
  double rnd = randdouble();
  double cum_prob = 0;
  
  for (int l = i+1; l <= j; ++l)
  {
    cum_prob += U1D_ij_il_case1(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      //set_single_stranded(l+1,j,structure);
      g_stack.push(bp1);
      return;
    }

    cum_prob += U1D_ij_il_case2(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      base_pair bp2(l+2,j,U1);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }
    
    cum_prob += U1D_ij_il_case3(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i,l,UP);
      base_pair bp2(l+1,j,U1D);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return;
    }
  }
}

void rnd_upm(int i, int j, int* structure)
{
  double rnd = randdouble();
  double cum_prob = 0;

  for (int l = i+2; l < j; ++l) 
  {
    cum_prob += UPM_ip1l_case1(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i+1,l,UP);
      base_pair bp2(l+2,j-1,U1);
      g_stack.push(bp2);
      g_stack.push(bp1);
      return ;
    }

    cum_prob += UPM_ip1l_case2(i,l,j);
    if (rnd < cum_prob)
    { 
      base_pair bp1(i+1,l,UP);
      base_pair bp2(l+1,j-1,U1D);
      g_stack.push(bp2);
      g_stack.push(bp1);
      return ;
    }
  }

  for (int l = i+3; l < j; ++l) 
  {
    cum_prob += UPM_ip2l_case1(i,l,j);
    if (rnd < cum_prob)
    {
      base_pair bp1(i+2,l,UP);
      base_pair bp2(l+2,j-1,U1);
      g_stack.push(bp2);
      g_stack.push(bp1);
      return ;
    }

    cum_prob += UPM_ip2l_case2(i,l,j);
    if (rnd < cum_prob)
    { 
      base_pair bp1(i+2,l,UP);
      base_pair bp2(l+1,j-1,U1D);
      g_stack.push(bp1);
      g_stack.push(bp2);
      return ;
    }
  }

  int h1 = -1;
  for (int h = i+3; h < j-1; ++h)
  {
    cum_prob += UPM_ijs2h(i,h,j);
    if (rnd < cum_prob)
    {
      h1 = h;
      break;
    }
  }
  
  rnd = randdouble();
  cum_prob = 0;
  for (int l = h1+1; l < j; ++l)
  {
    cum_prob += UPM_ijhl_case1(i,h1,l,j);
    if (rnd < cum_prob)
    {
        base_pair bp1(h1,l,UP);
        base_pair bp2(l+2,j-1,U1);
        g_stack.push(bp1);
        g_stack.push(bp2);
        return;
    }

    cum_prob += UPM_ijhl_case2(i,h1,l,j);
    if (rnd < cum_prob)
    {
        base_pair bp1(h1,l,UP);
        base_pair bp2(l+1,j-1,U1D);
        g_stack.push(bp1);
        g_stack.push(bp2);
        return;
    }
  }
}

void rnd_structure(int* structure, int len)
{
  srand(rand());
  base_pair first(1,len,U);
  g_stack.push(first);

  while (!g_stack.empty())
  {
    base_pair bp = g_stack.top();
 //   std::cout << bp;
    g_stack.pop();
  
    if (bp.type() == U)
      rnd_u(bp.i,bp.j,structure);
    else if (bp.type() == UD)
      rnd_ud(bp.i,bp.j,structure);
    else if (bp.type() == UP)
      rnd_up(bp.i,bp.j,structure);
    else if (bp.type() == U1)
      rnd_u1(bp.i,bp.j,structure);
    else if (bp.type() == U1D)
      rnd_u1d(bp.i,bp.j,structure);
    
  }
}

void set_single_stranded(int i, int j, int* structure)
{
  for(;i<=j;++i) 
    structure[i] = 0;
}

void set_base_pair(int i, int j, int* structure)
{
    bool cond = j-i > TURN && canPair(RNA[i],RNA[j]);
    //if (cond == false) printf("UP(%d %d) = %lf \n",i,j,up[i][j]);
    assert(cond);
    //  printf("PAIRING %d %d\n",i,j);
    structure[i] = j;
    structure[j] = i;
}

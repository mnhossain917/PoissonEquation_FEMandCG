// MatVec
void MatVec_Seq(double **Mat, double *Vec, double *MV_product, int Nrow, int Ncol)
{ for(int i=0;i<Nrow; i++)
     for(int j=0;j<Ncol; j++)
        {MV_product[i]=MV_product[i]+Mat[i][j]*Vec[j];}
}

// VecVec
double VecVec_Seq(double *Vec1, double *Vec2, int N)
{  double vv_product;
    for(int i=0;i<N; i++)
      {vv_product=vv_product+Vec1[i]*Vec2[i];}
      return vv_product;
}

// ScalVec
void ScalVec_Seq(double alpha, double *VecInp, double *VecOut, int N)
{  for(int i=0;i<N; i++)
      {VecOut[i]=alpha*VecInp[i];}
}

// VecVec
void Res_bAx_Seq(double *b, double **A, double *x, double *r, int N)
{ double *qq; qq=new double[N];
  MatVec_Seq(A, x, qq, N, N);
   for(int i=0;i<N; i++)
          {r[i]=b[i]-qq[i];}
}

// VecMinusVecs
void ContVec_Pl_contVec( double cont1, double *Vec1, double cont2, double *Vec2, double *VecOut, int N)
{  for(int i=0;i<N; i++)
      {VecOut[i]=cont1*Vec1[i]+cont2*Vec2[i];}
}

// duplicate vactor
void DublicateVec(double *b1, double *b2, int N)
{  for(int i=0;i<N; i++)
          {b2[i]=b1[i];}
}


// Solving using Conjugate Gradient
void CG_LinSys(double **A, double *b, double *x_k, int N)
{ double alpha_k, beta_k;
  double *p_k; p_k=new double[N];
  double *r_k; r_k=new double[N];
  double *Ap_k; Ap_k=new double[N];
  
  Res_bAx_Seq(b, A, x_k, r_k, N);
  DublicateVec(r_k,p_k,N);

  for(int i=0;i<100;i++)
  { double dotprd_r_k=VecVec_Seq(r_k, r_k, N);
    MatVec_Seq(A, p_k, Ap_k, N, N);

    alpha_k=VecVec_Seq(r_k, r_k, N)/ VecVec_Seq(p_k, Ap_k, N);

    ContVec_Pl_contVec( 1, x_k, alpha_k, p_k, x_k, N);
    ContVec_Pl_contVec( 1, r_k, -alpha_k, Ap_k, r_k, N);

    double dotprd_r_k1=VecVec_Seq(r_k, r_k, N);

    beta_k=dotprd_r_k1/dotprd_r_k;
    
    ContVec_Pl_contVec( 1, r_k, beta_k, p_k, r_k, N);
          
  }

}



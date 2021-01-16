#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//[[Rcpp::export]]
arma::mat mvrnormArma2(int n, arma::vec mu, double sigma) {
  int ncols = mu.n_elem;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * pow(sigma, 0.5);
}

//[[Rcpp::export]]
arma::mat armaInv(arma::mat x) { return arma::inv(x); }

// [[Rcpp::export]]
mat repeach(vec x, vec y) {
  int n = y.size();
  vec myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
    ind += p;
  }
  return myvector;
}

// [[Rcpp::export]]
vec vecrev( vec x ){
  
  int N = x.n_elem;
  uvec y = linspace<uvec>(N-1, 0, N);
  vec z = x.elem(y);
  
  return z;
}

// [[Rcpp::export]]
vec vecrbeta( vec aa, vec bb ){
  
  int tempn = aa.n_elem;
  vec tempvec = zeros<vec>(tempn);
  
  for(int i=0; i<tempn; ++i){
    
    double tempval = R::rbeta( aa(i), bb(i) );
    tempvec(i) = tempval;
  }
  
  return tempvec;
}

// [[Rcpp::export]]
int rmult( vec prob ){
  
  
  vec u = randu<vec>(1); 
  double u1 = conv_to<double>::from(u);
  vec cumprob = cumsum(prob);
  
  uvec ind = find( cumprob > u1 );
  int i = ind(0)+1;
  
  return i;
}

// [[Rcpp::export]]
vec BetaFun( vec alpha, vec beta ){
  
  int k = alpha.n_elem;
  vec val = zeros<vec>(k);
  
  for(int i=0; i < k; ++i){
    
    val(i) = R::beta(alpha(i), beta(i));
  }

  return val;
}

// [[Rcpp::export]]
int rJkl(double GAM, double PPI, double LLAMB){
  
  int bin_val = 0;
  
  vec u = randu<vec>(1); 
  double u1 = conv_to<double>::from(u);
  if( u1 < PPI ){ bin_val = 1; }
  
  int kl = 1;
  if( bin_val != 1 ){
    
    kl = R::rpois(LLAMB * GAM);
  }
  
  return kl;
}

// [[Rcpp::export]]
double dJkl(int num, double GAM, double PPI, double LLAMB){
  
  double dens = 0;
  if( num == 1 ){
    
    dens = PPI;
  }
  dens += (1-PPI) * R::dpois(num, LLAMB*GAM, 0);
  
  return dens;
}

// [[Rcpp::export]]
vec seqN(int N1, int N2){
  
  vec seq = linspace<vec>(N1, N2, N2-N1+1);
  return(seq);
}

// [[Rcpp::export]]
vec seqComp(vec seq1, vec seq2){
  
  vec seq;
  int n = seq1.n_elem;
  
  for(int i=0; i < n; ++i){
    
    int x = seq1(i);
    vec ind2 = conv_to<vec>::from( find(seq2 == x) );
    
    int cont = ind2.n_elem;
    
    if( cont == 0 ){
      
      vec temp = {0};
      temp(0) = x;
      seq = join_vert(seq, temp);
    }
  }
  
  return seq;
}

// [[Rcpp::export]]
mat vecTOmat( vec target, int n1, int n2 ){
  
  mat result = zeros<mat>(n1,n2);
  
  for(int j=0; j<n2; ++j){
    
    for(int i=0; i<n1; ++i){
      
      result(i,j) = target(j*n1 + i);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
vec matTOvec( mat target ){
  
  int n1 = target.n_rows;
  int n2 = target.n_cols;
  
  vec result = zeros<vec>(n1*n2);
  
  for(int j=0; j<n2; ++j){
    
    for(int i=0; i<n1; ++i){
      
      result(j*n1 + i) = target(i,j);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
double get_logmv( vec x, vec m, mat sig ){
  
  vec err = x-m;
  mat invsig = armaInv(sig);
  
  double val = -1.0/2 * log( det((2*datum::pi) * sig) );
  val += -1.0/2 * as_scalar( err.t() * invsig * err );
  
  return val;
}

// [[Rcpp::export]]
mat getZandV(int M, double alpha, double beta, int tempKp, double ppi, double llamb, mat ZETA,
             double a_u, double b_u, double a_w, double b_w, vec Zis,
             mat ZZ, mat VV, mat Eta, vec ums, vec Sigmas, vec mm, vec wms, mat Zmu){

  int n = ZETA.n_cols;
  double gam = alpha * beta / ( tempKp + beta - 1 );
  vec  kls = zeros<vec>(tempKp);
  for(int t=0; t < tempKp; ++t){

    kls(t) = rJkl(gam, ppi, llamb);
  }

  // main
  for(int l=0; l < tempKp; ++l){

    int r = ZZ.n_cols;

    // current feature
    for(int e=0; e < r; ++e){

      vec ZZe = ZZ.col(e);
      int mlk = sum(ZZe) - ZZe(l);

      if( mlk > 0 ){

        double rp = mlk/(tempKp - mlk + beta - 1);

        vec tempZ0 = conv_to<vec>::from( ZZ.row(l) );
        tempZ0(e) = 0;
        vec tempF0 = conv_to<vec>::from( VV.row(l) ) % tempZ0;

        vec Etae = conv_to<vec>::from( Eta.row(e) );
        vec El0 = conv_to<vec>::from( ZETA.row(l) ) - conv_to<vec>::from( tempF0.t() * Eta );

        double sigl = Sigmas(l);

        double tempSIG = conv_to<double>::from( 1.0/( 1.0/sigl * ( sum(pow(Etae,2)) ) + 1.0/ums(e) ) );
        double tempMU = tempSIG * 1.0/sigl * sum(El0 % Etae)  ;

        double rl = sqrt( tempSIG / ums(e) ) * exp( pow(tempMU,2)/(tempSIG*2) ) ;

        double rr = rp * rl;
        double prob;
        if( rr > pow(10,20) ){

          prob = 1;

        }else{
          prob = rr/(rr+1);
        }

        int zzle = 0;
        double u = conv_to<double>::from(randu<vec>(1));
        if( u < prob ){ zzle = 1; }

        rowvec ZZl = ZZ.row(l);
        ZZl(e) = zzle;
        ZZ.row(l) = ZZl;

        if(zzle ==1){

          double tempu = conv_to<double>::from(randn(1));
          rowvec VVl = VV.row(l);
          VVl(e) = tempMU + sqrt(tempSIG) * tempu;
          VV.row(l) = VVl;

        }else{

          rowvec VVl = VV.row(l);
          VVl(e) = 0;
          VV.row(l) = VVl;
        }

      }
    }

    // new feature
    int kl = kls(l);
    vec ZZll = conv_to<vec>::from( ZZ.row(l) );
    vec VVll = conv_to<vec>::from( VV.row(l) );

    uvec ind1 = find( conv_to<vec>::from( (sum(ZZ, 0) - ZZ.row(l)) ) == 0 );
    uvec ind2 = find( ZZll > 0 );

    vec index_zero = conv_to<vec>::from( intersect(ind1,ind2) );
    vec nindex_zero = seqComp( seqN(0, r-1) ,index_zero );

    int current_kl = index_zero.n_elem;

    mat tempZZ = ZZ;
    mat tempVV = VV;
    mat tempEta = Eta;
    mat tempZmu = Zmu;

    if( (kl > 0) || (current_kl > 0) ){

      vec ums_plus;
      vec Vl_plus;
      vec wms_plus;
      vec mm_plus;
      mat Zmu_plus;

      if( kl > 0 ){

        ums_plus = zeros<vec>(kl);
        wms_plus = zeros<vec>(kl);

        for(int ii=0; ii < kl; ++ii){

          ums_plus(ii) = conv_to<double>::from(1.0/randg<vec>(1, distr_param(a_u, 1.0/b_u)));
          wms_plus(ii) = conv_to<double>::from(1.0/randg<vec>(1, distr_param(a_w, 1.0/b_w)));
        }

        Vl_plus = randn(kl) % sqrt(ums_plus);
        mm_plus = randn(kl);

        Zmu_plus = zeros<mat>(kl, M);

        for(int iii=0; iii < M; ++iii){

          Zmu_plus.col(iii) = randn(kl) % sqrt(wms_plus) + mm_plus;
        }

      }

      vec current_Vl_plus;
      mat current_Eta_plus;
      mat current_Zmu;

      if(current_kl > 0){

        uvec index_zero1 = conv_to<uvec>::from(index_zero);

        current_Vl_plus = VVll.elem( index_zero1 );
        current_Eta_plus = tempEta.rows( index_zero1 );
        current_Zmu = tempZmu.rows( index_zero1 );
      }

      uvec nindex_zero1 = conv_to<uvec>::from( nindex_zero );
      mat finiteZZ = tempZZ.cols(nindex_zero1);
      mat finiteVV = tempVV.cols(nindex_zero1);
      mat finiteEta = tempEta.rows(nindex_zero1);


      double arr = (R::dpois(kl, gam, 0))/(R::dpois(current_kl, gam, 0)) * dJkl(current_kl, gam, ppi, llamb)/dJkl(kl, gam, ppi, llamb);

      vec El0 = conv_to<vec>::from( ZETA.row(l) - (finiteZZ.row(l) % finiteVV.row(l)) * finiteEta );

      double old_logarl = 0;
      if(current_kl>0){

        if(current_kl==1){

          for(int i1=0; i1 < n; ++i1){

            int Lind = Zis(i1)-1;

            double current_Vl_plus1 = conv_to<double>::from(current_Vl_plus);
            double current_Zmu1 = conv_to<double>::from( current_Zmu.col(Lind) );

            double tempINVSIG = 1.0/Sigmas(l) * pow( current_Vl_plus1, 2 ) + 1 ;
            double tempSIG = 1.0/tempINVSIG;
            double tempMUdSIGi = El0(i1)/Sigmas(l)*current_Vl_plus1 + current_Zmu1;

            old_logarl += tempSIG * pow( tempMUdSIGi, 2 )/2 + log(tempSIG)/2 - pow(current_Zmu1,2)/2 ;
          }

        }else{

          for(int i2=0; i2 < n; ++i2){

            int Lind = Zis(i2)-1;
            vec current_Zmu1 = current_Zmu.col(Lind);

            mat tempINVSIG = 1.0/Sigmas(l) * current_Vl_plus * current_Vl_plus.t() + diagmat(ones<vec>(current_kl));
            mat tempSIG = armaInv(tempINVSIG);
            vec tempMUdSIGi = El0(i2)/Sigmas(l) * current_Vl_plus + current_Zmu1;

            old_logarl += conv_to<double>::from(tempMUdSIGi.t() * tempSIG * tempMUdSIGi)/2 + log( det(tempSIG) )/2 -
              conv_to<double>::from( current_Zmu1.t() * current_Zmu1 )/2;
          }
        }

      }

      double new_logarl = 0;
      if( kl>0 ){

        if(kl==1){

          for(int i3=0; i3 < n; ++i3){

            int Lind = Zis(i3)-1;

            double Vl_plus1 = conv_to<double>::from(Vl_plus);
            double Zmu_plus1 = conv_to<double>::from( Zmu_plus.col(Lind) );

            double tempINVSIG = 1.0/Sigmas(l) * pow( Vl_plus1, 2 ) + 1 ;
            double tempSIG = 1.0/tempINVSIG;
            double tempMUdSIGi = El0(i3)/Sigmas(l)*Vl_plus1 + Zmu_plus1;

            new_logarl += tempSIG * pow( tempMUdSIGi, 2 )/2 + log(tempSIG)/2 - pow(Zmu_plus1,2)/2 ;
          }

        }else{

          for(int i4=0; i4 < n; ++i4){

            int Lind = Zis(i4)-1;
            vec Zmu_plus1 = Zmu_plus.col(Lind);

            mat tempINVSIG = 1.0/Sigmas(l) * Vl_plus * Vl_plus.t() + diagmat(ones<vec>(kl));
            mat tempSIG = armaInv(tempINVSIG);
            vec tempMUdSIGi = El0(i4)/Sigmas(l) * Vl_plus + Zmu_plus1;

            new_logarl += conv_to<double>::from(tempMUdSIGi.t() * tempSIG * tempMUdSIGi)/2 + log( det(tempSIG) )/2 -
              conv_to<double>::from( Zmu_plus1.t() * Zmu_plus1 )/2;
          }
        }

      }

      double arl = exp( new_logarl - old_logarl );
      double ar = arl * arr;


      double uu = conv_to<double>::from(randu<vec>(1));
      if(uu < ar){

        if(current_kl>0){

          uvec index_zero11 = conv_to<uvec>::from( index_zero );
          vec tempZl = conv_to<vec>::from( ZZ.row(l) );
          tempZl.elem( index_zero11 ) = zeros<vec>( index_zero11.n_elem );
          
          ZZ.row(l) = conv_to<rowvec>::from(tempZl);
        }

        if(kl>0){
          
          mat Eta_plus = zeros<mat>(kl,n);
          for(int jj=0; jj < n; ++jj){
            
            int Lind = Zis(jj)-1;
            vec tempmeann = Zmu_plus.col(Lind);
            
            Eta_plus.col(jj) = randn(kl) + tempmeann;
          }
          
          mat ZZ_plus = zeros<mat>(tempKp,kl);
          ZZ_plus.row(l) = conv_to<rowvec>::from( ones<vec>(kl) );
          
          mat VV_plus = zeros<mat>(tempKp,kl);
          VV_plus.row(l) = conv_to<rowvec>::from( Vl_plus );
          
          ZZ = join_horiz(ZZ, ZZ_plus);
          VV = join_horiz(VV, VV_plus);
          Eta = join_vert(Eta, Eta_plus);
          ums = join_vert(ums, ums_plus);
          wms = join_vert(wms, wms_plus);
          mm = join_vert(mm, mm_plus);
          Zmu = join_vert(Zmu, Zmu_plus);
        }

      }

    }

  }

  // arrange matrix
  int tempr = ZZ.n_cols;
  vec zero_ind = conv_to<vec>::from( find(conv_to<vec>::from(sum(ZZ,0)) == 0) );

  if( (zero_ind.n_elem>0) && (tempr>2) ){

    if( (tempr - zero_ind.n_elem) < 2 ){

      vec zero_ind1 = zero_ind;
      int tempnn = zero_ind1.n_elem;
      zero_ind = zero_ind1.subvec(0, (tempnn-1)-1 );
    }
    uvec nzero_ind = conv_to<uvec>::from( seqComp( seqN(0, tempr-1) ,zero_ind ) );

    ZZ = ZZ.cols(nzero_ind);
    VV = VV.cols(nzero_ind);
    Eta = Eta.rows(nzero_ind);
    ums = ums.elem(nzero_ind);
    mm = mm.elem(nzero_ind);
    wms = wms.elem(nzero_ind);
    Zmu = Zmu.rows(nzero_ind);
  }

  // tempKp + tempKp + n + 1 + 1 + 1 + M
  // ZZ + VV + Eta.t() + ums + mm + wms + Zmu.t()

  mat tempZV = join_vert(ZZ, VV);
  tempZV = join_vert(tempZV, Eta.t());
  tempZV = join_vert(tempZV, ums.t());
  tempZV = join_vert(tempZV, mm.t());
  tempZV = join_vert(tempZV, wms.t());
  tempZV = join_vert(tempZV, Zmu.t());

  return tempZV;
}


// [[Rcpp::export]]
vec getUMS(int r, mat ZZ, mat VV, double a_u, double b_u){

  vec ums = zeros<vec>(r);
  
  vec temp1 = a_u + conv_to<vec>::from( sum(ZZ,0)/2 );
  vec temp2 = b_u + conv_to<vec>::from( sum(pow(VV,2),0)/2 ) ;
  
  for(int i=0; i < r; ++i){
    
    double param1 = temp1(i);
    double param2 = temp2(i);
    
    ums(i) = conv_to<double>::from(1.0/randg<vec>(1, distr_param(param1, 1.0/param2)));
  }
  
  return ums;
}

// [[Rcpp::export]]
mat getEta(int r, mat LFLM, vec Sigmas, mat ZETA, vec Zis, mat Zmu ){

  int n = ZETA.n_cols;
  mat Eta = zeros<mat>(r,n);

  for(int i=0; i < n; ++i){

    mat tempvar = armaInv( LFLM.t() * diagmat(1.0/Sigmas) * LFLM + diagmat(ones<vec>(r)) );
    vec tempmean = LFLM.t() * diagmat(1.0/Sigmas) * ZETA.col(i) + Zmu.col(Zis(i)-1);
    tempmean = tempvar * tempmean;

    Eta.col(i) = conv_to<vec>::from(rmvnorm(1, tempmean, tempvar));
  }

  return Eta;
}

// [[Rcpp::export]]
double getAlpha(mat ZZ, int containBETA, int tempKp, double beta, double a_alpha, double b_alpha, double Hp){

  double alpha;
  uvec iind = find(conv_to<vec>::from(sum(ZZ,0)) == 0);
  int K_plus = ZZ.n_cols-iind.n_elem;
  
  double param1;
  double param2;
    
  if(containBETA == 1){
    
    double Hp_beta = 0;
    for(int i=0; i < tempKp; ++i){
      
      Hp_beta += beta/( (i+1)+beta-1);
    }
    
    param1 = a_alpha + K_plus;
    param2 = b_alpha + Hp_beta;
    
    alpha = conv_to<double>::from(randg<vec>(1, distr_param(param1, 1.0/param2)));
    
  }else{
    
    param1 = a_alpha + K_plus;
    param2 = b_alpha + Hp;
    
    alpha = conv_to<double>::from(randg<vec>(1, distr_param(param1, 1.0/param2)));
  }
  
  return alpha;
}

// [[Rcpp::export]]
double getBeta(mat ZZ, int containBETA, int tempKp, double a_beta, double b_beta, double beta, double alpha){

  double next_beta;

  if(containBETA == 1){

    double new_beta = conv_to<double>::from(randg<vec>(1, distr_param(a_beta, 1.0/b_beta)));
    double Hp_new_beta = 0;
    for(int i=0; i < tempKp; ++i){

      Hp_new_beta += new_beta/( (i+1)+new_beta-1);
    }

    double Hp_beta = 0;
    for(int j=0; j < tempKp; ++j){

      Hp_beta += beta/( (j+1)+beta-1);
    }

    vec mk = conv_to<vec>::from(sum(ZZ,0));
    uvec iind = find(mk == 0);
    int K_plus = ZZ.n_cols-iind.n_elem;

    uvec niind = find(mk != 0);
    vec mk1 = mk.elem(niind);

    double log_acc_rate = K_plus * log( new_beta/beta ) + ( -alpha * (Hp_new_beta - Hp_beta) );
    log_acc_rate += sum( log(BetaFun(mk1,new_beta+tempKp-mk1)) - log(BetaFun(mk1,beta+tempKp-mk1)) );

    vec u = randu<vec>(1); 
    double u1 = conv_to<double>::from(u);
    
    if( log(u1) < log_acc_rate ){
      
      next_beta = new_beta;
      
    }else{
      
      next_beta = beta;
    }

  }else{

    next_beta = 1;
  }

  return next_beta;
}


// [[Rcpp::export]]
mat getZmu(vec Zis, mat Eta, int r, vec wms, vec mm, int M){

  mat Zmu = zeros<mat>(r,M);

  for(int i=0; i<M; ++i){

    uvec iind = find(Zis==(i+1));
    int containi = iind.n_elem;
    
    mat Etai = Eta.cols(iind);
    
    if( containi > 0 ){
      
      mat tempvar = diagmat( 1.0/(containi * ones<vec>(r) + 1.0/wms) );
      vec tempmean = tempvar * (conv_to<vec>::from(sum(Etai,1)) + diagmat(1.0/wms) * mm);
      
      Zmu.col(i) = conv_to<vec>::from(rmvnorm(1, tempmean, tempvar));
      
    }else{
      
      mat tempvar = diagmat(wms);
      vec tempmean = mm;
      
      Zmu.col(i) = conv_to<vec>::from(rmvnorm(1, tempmean, tempvar));
    }
  }

  return Zmu;
}

// [[Rcpp::export]]
vec getmm(int M, int r, vec wms, mat Zmu){
  
  vec diagelm = M * 1.0/wms + ones<vec>(r);
  mat tempvar = diagmat(1.0/diagelm);
  
  vec tempmean = conv_to<vec>::from(sum(Zmu,1));
  tempmean = conv_to<vec>::from( diagmat(1.0/wms) * tempmean ); 
  tempmean = conv_to<vec>::from( tempvar * tempmean ); 
  
  vec mm = conv_to<vec>::from(rmvnorm(1, tempmean, tempvar));
  
  return mm;
}

// [[Rcpp::export]]
vec getwms(int M, mat Zmu, double a_w, double b_w, vec mm, int r){
  
  vec wms = zeros<vec>(r);
  
  mat mm1 = repmat(mm, 1, M);
  vec temp = sum(pow(Zmu - mm1, 2), 1);
  
  double param1 = a_w + M/2;
  double param2 = b_w;
  
  for(int i=0; i < r; ++i){
    
    param2 = b_w + temp(i)/2;
    
    wms(i) = conv_to<double>::from(1.0/randg<vec>(1, distr_param(param1, 1.0/param2)));
  }
  
  return wms;
}


// [[Rcpp::export]]
vec UPDATE__til_wl( vec njs, vec rji, double aw, double bw ){
  
  int J = njs.n_elem;
  vec rji1 = rji - 1;
  
  vec ind = cumsum( conv_to<vec>::from(njs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  vec til_wl = zeros<vec>(J-1);
  for(int l=0; l<(J-1); ++l){
    
    double tempa = aw;
    double tempb = bw;
    
    vec temprjis1 = rji1.subvec( ind(l+1), ind(J)-1 );
    uvec l1ind = find(temprjis1 == l+1);
    tempa += l1ind.n_elem;
    
    for(int h=0; h<(l+1); ++h){
      
      uvec hind = find(temprjis1 == h);
      tempb += hind.n_elem;
    }
    
    til_wl(l) = R::rbeta(tempa, tempb);
  }
  
  return til_wl;
}

// [[Rcpp::export]]
mat UPDATE__wjl( int J, double til_w0, vec til_wl ){
  
  mat wjl = zeros<mat>(J,J);
  wjl(0,0) = 1;
  for(int j=1; j<J; ++j){
    
    
    vec temp_tilw = til_wl.subvec(0, j-1);
    vec temp_tilw0 = ones<vec>(1) * til_w0;
    temp_tilw = join_vert(temp_tilw0, temp_tilw);
    
    vec temp1 = vecrev(cumprod(1-vecrev(temp_tilw)));
    temp1 = temp1.subvec(1,j);
    vec temp11 = {1};
    temp1 = join_vert(temp1, temp11) % temp_tilw;
    
    if( (j+1) < J ){
      
      vec temp0 = zeros<vec>(J-j-1);
      temp1 = join_vert( temp1, temp0 );
    }
    
    wjl.row(j) = conv_to<rowvec>::from(temp1);
  }
  
  return wjl;
}

// [[Rcpp::export]]
mat UPDATE__til_pilk( vec njs, int KK, vec betak, vec alpha0j, vec rji, vec zji, double eps ){
  
  vec rji1 = rji - 1;
  vec zji1 = zji - 1;
  
  int J = njs.n_elem;
  
  vec ind = cumsum( conv_to<vec>::from(njs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  mat til_pilk = zeros<mat>(J, KK-1);
  for(int l=0; l<J; ++l){
    
    double temp_alphal = alpha0j[l];
    uvec r_lind = find( rji1 == l );
    
    if( r_lind.n_elem > 0 ){
      
      for(int k=0; k<(KK-1); ++k){
        
        double tempa = temp_alphal * betak(k);
        double tempb = temp_alphal * (1-sum(betak.subvec(0,k)));
        
        if( tempb < 0 ){
          
          tempb = -1.0 * tempb;
        }
        
        vec zji1l = zji1.elem(r_lind);
        uvec zji1lk = find(zji1l == k);
        tempa += zji1lk.n_elem;
        
        for(int kp=(k+1); kp<KK; ++kp){
          
          uvec zji1lkp = find(zji1l == kp);
          tempb += zji1lkp.n_elem;
        }
        
        double tempval = R::rbeta( tempa, tempb );
        
        if(tempval < eps){ tempval = eps; };
        if(tempval > 1-eps){ tempval = 1-eps; };
        
        til_pilk(l,k) = tempval; 
      }
      
    }else{
      
      for(int k=0; k<(KK-1); ++k){
        
        double tempa = temp_alphal * betak(k);
        double tempb = temp_alphal * (1-sum(betak.subvec(0,k)));
        
        if( tempb < 0 ){
          
          tempb = -1.0 * tempb;
        }
        
        double tempval = R::rbeta( tempa, tempb );
        
        if(tempval < eps){ tempval = eps; };
        if(tempval > 1-eps){ tempval = 1-eps; };
        
        til_pilk(l,k) = tempval; 
      }
    }
  }
  
  return til_pilk;
}

// [[Rcpp::export]]
mat UPDATE__pilk( int J, int KK, mat til_pilk ){
  
  mat pilk = zeros<mat>(J,KK);
  
  for(int j=0; j<J; ++j){
    
    vec temp = conv_to<vec>::from( til_pilk.row(j) );
    
    vec temp1 = cumprod((1-temp));
    temp1 = temp1.subvec(0,KK-3);
    vec temp11 = {1};
    
    temp1 = join_vert(temp11, temp1) % temp;
    vec temp111 = {prod(1-temp)};
    
    vec pis = join_vert( temp1, temp111 );
    pis = 1.0/sum(pis) * pis;
    
    pilk.row(j) = conv_to<rowvec>::from( pis );
  }
  
  return pilk;
}


// [[Rcpp::export]]
vec UPDATE__rji( vec njs, mat wjl, mat pilk, vec zji, mat XXdat, mat thetakstar ){
  
  int J = njs.n_elem;
  int p = thetakstar.n_cols;
  
  vec ind = cumsum( conv_to<vec>::from(njs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  vec rji = ones<vec>( sum(njs) );
  for(int j=0; j<J; ++j){
    
    int nj = njs(j);
    
    for(int i=ind(j); i<ind(j+1); ++i){
      
      int k_z1 = zji(i)-1;
      
      vec logprob = zeros<vec>(j+1);
      for(int l=0; l<(j+1); ++l){
        
        vec xi = conv_to<vec>::from(XXdat.row(i));
        vec meank = conv_to<vec>::from(thetakstar.row(k_z1));
        
        logprob(l) = log(wjl(j,l)) + log(pilk(l,k_z1)) + get_logmv( xi, meank, diagmat( ones<vec>(p) )  );
      }
      
      logprob = logprob - max(logprob);
      vec prob = exp(logprob);
      prob = prob/sum(prob);
      
      rji(i) = rmult(prob);
    }
  }
  
  return rji;
}


// [[Rcpp::export]]
vec UPDATE__zji( vec njs, mat pilk, vec rji, mat XXdat, mat thetakstar ){
  
  int J = njs.n_elem;
  int p = thetakstar.n_cols;
  int KK = thetakstar.n_rows;
  
  vec ind = cumsum( conv_to<vec>::from(njs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  vec zji = ones<vec>( sum(njs) );
  for(int j=0; j<J; ++j){
    
    int nj = njs(j);
    
    for(int i=ind(j); i<ind(j+1); ++i){
      
      vec logprob = zeros<vec>(KK);
      for(int k=0; k<KK; ++k){
        
        vec xi = conv_to<vec>::from(XXdat.row(i));
        vec meank = conv_to<vec>::from(thetakstar.row(k));
        
        int rind = rji(i)-1;
        
        logprob(k) = log(pilk(rind,k)) + get_logmv( xi, meank, diagmat( ones<vec>(p) )  );
      }
      
      logprob = logprob - max(logprob);
      vec prob = exp(logprob);
      prob = prob/sum(prob);
      
      zji(i) = rmult(prob);
    }
  }
  
  return zji;
}

// [[Rcpp::export]]
vec UPDATE__alpha0j( mat pilk, double c0, double d0 ){
  
  int J = pilk.n_rows;
  int KK = pilk.n_cols;
  
  vec alpha0j = zeros<vec>( J );
  for(int j=0; j<J; ++j){
    
    double temp = pilk(j,KK-1);
    if( temp == 0 ){ temp = 0.0001; }
    
    double param1 = c0 + KK-1;
    double param2 = d0 - log(temp);
    
    vec delta = randg<vec>(1, distr_param(param1, 1.0/param2));
    double delta1 = as_scalar(delta);
    
    alpha0j(j) = delta1;
  }
  
  return alpha0j;
}


// [[Rcpp::export]]
vec UPDATE__betak(int KK, vec njs, double eps, vec zji, double gam) { // N -> M, delta -> nu
  
  vec temp = unique(zji);
  int nstar = temp.n_elem;
  
  vec Ml = zeros<vec>(KK);
  for(int i=0; i<nstar; ++i){
    
    uvec temp2 = find(zji == temp(i));
    int tempn = temp2.n_elem;
    
    Ml(temp(i)-1) = tempn;
  }
  
  vec tempa = 1+Ml.subvec(0, KK-2);
  vec tempb = vecrev( cumsum(vecrev(Ml)) ).subvec(1, KK-1) + gam;
  
  vec Vs = vecrbeta( tempa, tempb );
  
  vec veceps = {eps};
  vec vec1eps = {1-eps};
  
  if( sum(Vs<eps)>0 ){
    
    uvec nzind = find(Vs>eps);
    double tempval = min( join_vert(Vs.elem(nzind), veceps) );
    
    uvec zind = find(Vs<eps);
    vec tempvalvec = ones( zind.n_elem ) * tempval;
    
    Vs.elem(zind) = tempvalvec;
  }
  
  if( sum(Vs>1-eps)>0 ){
    
    uvec noind = find(Vs<1-eps);
    double tempval1 = max( join_vert(Vs.elem(noind), vec1eps) );
    
    uvec oind = find(Vs>1-eps);
    vec tempvalvec1 = ones( oind.n_elem ) * tempval1;
    
    Vs.elem(oind) = tempvalvec1;
  }
  
  vec temp1 = cumprod((1-Vs));
  temp1 = temp1.subvec(0,KK-3);
  vec temp11 = {1};
  
  temp1 = join_vert(temp11, temp1) % Vs;
  vec temp111 = {prod(1-Vs)};
  
  vec pis = join_vert( temp1, temp111 );
  pis = pis/sum(pis);
  
  return pis;
}


// [[Rcpp::export]]
double UPDATE__gamma( vec betak, double gam01, double gam02 ){
  
  int KK = betak.n_elem;
  
  double temp = betak(KK-1);
  if( temp == 0 ){ temp = 0.0001; }
  
  double param1 = gam01 + KK-1;
  double param2 = gam02 - log(temp);
  
  vec delta = randg<vec>(1, distr_param(param1, 1.0/param2));
  double delta1 = as_scalar(delta);
  
  return delta1;
}

// [[Rcpp::export]]
vec getSigmas( int tempKp, double a_sigma, double b_sigma, mat Eta, mat ZZ, mat VV, mat ZETA ){

  int n = Eta.n_cols;
  vec Sigmas = zeros<vec>(tempKp);
  
  mat LFLM = ZZ % VV;
  mat tempE = ZETA - LFLM * Eta;
  vec temprate = b_sigma + sum(pow(tempE,2), 1)/2;
  
  double param1 = a_sigma + n/2;
  double param2 = b_sigma;
  
  for(int i=0; i < tempKp; ++i){
    
    param2 = temprate(i);
    
    Sigmas(i) = conv_to<double>::from(1.0/randg<vec>(1, distr_param(param1, 1.0/param2)));
  }
  
  return Sigmas;
}









// [[Rcpp::export]]
mat getTHETA0( mat nobs, mat B, vec Sigmas, int tempK, double psi, mat LFLM, mat Eta, mat Ydat ){   //////
  
  int n = Eta.n_cols;
  mat THETA = zeros<mat>(tempK,n);
  
  vec ind = cumsum( conv_to<vec>::from(nobs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  vec Ydat1 = conv_to<vec>::from( Ydat );
  mat LFLMk = LFLM.rows(0,tempK-1);
  
  for(int i=0; i < n; ++i){
    
    mat Bi = B.rows(ind(i), ind(i+1)-1);
    vec Ydat1i = Ydat1.subvec(ind(i), ind(i+1)-1);
    vec Etai = Eta.col(i);
    
    vec Sigmask = Sigmas.subvec(0,tempK-1);
    mat inv_Sigma = diagmat(1.0/Sigmask);
    
    mat tempsig = inv_Sigma + 1.0/psi * Bi.t() * Bi;
    tempsig = armaInv(tempsig);
    
    vec tempmean = 1.0/psi * Bi.t() * Ydat1i + inv_Sigma * LFLMk * Etai;
    tempmean = tempsig * tempmean;
    
    THETA.col(i) = conv_to<vec>::from( rmvnorm(1, tempmean, tempsig) );
  }
  
  return THETA;
}

// [[Rcpp::export]]
mat getTHETA( int cc, vec kkind, mat nobs, mat B, vec Sigmas, vec KKs, vec psi, mat LFLM, mat Eta, mat Ydat ){  //////

  int Ksum = sum(KKs);
  int n = Eta.n_cols;
  mat THETA = zeros<mat>(Ksum,n);

  vec Ns = conv_to<vec>::from( sum(nobs,0) );
  
  vec Nsind = cumsum(Ns);
  vec Nsind0 = {0};
  Nsind = join_vert(Nsind0, Nsind);
  
  vec KKsind = cumsum(KKs);
  vec KKsind0 = {0};
  KKsind = join_vert(KKsind0, KKsind);
  
  for(int h=0; h < cc; ++h){

    mat tempTHETAh = zeros<mat>(KKs(h),n);
    
    vec ind = cumsum(nobs.col(h));
    vec ind0 = {0};
    ind = join_vert(ind0, ind);

    vec Ydath = Ydat.col(h);
    mat LFLMh = LFLM.rows(kkind(h),kkind(h+1)-1);
    vec Sigmash = Sigmas.subvec(kkind(h),kkind(h+1)-1);
    mat inv_Sigma = diagmat(1.0/Sigmash);
    mat Bh = B.cols(KKsind(h), KKsind(h+1)-1);

    for(int i=0; i < n; ++i){

      mat Bhi = Bh.rows(ind(i), ind(i+1)-1);
      vec Ydathi = Ydath.subvec(ind(i), ind(i+1)-1);
      vec Etai = Eta.col(i);

      mat tempsig = inv_Sigma + 1.0/psi(h) * Bhi.t() * Bhi;
      tempsig = armaInv(tempsig);

      vec tempmean = 1.0/psi(h) * Bhi.t() * Ydathi + inv_Sigma * LFLMh * Etai;
      tempmean = tempsig * tempmean;

      tempTHETAh.col(i) = conv_to<vec>::from( rmvnorm(1, tempmean, tempsig) );
    }
    
    THETA.rows(kkind(h),kkind(h+1)-1) = tempTHETAh;
  }

  return THETA;
}

// [[Rcpp::export]]
vec getpsi( double a_psi, double b_psi, int cc, mat nobs, vec KKs, mat THETAcs, mat Ydat, mat B ){

  int n = nobs.n_rows;
  vec psi = zeros<vec>(cc);

  vec Ns = conv_to<vec>::from( sum(nobs,0) );

  vec Nsind = cumsum(Ns);
  vec Nsind0 = {0};
  Nsind = join_vert(Nsind0, Nsind);

  vec KKsind = cumsum(KKs);
  vec KKsind0 = {0};
  KKsind = join_vert(KKsind0, KKsind);

  for(int h=0; h < cc; ++h){

    vec Ydath = Ydat.col(h);

    int N = Ns(h);
    mat Bh = B.cols(KKsind(h), KKsind(h+1)-1);
    mat Thetah = THETAcs.rows(KKsind(h), KKsind(h+1)-1);

    vec tempMEAN = zeros<vec>(N);
    vec ind = cumsum(nobs.col(h));
    vec ind0 = {0};
    ind = join_vert(ind0, ind);

    for(int i=0; i < n; ++i){

      tempMEAN.subvec(ind(i),ind(i+1)-1) = Bh.rows(ind(i), ind(i+1)-1) * Thetah.col(i);
    }

    double param1 = a_psi + N/2;
    double param2 = b_psi + sum( pow(Ydath - tempMEAN,2) /2);

    double psih = conv_to<double>::from( 1.0/randg<vec>(1, distr_param(param1, 1.0/param2)) );
    psi(h) = psih;
  }

  return psi;
}


// [[Rcpp::export]]
double getpsi0( double a_psi, double b_psi, mat nobs, mat THETA, mat Ydat, mat B ){

  int n = nobs.n_rows;
  int N = sum( conv_to<vec>::from(nobs) );
  vec tempMEAN = zeros<vec>(N);
  vec ind = cumsum( conv_to<vec>::from(nobs) );
  vec ind0 = {0};
  ind = join_vert(ind0, ind);
  
  vec Ydat1 = conv_to<vec>::from( Ydat );
  
  for(int i=0; i < n; ++i){
    
    tempMEAN.subvec(ind(i),ind(i+1)-1) = B.rows(ind(i), ind(i+1)-1) * THETA.col(i);
  }
  
  double param1 = a_psi + N/2;
  double param2 = b_psi + sum( pow(Ydat1 - tempMEAN,2)/2);
  
  double psi = conv_to<double>::from( 1.0/randg<vec>(1, distr_param(param1, 1.0/param2)) );
  
  return psi; 
}
















// [[Rcpp::export]]
mat getgis( int tempKp, vec qind1, vec qind, vec qdims, mat Zdatstar, mat gis, mat LFLM, mat Eta, vec Sigmas, int qq ){
  
  int n = Eta.n_cols;
  mat next_gis = gis;
  
  vec tempKpqj = tempKp*ones<vec>(qind1.n_elem) + qind1;
  
  for(int ww=0; ww < qq; ++ww){
    
    int tempq = qdims(ww);
    mat tempz = Zdatstar.rows( qind(ww), qind(ww+1)-1 );
    mat tempgis = gis.rows( qind1(ww), qind1(ww+1)-1 );
    
    for(int i=0; i < n; ++i){
      
      int ind1 = conv_to<int>::from( find(tempz.col(i) == 1) );
      
      if( ind1 == (tempq-1) ){
        
        for(int j=0; j < (tempq-1); ++j){
          
          double tempmean = conv_to<double>::from( LFLM.row(tempKpqj(ww)+j) * Eta.col(i) )  ;
          double tempvar = Sigmas(tempKpqj(ww)+j);
          
          tempgis(j,i) = r_truncnorm( tempmean, sqrt(tempvar), -datum::inf, 0 );
        }
        
      }else{
        
        vec tempvec = tempgis.col(i);
        vec xx = {0};
        vec xx1 = {0};
        xx(0) = ind1;
        uvec tempind = conv_to<uvec>::from( seqComp( seqN(0, tempvec.n_elem-1), xx ) );
        vec tempvec1 = tempvec.elem(tempind);
        tempvec1 = join_vert(tempvec1, xx1);
        
        double lowbod = max(tempvec1);
        
        for(int j=0; j < (tempq-1); ++j){
          
          if( j == ind1 ){
            
            double tempmean = conv_to<double>::from( LFLM.row(tempKpqj(ww)+j) * Eta.col(i) )  ;
            double tempvar = Sigmas(tempKpqj(ww)+j);
            
            tempgis(j,i) = r_truncnorm( tempmean, sqrt(tempvar), lowbod, datum::inf );
            
          }else{
            
            double tempmean = conv_to<double>::from( LFLM.row(tempKpqj(ww)+j) * Eta.col(i) )  ;
            double tempvar = Sigmas(tempKpqj(ww)+j);
            
            tempgis(j,i) = r_truncnorm( tempmean, sqrt(tempvar), -datum::inf, tempgis(ind1,i) );
          }
        }
        
      }
    }
    
    next_gis.rows( qind1(ww), qind1(ww+1)-1) = tempgis;
  }
  
  return next_gis;
}



// [[Rcpp::export]]
mat Dahlclust(mat PSLL){

  int niter = PSLL.n_rows;
  int n = PSLL.n_cols;

  cube ppmall = zeros<cube>(n,n,niter);
  mat ppm = zeros<mat>(n,n);
  for(int i=0; i < niter; ++i){

    mat temp = diagmat( ones<vec>(n) );
    vec tempL = conv_to<vec>::from( PSLL.row(i) );

    for(int a=0; a < (n-1); ++a){

      for(int b=a+1; b < n; ++b){

        double tempval = 0;
        if( tempL(a) == tempL(b) ){

          tempval = 1;
        }

        temp(a,b) = tempval;
        temp(b,a) = tempval;
      }
    }

    ppm += temp/niter;
    ppmall.slice(i) = temp;
  }

  vec sqval = zeros<vec>(niter);
  for(int i=0; i < niter; ++i){

    sqval(i) = sum(sum(pow(ppmall.slice(i) - ppm, 2),0));
  }

  uword minind = index_min(sqval);

  vec cluster = conv_to<vec>::from( PSLL.row(minind) );
  mat result = join_horiz(cluster, ppm);

  return(result);
}


// [[Rcpp::export]]
vec getCLUST( vec burnthin, int numiter, double eps,
              vec nts, int tempKpq, int tempKp, int tempK, double ppi, double llamb, int containCOV, int containCAT,
              double Hp, vec kkind, mat nobs, mat B, vec KKs, mat Ydat, mat Xdat, mat Zdat, int qq, vec qdims, vec qind, vec qind1,
              mat ZZ, mat VV, vec ums, mat Eta, double alpha, double beta, 
              double til_w0, vec til_wl, mat wjl, mat til_pilk, mat pilk, vec rrji, vec zzji, vec alpha0j, 
              mat Zmu, vec mm, vec wms, vec til_betak, vec betak, double gam, 
              vec Sigmas, mat tempTHETA, mat ZETA, vec psi, double psi0, mat gis ){

  int nsum = sum(nts);
  int TT = nts.n_elem;
  int M = Zmu.n_cols;
  int cc = Ydat.n_cols;
  int containBETA = 1;
  
  double a_u = 1;
  double b_u = 1; 
  double a_w = 1;
  double b_w = 1;
  double a_alpha = 1; 
  double b_alpha = 1;
  double aw = 1;
  double bw = 1;
  double c0 = 1; 
  double d0 = 1;
  double gam01 = 1; 
  double gam02 = 1;
  double a_beta = 1; 
  double b_beta = 1;
  double a_sigma = 1; 
  double b_sigma = 1; 
  double a_psi = 1; 
  double b_psi = 1;
  
  mat ZZ1 = ZZ;
  mat VV1 = VV;
  vec ums1 = ums;
  mat Eta1 = Eta;
  double alpha1 = alpha;
  double beta1 = beta;
  
  vec til_wl1 = til_wl;
  mat wjl1 = wjl;
  mat til_pilk1 = til_pilk; 
  mat pilk1 = pilk;
  vec rrji1 = rrji;
  vec zzji1 = zzji;
  vec alpha0j1 = alpha0j;
  vec mm1 = mm;
  vec wms1 = wms;
  mat Zmu1 = Zmu;
  vec betak1 = betak;
  double gam1 = gam;
  
  vec Sigmas1 = Sigmas;
  mat ZETA1 = ZETA;
  mat tempTHETA1 = tempTHETA;
  vec psi1 = psi;
  double psi01 = psi0;
  mat gis1 = gis;

  mat PS_LL = zeros<mat>(nsum, numiter);
  for(int i=0; i<numiter; ++i){

    mat tempzv = getZandV( M, alpha1, beta1, tempKpq, ppi, llamb, ZETA1, a_u, b_u, a_w, b_w,
                           zzji1, ZZ1, VV1, Eta1, ums1, Sigmas1, mm1, wms1, Zmu1 );
    ZZ1 = tempzv.rows(0, tempKpq-1);
    VV1 = tempzv.rows(tempKpq, 2*tempKpq-1);
    mat Eta2 = tempzv.rows(2*tempKpq, 2*tempKpq+nsum-1);
    Eta1 = Eta2.t();
    ums1 = conv_to<vec>::from( tempzv.row(2*tempKpq+nsum) );
    mm1 = conv_to<vec>::from( tempzv.row(2*tempKpq+nsum+1) );
    wms1 = conv_to<vec>::from( tempzv.row(2*tempKpq+nsum+2) );
    mat Zmu2 = tempzv.rows(2*tempKpq+nsum+3, 2*tempKpq+nsum+M+2);
    Zmu1 = Zmu2.t();

    mat LFLM1 = ZZ1 % VV1;
    int r = Eta1.n_rows;

    ums1 = getUMS(r, ZZ1, VV1, a_u, b_u);
    Eta1 = getEta(r, LFLM1, Sigmas1, ZETA1, zzji1, Zmu1 );
    alpha1 = getAlpha(ZZ1, containBETA, tempKpq, beta1, a_alpha, b_alpha, Hp);
    beta1 = getBeta(ZZ1, containBETA, tempKpq, a_beta, b_beta, beta1, alpha1);
    
    til_wl1 = UPDATE__til_wl( nts, rrji1, aw, bw );
    wjl1 = UPDATE__wjl( TT, til_w0, til_wl1 );
    til_pilk1 = UPDATE__til_pilk( nts, M, betak1, alpha0j1, rrji1, zzji1, eps );
    pilk1 = UPDATE__pilk( TT, M, til_pilk1 );
    rrji1 = UPDATE__rji( nts, wjl1, pilk1, zzji1, Eta1.t(), Zmu1.t() );
    zzji1 = UPDATE__zji( nts, pilk1, rrji1, Eta1.t(), Zmu1.t() );
    alpha0j1 = UPDATE__alpha0j( pilk1, c0, d0 );
    betak1 = UPDATE__betak(M, nts, eps, zzji1, gam1);
    gam1 = UPDATE__gamma( betak1, gam01, gam02 );
    Zmu1 = getZmu(zzji1, Eta1, r, wms1, mm1, M);
    mm1 = getmm(M, r, wms1, Zmu1);
    wms1 = getwms(M, Zmu1, a_w, b_w, mm1, r);
    
    Sigmas1 = getSigmas( tempKpq, a_sigma, b_sigma, Eta1, ZZ1, VV1, ZETA1 );

    if( containCAT == 1 ){

      gis1 = getgis( tempKp, qind1, qind, qdims, Zdat, gis1, LFLM1, Eta1, Sigmas1, qq );
    }

    if(cc>1){

      tempTHETA1 = getTHETA( cc, kkind, nobs, B, Sigmas1, KKs, psi1, LFLM1, Eta1, Ydat );
      psi1 = getpsi( a_psi, b_psi, cc, nobs, KKs, tempTHETA1, Ydat, B );

      if( containCOV == 1 ){

        ZETA1 = join_vert( tempTHETA1, Xdat );

      }else{

        ZETA1 = tempTHETA1;
      }

    }else{

      tempTHETA1 = getTHETA0( nobs, B, Sigmas1, tempK, psi01, LFLM1, Eta1, Ydat );
      psi01 = getpsi0( a_psi, b_psi, nobs, tempTHETA1, Ydat, B );

      if( containCOV == 1 ){

        ZETA1 = join_vert( tempTHETA1, Xdat );

      }else{

        ZETA1 = tempTHETA1;
      }
    }

    if( containCAT == 1 ){

      ZETA1 = join_vert( ZETA1, gis1 );
    }

    PS_LL.col(i) = zzji1;
  }

  mat PS_LL_burnthin = PS_LL.cols( conv_to<uvec>::from(burnthin) );
  mat clustRESULT = Dahlclust( PS_LL_burnthin.t() );
  vec RESULT = clustRESULT.col(0);

  return RESULT;
}




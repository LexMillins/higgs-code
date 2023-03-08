void Poisson(TH1D* & h_obs, TH1D* & h_Sig, TH1D* & h_Bkgd){

  const double lower = 0.0;
  const double higher = 100.0;
  const double Ns_Sig_SR1 = 17.1;
  const double Ns_Sig_SR2 = 6.3;
  const unsigned int bins = 30;
  const double Ns_Bkgd_SR1 = 5.2; 
  const double Ns_Bkgd_SR2 = 0.9;


  h_obs = new TH1D("h_obs",";Number of Observed Events SR1*SR2; Counts", bins, lower, higher);
  
  h_Sig = new TH1D("h_Sig", ";;", bins, lower, higher);

  h_Bkgd = new TH1D("h_Bkgd", ";;", bins, 0.0, 100.0);

  const unsigned int nSignal = 20000;

  TRandom3* m_rand = new TRandom3(100);

  for (unsigned int n =0; n< nSignal; ++n){

    const double m1 = m_rand->Poisson(Ns_Sig_SR1);

    const double m2 = m_rand->Poisson(Ns_Sig_SR2);
    
    h_Sig->Fill(m1+m2);
    
    const double m3 = m_rand->Poisson(Ns_Bkgd_SR1);

    const double m4 = m_rand->Poisson(Ns_Bkgd_SR2);

    h_Bkgd->Fill(m3+m4);

    const double m = (m1+m2)*(m3+m4);

    h_obs->Fill(m);

  }
}

double LLH(const double mu_ggF, const double mu_VBF, const double mu_b){
  
  const double n_ggF_SR1 = 16.2;
  const double n_VBF_SR1 = 0.9;
  const double n_b_SR1 = 5.2;
  const double n_ggF_SR2 = 2.1;
  const double n_VBF_SR2 = 4.2;
  const double n_b_SR2 = 0.9;

  const double L = (mu_ggF*n_ggF_SR1 + mu_VBF*n_VBF_SR1 + mu_b*n_b_SR1 -1)*log(n_ggF_SR1 + n_VBF_SR1+ n_b_SR1) - (n_ggF_SR1 + n_VBF_SR1 + n_b_SR1) + (mu_ggF*n_ggF_SR2 + mu_VBF*n_VBF_SR2 + mu_b*n_b_SR2 -1)*log(n_ggF_SR2 + n_VBF_SR2 + n_b_SR2) - (n_ggF_SR2 + n_VBF_SR2 + n_b_SR2);
    
 
  return -1*L;
    
}

int Minimiser(const char * algoName, int printlevel){
  ROOT::Math::Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer(algoName);

  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(0.001);
  min->SetPrintLevel(printlevel);

  ROOT::Math::Functor f(&LLH, 3);

  double variable[3] = {1.,1., 1.};
  double step[3] = {0.01, 0.01, 0.01};

  min->SetFunction(f);

  min->SetVariable(0, "mu_ggF", variable[0], step[0]);
  min->SetVariable(1, "mu_VBF", variable[1], step[1]);
  min->SetVariable(2, "mu_b", variable[2], step[2]);

  min->Minimize();

  const double *vals = min->X();
  std::cout << "Minimum: f(" << vals[0] << "," << vals[1] << "," << vals[2] << ")" << min->MinValue() << std::endl;

  return 0;

}


void Fit(){
    using namespace RooFit;
 
    TH1D* h_obs = NULL;
    TH1D* h_Sig = NULL;
    TH1D* h_Bkgd = NULL;

    Poisson(h_obs, h_Sig, h_Bkgd);
    
    h_obs->Draw();
}


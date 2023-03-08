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

    const double m = m1 + m2 + m3 +m4;

    h_obs->Fill(m);

  }
}


void Fit(){
    using namespace RooFit;
 
    TH1D* h_obs = NULL;
    TH1D* h_Sig = NULL;
    TH1D* h_Bkgd = NULL;

    Poisson(h_obs, h_Sig, h_Bkgd);
    
    h_obs->Draw();
}

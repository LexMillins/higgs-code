void Poisson(TH1D* & h_obs){

  const double lower = 0.0;
  const double higher = 50.0;
  const double mean = 22.3;
  const unsigned int bins = 20;

  h_obs = new TH1D("h_obs",";;", bins, lower, higher);

  const unsigned int nSignal = 5000;

  TRandom3* m_rand = new TRandom3(100);

  for (unsigned int n =0; n< nSignal; ++n){

    const double m = m_rand->Poisson(mean);

    h_obs->Fill(m);
    
  }
}


void Fit(){
    using namespace RooFit;
    std::cout << "Hello world" << std::endl;

    TH1D* h_obs = NULL;

    Poisson(h_obs);

    h_obs->Draw();
}

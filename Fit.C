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

    RooRealVar n_meas("N", "", 20.0, h_obs->GetBinLowEdge(h_obs->GetNbinsX()+1), "");

    RooDataHist n_meas_hist("n_meas_hist", "", n_meas, h_obs);
    RooDataHist n_hist_signal("n_hist_signal", n_meas, h_Sig);
    RooDataHist n_hist_bkgd("n_hist_bkgd", n_meas, h_Bkgd);

    RooHistPdf pdf_signal("pdf_signal", "", n_meas, n_hist_signal);
    RooHistPdf pdf_bkgd("pdf_bkgd", "", n_meas, n_hist_bkgd);

    RooRealVar mu_signal("mu_signal", "", 1.0, -1000.0, 1000.0);
    RooRealVar mu_bkgd("mu_bkgd", 1.0, -1000.0, 1000.0);

    mu_bkgd.setConstant(kTRUE);
    mu_bkgd.setVal(0);

    RooRealVar N_Signal_pred("n_Signal_pred", "", h_Sig->Integral());
    RooRealVar N_Bkgd_pred("n_Bkgd_pred", "", h_Bkgd->Integral());

    RooFormulaVar N_Signal("N_Signal", "mu_signal*n_Signal_pred", RooArgSet(mu_signal, n_Signal_pred));
    RooFormulaVar N_Bkgd("N_Bkgd", "mu_bkgd*n_Bkgd_pred", RooArgSet(mu_bkgd, n_Bkgd_pred));

    RooExtendPdf epdf_Signal("epdf_Signal", "", pdf_signal, N_Signal);
    RooExtendPdf epdf_Bkgd("epdf_Bkgd", "", pdf_bkgd, N_Bkgd);

    RooAddPdf pdf_Total("pdf_Total", "", RooArgList(epdf_Signal, epdf_Bkgd));

    RooFitResult* fit_result = pdf_Total.fitTo(n_meas_hist, Save());

    TCanvas* c = new  TCanvas("c", "", 800, 800);

    RooPlot* frame = n_meas.frame(Title(""));

    n_meas_hist.plotOn(frame, Name("Measured"));
    
    pdf_Total.plotOn(frame, Normalization(1.0, RooAbsReal::RelativeExpected), Name("Total"));

    pdf_Total.plotOn(frame, Components(epdf_Signal), LineColor(kRed), Normalization(1.0, RooAbsReal::RelativeExpected), Name("Signal"));

    pdf_Total.plotOn(frame, Components(epdf_Bkgd), LineColor(kBlue), Normalization(1.0, RooAbsReal::RelativeExpected), Name("Background"));

    frame->Draw();

    const double nll = fit_result->minNll();

    std::cout << "nll = " << nll << std::endl;
    
    fit_result->Print();

}

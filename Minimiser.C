#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

using namespace RooFit;

void Minimiser(){
  
  RooRealVar mu_ggF("mu_ggF", "ggF signal strength", 1, -1.0, 2.0);
  RooRealVar mu_VBF("mu_VBF", "VBF signal strength", 1, -1.0, 2.0);

  // Expected number of events
  RooRealVar n_ggF_SR1("n_ggF_SR1", "", 16.2);
  RooRealVar n_VBF_SR1("n_VBF_SR1", "", 0.9);
  RooRealVar n_b_SR1("n_b_SR1", "", 5.2);
  RooRealVar n_ggF_SR2("n_ggF_SR2", "", 2.1);
  RooRealVar n_VBF_SR2("n_VBF_SR2", "", 4.2);
  RooRealVar n_b_SR2("n_b_SR2", "", 0.9);

  // Data 
  RooRealVar N_SR1_obs("N_SR1_obs", "", 22.3);
  RooRealVar N_SR2_obs("N_SR2_obs", "", 7.2);

  // Model
  RooFormulaVar N_SR1("N_SR1", "mu_ggF*n_ggF_SR1 + mu_VBF*n_VBF_SR1 + n_b_SR1", RooArgSet(mu_ggF, n_ggF_SR1, mu_VBF, n_VBF_SR1, n_b_SR1));
  RooFormulaVar N_SR2("N_SR2", "mu_ggF*n_ggF_SR2 + mu_VBF*n_VBF_SR2 + n_b_SR2", RooArgSet(mu_ggF, n_ggF_SR2, mu_VBF, n_VBF_SR2, n_b_SR2));  

  // Build PDFs
  RooRealVar N("N", "", 0, 10000);
  RooPoisson p1("p1", "", N, N_SR1);
  RooPoisson p2("p2", "", N, N_SR2);

  //RooAddPdf model_SR1("model_SR1", "", RooArgSet(p1));
  //RooAddPdf model_SR2("model_SR2", "", RooArgSet(p2));

  RooDataSet *data_SR1 = p1.generate(RooArgSet(N), 1000);
  RooDataSet *data_SR2 = p2.generate(RooArgSet(N), 1000);

  RooCategory sample("sample", "sample");
  sample.defineType("sr1");
  sample.defineType("sr2");

  RooDataSet combData("combData", "", N, Index(sample), Import("sr1", *data_SR1), Import("sr2", *data_SR2));

  RooSimultaneous simPdf("simPdf", "", sample);

  simPdf.addPdf(p1, "sr1");
  simPdf.addPdf(p2, "sr2");

   
  // Minimise likelihood

  RooAbsReal * nll = simPdf.createNLL(combData);

  RooMinimizer m(*nll);

  m.setVerbose(true);

  //mu_VBF.setConstant(kTRUE);
  
  m.migrad();

  //model.getParameters(N)->Print("s");

  m.setVerbose(false);
  
  m.hesse();

  m.minos();

  /*
  RooPlot *frame = mu_ggF.frame(Bins(10), Range(0.0, 2.0));
  nll->plotOn(frame, ShiftToZero());

  new TCanvas("", "", 600, 600);

  frame->Draw();
  */
  

  RooPlot *frame = m.contour(mu_ggF, mu_VBF, 1, 2, 3);
  
  new TCanvas("", "", 600, 600);

  //frame->SetMinimum(-2.0);
  //frame->SetMaximum(2.0);
  frame->SetTitle("Signal strength contour plot");

  auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  frame->Draw();
 
}

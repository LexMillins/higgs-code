#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

using namespace RooFit;

void EW_LLH(){

  //Variables

  RooRealVar mH("mH", "Higgs mass", 125, 20, 150);


  // ----------------------- W mass -----------------------------
  //Parameters
  RooRealVar M_W_ini("M_W_ini", "", 80.3799);
  RooRealVar c1("c1", "", 0.05429);
  RooRealVar c2("c2", "", 0.008939);
  RooRealVar c3("c3", "", 0.000089);
  RooRealVar c4("c4", "", 0.000161);
  RooRealVar c5("c5", "", 1.070);
  RooRealVar c6("c6", "", 0.5256);
  RooRealVar c7("c7", "", 0.0678);
  RooRealVar c8("c8", "", 0.00179);
  RooRealVar c9("c9", "", 0.0000659);
  RooRealVar c10("c10", "", 0.0737);
  RooRealVar c11("c11", "", 114.9);

  RooRealVar mt("mt", "", 172.4);
  RooRealVar mt_err("mt_err", "", 0.1724);
  RooRealVar mZ("mZ", "", 91.1875);
  RooRealVar mW_meas("mW_meas", "", 80.399);
  RooRealVar mW_err("mW_err", "", 0.80399);
  RooRealVar alpha_s("alpha_s", "", 0.1176);
  RooRealVar delta_alpha("delta_alpha", "", 0.059107);

  RooFormulaVar dt("dt", "(mt/174.3)**2 -1", RooArgSet(mt));
  RooFormulaVar dZ("dZ", "mZ/91.1875 -1", RooArgSet(mZ));
  RooFormulaVar dalpha("dalpha", "delta_alpha/0.05907 -1", RooArgSet(delta_alpha));
  RooFormulaVar dalpha_s("dalpha_s", "alpha_s/0.119 -1", RooArgSet(alpha_s));

  
  RooFormulaVar dH("dH", "log(mH/100)", RooArgSet(mH));
  RooFormulaVar dh("dh", "(mH/100)**2", RooArgSet(mH));

  RooFormulaVar H_terms("H_terms", "-c1*dH -c2*(dH**2) + c3*(dH**4)", RooArgList(c1, c2, dH, c3));
  RooFormulaVar term5("term5", "c4*(dh-1)", RooArgList(c4, dh));
  RooFormulaVar term6("term6", "-c5*dalpha", RooArgList(c5, dalpha));
  RooFormulaVar t_terms("t_terms", "c6*dt - c7*dt**2", RooArgList(c6, dt, c7));
  RooFormulaVar mix_terms("mix_terms", "-c8*dH*dt + c9*dh*dt", RooArgList(c8, dH, dt, c9, dh));
  RooFormulaVar final_term("final_term", "-c10*dalpha_s + c11*dZ", RooArgList(c10, dalpha_s, c11, dZ));

  //Build gaussian from formula vars
 
  RooGenericPdf mW("mW", "exp(-((M_W_ini + H_terms + term5 + term6 + t_terms + mix_terms + final_term - mW_meas)**2)/(2*((mW_err)**2)))", RooArgList(M_W_ini, H_terms, term5, term6, t_terms, mix_terms, final_term, mW_meas, mW_err));

  // Build PDFs

  RooRealVar M("M", "", 0, 1000);
  RooGaussian gt("gt", "", M, mt, mt_err);

  RooProdPdf pdf("pdf", "", mW, gt);

  RooDataSet *data = pdf.generate(M, 10000);


  // Minimise likelihood

  RooAbsReal * nll = mW.createNLL(*data);

  // ----------------- Weak Mixing Angle ----------------



  RooMinimizer m(*nll);

  m.setVerbose(true);
  
  m.migrad();

  //model.getParameters(mH)->Print("s");

  m.setVerbose(false);
  
  m.hesse();

  RooPlot *frame = mH.frame(Bins(10), Range(20.0, 150.0));
  nll->plotOn(frame, ShiftToZero());

  new TCanvas("", "", 600, 600);

  frame->Draw();
  
  /*
  RooPlot *frame = m.contour(mH, M, 1, 2, 3);

  new TCanvas("", "", 600, 600);
  //gPad->SetLeftMragin(0.15);
  frame->Draw();
  */
  
}

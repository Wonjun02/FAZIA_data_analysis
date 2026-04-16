#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

void draw_cut()
{
  gStyle->SetOptStat(0);
  int Energy = 60; // 분석한 에너지 설정

  // 1. 파일 열기
  TString inFileName = Form("d_cut_histCorr_%dMeV.root", Energy);
  TFile *file = new TFile(inFileName, "READ");
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open %s\n", inFileName.Data());
    return;
  }

  // 2. draw_Si2_CsI_all.C와 동일한 padMap 설정 (row 1 ~ row 4 배치)
  int padMap[16][2] = {                                                                                                                                                                                                                               
      {1,1}, {1,2}, {2,1}, {2,2},   // row 1
      {1,4}, {1,3}, {2,4}, {2,3},   // row 2
      {4,1}, {4,2}, {3,1}, {3,2},   // row 3
      {4,4}, {4,3}, {3,4}, {3,3}    // row 4
  };

  // 3. 4개의 결과에 대한 캔버스 생성
  TCanvas* cCorr      = new TCanvas("cCorr", "Correlation", 1400, 1200);
  TCanvas* cCorrCut   = new TCanvas("cCorrCut", "Cut_Correlation", 1400, 1200);
  TCanvas* cProjX     = new TCanvas("cProjX", "Cut_x_Projection", 1400, 1200);
  TCanvas* cProjY     = new TCanvas("cProjY", "Cut_y_Projection", 1400, 1200);

  TCanvas* canvases[4] = {cCorr, cCorrCut, cProjX, cProjY};
  for(int i=0; i<4; i++) canvases[i]->Divide(4, 4);

  // 4. 4x4 배열 루프 돌며 그리기
  for (int pad = 1; pad <= 16; pad++)
  {
    int q = padMap[pad-1][0];
    int t = padMap[pad-1][1];
    TString suffix = Form("_q%d_t%d", q, t);

    // 각 캔버스에 히스토그램 로드 및 그리기
    // 1) Original Correlation
    cCorr->cd(pad);
    gPad->SetLogz();
    TH2D* h1 = (TH2D*)file->Get("histCorr" + suffix);
    if(h1) h1->Draw("colz");

    // 2) Cut Correlation
    cCorrCut->cd(pad);
    gPad->SetLogz();
    TH2D* h2 = (TH2D*)file->Get("histCorr_cut" + suffix);
    if(h2) h2->Draw("colz");

    // 3) X Projection (Cut)
    cProjX->cd(pad);
    TH1D* h3 = (TH1D*)file->Get("histProj_x_cut" + suffix);
    if(h3) h3->Draw();

    // 4) Y Projection (Cut)
    cProjY->cd(pad);
    TH1D* h4 = (TH1D*)file->Get("histProj_y_cut" + suffix);
    if(h4) h4->Draw();
  }

  // 5. 저장 경로 설정 및 폴더 생성
  TString path = Form("./d_cut_results_%dMeV/", Energy);
  gSystem->Exec("mkdir -p " + path);

  // 6. PDF 저장
  cCorr->SaveAs(path + "histCorr_all.pdf");
  cCorrCut->SaveAs(path + "histCorr_cut_all.pdf");
  cProjX->SaveAs(path + "histProj_x_cut_all.pdf");
  cProjY->SaveAs(path + "histProj_y_cut_all.pdf");

  std::cout << "All cut result PDFs have been saved to: " << path << std::endl;
}

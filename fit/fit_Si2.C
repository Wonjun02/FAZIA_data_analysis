#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <iostream>

void fit_Si2()
{
  // ================= [1. 대상 선택 및 설정] =================
  int Energy = 70; // 분석한 에너지 설정
  int target_q = 4; // 원하는 Quadrant 번호
  int target_t = 1; // 원하는 Telescope 번호

  //for 60 MeV
  //Q1T2
  //double mean = 143.3, sigma = 5.947;
  //double mean = 127.6, sigma = 6.104;
  //Q4T4
  //double mean = 176.2, sigma = 6.487;
  //double mean = 158.8, sigma = 5.63;
  //for 70 MeV
  //Q1T2
  //double mean = 172.5, sigma = 9.965;
  //double mean = 155.0, sigma = 8.012;
  //Q4T4
  //double mean = 208.3, sigma = 10.00;
  //double mean = 190.0, sigma = 7.078;
  
  double mean = 150.3, sigma = 12.;
  double fit_min = 2.75, fit_max = 6.; 
  
  double CsI_min = mean - 0.5*sigma, CsI_max = mean + 0.5*sigma;
  //double CsI_min = mean - sigma, CsI_max = mean + sigma;
  //double CsI_min = 30., CsI_max = mean/2.;

  // Si2 피크 피팅을 수행할 Y축 범위 설정
  // !! 실제 Y축(Si2) 데이터의 피크 위치에 맞게 적절히 조정하세요 !!
  //double fit_min = 0., fit_max = 5.; 
  //double fit_min2 = 1.287-1.517, fit_max2 = 1.287+1.517; 
  //double fit_min2 = 2.199-1.423, fit_max2 = 2.199+1.423; 

  // ================= [2. 파일 및 2D 히스토그램 로드] =================
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("e");
  gROOT->ForceStyle();

  TString inFileName = Form("d_cut_histCorr_%dMeV.root", Energy);
  TFile *file = new TFile(inFileName, "READ");
  if (!file || file->IsZombie()) {
    printf("Error: Cannot open %s\n", inFileName.Data());
    return;
  }

  TString suffix = Form("_q%d_t%d", target_q, target_t);
  //TH2D* hCorr = (TH2D*)file->Get("histCorr" + suffix);
  TH2D* hCorr = (TH2D*)file->Get("histCorr_cut" + suffix);
  
  if(!hCorr) {
    printf("Error: Cannot find histCorr%s\n", suffix.Data());
    return;
  }

  // ================= [3. Si2 Projection 생성] =================
  // 지정된 CsI 구간의 실제 Bin을 찾아 Y축 투영 (draw_cal.C 방식)
  int bin_min = hCorr->GetXaxis()->FindBin(CsI_min);
  int bin_max = hCorr->GetXaxis()->FindBin(CsI_max);
  
  TH1D* hProjY = hCorr->ProjectionY(Form("hProjY_cut%s", suffix.Data()), bin_min, bin_max);
  hProjY->SetTitle(Form("Si2 Fit: Qua %d Tel %d (CsI %.1f-%.1f)", target_q, target_t, CsI_min, CsI_max));

  // ================= [4. 초기값 자동 추출 및 피팅] =================
  TCanvas* cFit = new TCanvas("cFit", Form("Fit_Si2_q%d_t%d", target_q, target_t), 800, 600);

  // 탐색 구간(fit_min ~ fit_max) 내에서 최대값(피크) 찾기
  hProjY->GetXaxis()->SetRangeUser(fit_min, fit_max);
  double size0 = hProjY->GetMaximum();
  double mean0 = hProjY->GetXaxis()->GetBinCenter(hProjY->GetMaximumBin());
  hProjY->GetXaxis()->SetRange(0, 0); // 시각화를 위해 X축 표시 범위 원상복구

  // 단일 가우시안 피팅 함수 정의
  TF1 *fitFunc = new TF1("fitFunc", "gaus(0)", fit_min, fit_max);
  //TF1 *fitFunc = new TF1("fitFunc", "pol3(0)", fit_min, fit_max);
  //TF1 *fitFunc = new TF1("fitFunc", "gaus(0)+pol2(3)", fit_min, fit_max);
  //TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3)", fit_min, fit_max);
  
  // 초기값 설정: [0] Constant (Height), [1] Mean, [2] Sigma
  //fitFunc->SetParameters(size0, mean0, 2.0); // Sigma 초기값은 임의로 2.0으로 부여
  //fitFunc->SetParameters(size0, mean0, 2.0,size0,mean0,2.0); // Sigma 초기값은 임의로 2.0으로 부여
  fitFunc->SetParameters(size0, mean0, 2.0,0.1,0.1,0.1); // Sigma 초기값은 임의로 2.0으로 부여

  // 파라미터 한계치 설정 (잘못된 영역으로 튀는 현상 방지)
  fitFunc->SetParLimits(0,size0*0.5, size0 *1.05); // 높이 상하한선
  fitFunc->SetParLimits(1, fit_min, fit_max);         //sig
  fitFunc->SetParLimits(2, 0.1, 3);                // Sigma 제한

  TF1 *f1 = new TF1("gaus", "gaus", fit_min, fit_max);
  //TF1 *f2 = new TF1("gaus", "gaus", fit_min, fit_max);
  //TF1 *f2 = new TF1("bg", "pol2", fit_min, fit_max);
  // 피팅 실행 (R: 지정 범위, M: 개선된 결과, I: 적분 사용)
  hProjY->Fit(fitFunc, "RMI");

  // ================= [5. 시각화] =================
  hProjY->Draw("hist");
  fitFunc->SetLineColor(kRed);
  
  f1 -> SetParameters(fitFunc -> GetParameter(0), fitFunc -> GetParameter(1), fitFunc -> GetParameter(2));
  //f2 -> SetParameters(fitFunc -> GetParameter(3), fitFunc -> GetParameter(4), fitFunc -> GetParameter(5));
  
  f1 -> SetLineColor(kGreen); f1 -> SetLineStyle(2);
  //f2 -> SetLineColor(kMagenta); f2 -> SetLineStyle(2);
  
  f1->Draw("same");
  //f2->Draw("same");
  fitFunc->Draw("same");
  
  cFit->Update();

  std::cout << Form("Target q%d_t%d: CsI range (%.1f ~ %.1f) projected and fitted.", target_q, target_t, CsI_min, CsI_max) << std::endl;
}

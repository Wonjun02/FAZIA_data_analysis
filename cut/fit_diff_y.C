#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TString.h>

void fit_diff_y()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1111); 
  
  // --- 사용자 설정 파라미터 ---
  int Energy = 70; 
  double fitMin = -4.0;      // 가우시안 피팅을 시작점
  double fitMax = 4.0;       // 가우시안 피팅 끝점
  double searchMin = -3.0;   // 초기 평균값 범위
  double searchMax = 3.0;    // 초기 평균값 범위
  // -------------------------

  TString inFileName = Form("Diff_y_%dMeV.root", Energy);
  TString outFileName = Form("Diff_fitting_%dMeV.root", Energy);

  TFile *inFile = new TFile(inFileName, "READ");
  if (!inFile || inFile->IsZombie()) {
    printf("Error: Cannot open input file %s\n", inFileName.Data());
    return;
  }

  TFile *outFile = new TFile(outFileName, "RECREATE");

  for (int q = 1; q <= 4; q++) {
    for (int t = 1; t <= 4; t++) {
      
      TString histName = Form("histDiff_y_q%d_t%d", q, t);
      TH1D *histDiff_y = (TH1D*)inFile->Get(histName);
      if (!histDiff_y) continue;

      // 가우시안 함수 정의
      TF1 *gaus = new TF1(Form("gaus_q%d_t%d", q, t), "gaus", fitMin, fitMax);
      
      // find peak
      histDiff_y->GetXaxis()->SetRangeUser(searchMin, searchMax);
      double peakPos = histDiff_y->GetXaxis()->GetBinCenter(histDiff_y->GetMaximumBin());
      double peakSize = histDiff_y->GetMaximum();
      
      histDiff_y->GetXaxis()->SetRange(0, 0);

      gaus->SetParameters(peakSize, peakPos, 1.0);

      // "S": 피팅 결과를 객체에 저장하여 캔버스 없이도 통계 박스가 유지되도록 함
      // "R": 함수 범위(fitMin~fitMax) 내에서 피팅
      // "M": 피팅 알고리즘 개선, "Q": 로그 출력 최소화
      histDiff_y->Fit(gaus, "RMIQS"); 

      histDiff_y->SetStats(kTRUE);

      outFile->cd();
      histDiff_y->Write(); 
      gaus->Write();

      delete gaus;
    }
  }
  
  outFile->Close();
  inFile->Close();

  printf("\n==========================================\n");
  printf("  - File Saved: %s\n", outFileName.Data());
  printf("==========================================\n");
}

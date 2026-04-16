void fit_aTel()
{
  // ================= [1. 대상 선택 및 설정] =================
  int Energy = 70; 
  int q = 4;  // 원하는 Quadrant 번호 (1~4)
  int t = 1;  // 원하는 Telescope 번호 (1~4)

  // 피팅을 수행할 전체 x축 범위
  double xFitMin = 115.0, xFitMax = 200.0;

  // 가우시안 피크를 찾을 탐색 구간 (Search Window)
  double g1_search[2] = {145.0, 170.0}, g2_search[2] = {147.0, 152.0}; 

  // ================= [2. 파일 및 히스토그램 로드] =================
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("e");
  gROOT->ForceStyle();

  TFile *file = new TFile(Form("d_cut_histCorr_%dMeV.root", Energy), "READ");
  if (file->IsZombie()) return;

  // 파일 내 이름 규칙에 따라 히스토그램 로드
  TString histName = Form("histProj_x_cut_q%d_t%d", q, t);
  TH1D *histProj = (TH1D*)file->Get(histName);

  if (!histProj) {
    std::cout << "히스토그램 " << histName << "을 찾을 수 없습니다!" << std::endl;
    return;
  }

  // ================= [3. 초기값 자동 추출 및 피팅] =================
  auto c1 = new TCanvas("c1", "Specific Projection Fit", 800, 600);
  
  // 구간 1 피크 찾기
  histProj->GetXaxis()->SetRangeUser(g1_search[0], g1_search[1]);
  double size1 = histProj->GetMaximum();
  double mean1 = histProj->GetXaxis()->GetBinCenter(histProj->GetMaximumBin());
  
  // 구간 2 피크 찾기
  histProj->GetXaxis()->SetRangeUser(g2_search[0], g2_search[1]);
  double size2 = histProj->GetMaximum();
  double mean2 = histProj->GetXaxis()->GetBinCenter(histProj->GetMaximumBin());
  
  histProj->GetXaxis()->SetRange(0, 0); // 범위 원상복구

  // 피팅 함수 정의 (gaus+gaus+pol2)
  TF1 *totalFit = new TF1("totalFit", "gaus(0) + gaus(3) + pol2(6)", xFitMin, xFitMax);
  //TF1 *totalFit = new TF1("totalFit", "gaus(0) + gaus(3) + [8]*(x-[7])*(x-[7])+[6]", xFitMin, xFitMax);
  
  // 초기값 설정
  totalFit->SetParameters(size1, mean1, 1.2, size2, mean2, 1.2, 0, 0, 0);

  // 파라미터 한계치 설정 (튀는 현상 방지)
  totalFit->SetParLimits(0, size1 * 0.2, size1 * 1.05); // 높이 상한선 제어
  totalFit->SetParLimits(1, g1_search[0], g1_search[1]); //mean 제한
  totalFit->SetParLimits(2, 0, 10); //sigma 제한
  totalFit->SetParLimits(3, size2 * 0.2, size2 * 1.05);
  totalFit->SetParLimits(4, g2_search[0], g2_search[1]);
  totalFit->SetParLimits(5, 5, 12);
  totalFit->SetParLimits(6, 1000, 10000); 
  //totalFit->SetParLimits(8, 1e-9, 100); 

  // 피팅 실행
  histProj->Fit(totalFit, "RMI");

  // ================= [4. 시각화] =================
  histProj->SetTitle(Form("Fitting Result: Qua %d Tel %d (%d MeV)", q, t, Energy));
  histProj->Draw("hist");

  // 컴포넌트별 그리기 (G1:초록, G2:보라, BG:회색)
  TF1 *gaus1 = new TF1("gaus1", "gaus", xFitMin, xFitMax);
  TF1 *gaus2 = new TF1("gaus2", "gaus", xFitMin, xFitMax);
  TF1 *bg = new TF1("bg", "pol2", xFitMin, xFitMax);
  gaus1->SetParameters(totalFit->GetParameter(0), totalFit->GetParameter(1), totalFit->GetParameter(2));
  gaus2->SetParameters(totalFit->GetParameter(3), totalFit->GetParameter(4), totalFit->GetParameter(5));
  bg->SetParameters(totalFit->GetParameter(6), totalFit->GetParameter(7), totalFit->GetParameter(8));
  
  gaus1->SetLineColor(kGreen+2); gaus1->SetLineStyle(2);
  gaus2->SetLineColor(kMagenta); gaus2->SetLineStyle(2);
  bg->SetLineColor(kGray+2);  bg->SetLineStyle(3);
  
  gaus1->Draw("same");
  gaus2->Draw("same");
  bg->Draw("same");
  totalFit->Draw("same");

  gPad->Update();
}

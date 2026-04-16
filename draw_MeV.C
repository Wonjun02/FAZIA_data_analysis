#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

void draw_MeV()
{
    // 통계 박스 옵션 설정 (1111: Entries, Mean, RMS, Underflow, Overflow 모두 표시)
    gStyle->SetOptStat(1111);

    int Energy = 60; // 분석 에너지 설정

    // 1. 파일 열기
    TString inFileName = Form("hist_%dMeV_MeV.root", Energy);
    TFile *file = new TFile(inFileName, "READ");
    if (!file || file->IsZombie()) {
        printf("Error: Cannot open %s\n", inFileName.Data());
        return;
    }

    // 2. 저장 경로 설정 및 폴더 생성
    TString path = Form("./MeV_results_%dMeV/", Energy);
    gSystem->Exec("mkdir -p " + path);

    // 3. 모든 텔레스코프 조합(4x4) 루프
    for (int q = 1; q <= 4; q++) {
        for (int t = 1; t <= 4; t++) {
            
            TString histName = Form("histCorr_q%d_t%d", q, t);
            TH2D* h = (TH2D*)file->Get(histName);

            // 히스토그램이 없으면(칼리브레이션 안 한 경우) pass
            if (!h) {
                continue;
            }

            // 4. 개별 캔버스 생성 및 설정
            TCanvas* c = new TCanvas(Form("c_%s", histName.Data()), histName, 900, 750);
            c->SetRightMargin(0.15); // 통계 박스와 팔레트 공간 확보
            //c->SetGrid();
            c->SetTicks();
            c->SetLogz(); //

            // 5. 그리기 및 스타일 조정
            h->SetStats(1); // 통계 박스 활성화
            h->Draw("colz");
            h->SetTitle(Form("Q%dT%d (%dMeV)", q, t, Energy));
            
            // 축 레이블 크기 조정
            h->GetXaxis()->SetTitleSize(0.04);
            h->GetYaxis()->SetTitleSize(0.04);
            h->GetXaxis()->SetLabelSize(0.03);
            h->GetYaxis()->SetLabelSize(0.03);

            // 6. 개별 PDF로 저장
            c->SaveAs(path + histName + "_MeV.pdf");
            
            //delete c; 
            
            std::cout << ">>> Saved: " << histName << "_MeV.pdf with Stats(1111)" << std::endl;
        }
    }

    //file->Close();
    std::cout << "\nAll existing MeV histograms have been saved to " << path << std::endl;
}

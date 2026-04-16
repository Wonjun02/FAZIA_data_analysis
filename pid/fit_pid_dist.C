#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>

using namespace std;

void fit_pid_dist()
{
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat("e");
    gROOT->ForceStyle();

    int Energy = 70;
    double xFitMin = 0.0, xFitMax = 4.0;
    double g1_search[2] = {0.7, 1.7};
    //double g2_search[2] = {1.4, 3.0};
    double g2_search[2] = {1.8, 3.0};

    TFile *fIn = new TFile(Form("pid_dist_%dMeV.root", Energy), "READ");
    if (!fIn || fIn->IsZombie()) return;

    TFile *fOut = new TFile(Form("fit_pid_dist_%dMeV.root", Energy), "RECREATE");
    TString plotPath = Form("./PID_Fit_Results_%dMeV/", Energy);
    gSystem->Exec("mkdir -p " + plotPath);

    for (int q = 1; q <= 4; q++) {
        for (int t = 1; t <= 4; t++) {
            TString histName = Form("hist_PID_q%d_t%d", q, t);
            TH1D *h_pid = (TH1D*)fIn->Get(histName);
            if (!h_pid) continue;

            // 1. 초기값 추출
            h_pid->GetXaxis()->SetRangeUser(g1_search[0], g1_search[1]);
            double size1 = h_pid->GetMaximum();
            double mean1 = h_pid->GetXaxis()->GetBinCenter(h_pid->GetMaximumBin());
            
            h_pid->GetXaxis()->SetRangeUser(g2_search[0], g2_search[1]);
            double size2 = h_pid->GetMaximum();
            double mean2 = h_pid->GetXaxis()->GetBinCenter(h_pid->GetMaximumBin());
            h_pid->GetXaxis()->SetRange(0, 0);

            // 2. 피팅 함수 정의 (이름 중복 방지를 위해 루프 인덱스 포함)
            TString fTotalName = Form("Fit_q%d_t%d", q, t);
            TF1 *totalFit = new TF1(fTotalName, "gaus(0) + gaus(3) + pol2(6)", xFitMin, xFitMax);
            totalFit->SetParameters(size1, mean1, 0.1, size2, mean2, 0.15, 0, 0, 0);
            totalFit->SetParLimits(0, size1*0.2, size1*1.05);
            totalFit->SetParLimits(1, g1_search[0], g1_search[1]);
            totalFit->SetParLimits(3, size2*0.2, size2*1.05);
            totalFit->SetParLimits(4, g2_search[0], g2_search[1]);

            // "0" 옵션: 피팅 결과를 자동으로 그리지 않음 (중복 창 방지)
            h_pid->Fit(totalFit, "RMIQ0");

            // 3. 수동으로 캔버스 생성 및 그리기
            TCanvas *c = new TCanvas(Form("c_fit_q%d_t%d", q, t), Form("Fit Q%dT%d", q, t), 800, 600);
            h_pid->Draw("hist");

            // 4. 개별 컴포넌트 함수 생성 (이름 고유화)
            TF1 *fG1 = new TF1(Form("fG1_q%dt%d", q, t), "gaus", xFitMin, xFitMax);
            TF1 *fG2 = new TF1(Form("fG2_q%dt%d", q, t), "gaus", xFitMin, xFitMax);
            TF1 *fBG = new TF1(Form("fBG_q%dt%d", q, t), "pol2", xFitMin, xFitMax);

            fG1->SetParameters(totalFit->GetParameter(0), totalFit->GetParameter(1), totalFit->GetParameter(2));
            fG2->SetParameters(totalFit->GetParameter(3), totalFit->GetParameter(4), totalFit->GetParameter(5));
            fBG->SetParameters(totalFit->GetParameter(6), totalFit->GetParameter(7), totalFit->GetParameter(8));

            fG1->SetLineColor(kGreen+2); fG1->SetLineStyle(2);
            fG2->SetLineColor(kMagenta); fG2->SetLineStyle(2);
            fBG->SetLineColor(kGray+2);  fBG->SetLineStyle(3);

            // 5. 현재 캔버스에 모두 그리기
            fG1->Draw("same");
            fG2->Draw("same");
            fBG->Draw("same");
            totalFit->Draw("same");

            c->Update();
            c->Modified();

            // 6. 파일 및 PDF 저장
            fOut->cd();
            c->Write();
            h_pid->Write(Form("fit_hist_pid_q%d_t%d", q, t));
            c->SaveAs(plotPath + Form("fit_pid_q%d_t%d_%dMeV.pdf", q, t, Energy));
            
            cout << ">>> Finished: Q" << q << "T" << t << endl;
        }
    }

    fOut->Close();
    cout << "\nAnalysis Complete. Result saved in fit_pid_dist_" << Energy << "MeV.root" << endl;
}

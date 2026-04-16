#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

void pid_dist()
{
    gStyle->SetOptStat(0);
    int Energy = 60; // 현재 분석 중인 에너지이자 x의 최대값 기준
    double x_min_cut = 10.0;
    double x_max_cut = (double)Energy-20.; // Energy 변수를 x 최대값으로 사용

    // 1. Reference 피팅 함수 로드
    TFile *fRef = new TFile("ref_PID.root", "READ");
    TF1 *ft_p = (TF1*)fRef->Get("ft_proton");
    TF1 *ft_d = (TF1*)fRef->Get("ft_deutron");

    // 2. 실험 데이터 로드
    TString inputName = Form("hist_%dMeV_MeV.root", Energy);
    TFile* fIn = new TFile(inputName, "READ");
    if (!fIn || fIn->IsZombie()) {
        cout << "Error: Cannot find " << inputName << endl;
        return;
    }

    // 출력 파일 및 폴더 설정
    TString outputFileName = Form("pid_dist_%dMeV.root", Energy);
    TFile *fOut = new TFile(outputFileName, "RECREATE");

    TString plotDirPath = Form("./PID_Check_Plots_%dMeV/", Energy);
    gSystem->Exec("mkdir -p " + plotDirPath);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int q = i + 1;
            int t = j + 1;

            TString histName = Form("histCorr_q%d_t%d", q, t);
            TH2D* hist_Corr = (TH2D*)fIn->Get(histName);
            if (!hist_Corr) continue;

            hist_Corr->SetDirectory(0); 

            // PID 히스토그램 생성
            TString hPidName = Form("hist_PID_q%d_t%d", q, t);
            TString hPidTitle = Form("PID_q%dt%d (%dMeV, X:%.0f-%.0f);PID Index;Counts", q, t, Energy, x_min_cut, x_max_cut);
            TH1D* h_pid = new TH1D(hPidName, hPidTitle, 500, 0, 5);
            
            h_pid->SetDirectory(0);
            h_pid->SetOption("hist");
            h_pid->SetLineWidth(2);

            TCanvas *c_check = new TCanvas(Form("c_check_q%d_t%d", q, t), "PID Mapping Check", 800, 600);
            c_check->SetLogz();
            hist_Corr->Draw("colz");

            //분석 data
            int nx = hist_Corr->GetNbinsX();
            int ny = hist_Corr->GetNbinsY();

            for (int ix = 1; ix <= nx; ix++) {
                double x_val = hist_Corr->GetXaxis()->GetBinCenter(ix);
                
                //pid진행할 x범위 
                if (x_val < x_min_cut || x_val > x_max_cut) continue;

                double y_p = ft_p->Eval(x_val);
                double y_d = ft_d->Eval(x_val);
                double dist_pd = y_d - y_p;

                if (dist_pd <= 0) continue; 

                for (int iy = 1; iy <= ny; iy++) {
                    double content = hist_Corr->GetBinContent(ix, iy);
                    if (content <= 0) continue;

                    double y_val = hist_Corr->GetYaxis()->GetBinCenter(iy);
                    if (y_val < (y_p - dist_pd)) continue; //noise를 포함하지 않기 위해 pid의 y값 제한

                    double pid_idx = 1.0 + (y_val - y_p) / dist_pd;
                    if (pid_idx > 5.0) continue; 

                    h_pid->Fill(pid_idx, content);
                }
            }

            // --- 가이드라인 시각화 ---
            for(int k=0; k<6; k++) {
                double target_pid = (double)k;
		// x is variable ->x[0] 1개만 존재, p는 parameter, p는 외부변수를 가져오기에  # of p =0
                TF1* fLine = new TF1(Form("fLine_%d_q%dt%d", k, q, t), [=](double* x, double* p){
                    double yp = ft_p->Eval(x[0]); //x에 따른 proton line 위치 계산 
                    double yd = ft_d->Eval(x[0]); //x에 따른 deutron line 위치 계산
                    return yp + (target_pid - 1.0) * (yd - yp); //target pid 구간 정의
                }, x_min_cut, x_max_cut, 0); 
                
                fLine->SetLineColor(kRed + k);
                fLine->SetLineStyle(k == 1 || k == 2 ? 1 : 7);
                fLine->Draw("same");
            }

            fOut->cd();
            h_pid->Write();
            c_check->Write();

            c_check->SaveAs(plotDirPath + Form("check_q%d_t%d_%dMeV.pdf", q, t, Energy));
            
            cout << ">>> Q" << q << "T" << t << " Processed (" << Energy << "MeV, X:" << x_min_cut << "-" << x_max_cut << ")" << endl;
        }
    }

    fOut->cd();
    fOut->Close();
    cout << "\nPID distribution for " << Energy << "MeV saved to: " << outputFileName << endl;
}

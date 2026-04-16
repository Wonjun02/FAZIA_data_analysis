#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <iostream>

using namespace std;

void test_cut()
{
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1111);

    int Energy = 70; 
    
    TFile::Open(Form("Diff_y_%dMeV.root", Energy));
    TFile* fFitResult = new TFile(Form("Diff_fitting_%dMeV.root", Energy), "READ");
    if (!fFitResult || fFitResult->IsZombie()) {
        cout << "Error: Cannot find Diff_fitting_" << Energy << "MeV.root" << endl;
        return;
    }

    TF1* fTG[4][4];
    TF1* fGaus[4][4];
    TH1D* hDiff_y[4][4];
    double cut_limit_up[4][4];
    double cut_limit_down[4][4];

    // Projection histograms 추가
    TH1D* histProj_x_cut[4][4];
    TH1D* histProj_y_cut[4][4];

    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            int q = i + 1;
            int t = j + 1;
            // fTG from Diff_y_60MeV.root, gaus from Diff_y_60MeV
            fTG[i][j] = (TF1*)gROOT->GetFile(Form("Diff_y_%dMeV.root", Energy))->Get(Form("fTG_q%d_t%d", q, t));
            fGaus[i][j] = (TF1*)fFitResult->Get(Form("gaus_q%d_t%d", q, t));
            hDiff_y[i][j] = (TH1D*)fFitResult->Get(Form("histDiff_y_q%d_t%d", q, t));

            if (fGaus[i][j]) {
                // mean - 1*sigma > data -> cut
                double mean = fGaus[i][j]->GetParameter(1);
                double sigma = fGaus[i][j]->GetParameter(2);
                cut_limit_down[i][j] = mean - sigma;
                cut_limit_up[i][j] = mean + sigma;
            } else {
                cut_limit_down[i][j] = -999.0; //
                cut_limit_up[i][j] = 999.0; //
            }

            // draw_PID.C와 동일한 범위로 프로젝션 히스토그램 초기화
            histProj_x_cut[i][j] = new TH1D(Form("histProj_x_cut_q%d_t%d", q, t), Form("cut_x_Projection_q%d_t%d;sQ3max;counts", q, t), 250, 50, 300);
            histProj_y_cut[i][j] = new TH1D(Form("histProj_y_cut_q%d_t%d", q, t), Form("cut_y_Projection_q%d_t%d;sQ2max;counts", q, t), 100, 0, 100);
        }
    }

    // data path
    TChain chain("faziatree");
    //int startRun = 43, endRun = 381; // 50MeV_1 /day2
    //int startRun = 382, endRun = 545; // 50MeV_2 /day2
    //int startRun = 546, endRun = 589; // 50MeV_3 /day2
    //int startRun = 43, endRun = 589; // total 50MeV /day2
    //int startRun = 590, endRun = 1095; // 40MeV /day2
    int startRun = 340, endRun = 1032; // 70MeV /day1
    //int startRun = 1033, endRun = 1699; // 60MeV /day1

    for (int run=startRun; run<=endRun; run++) 
        chain.Add(Form("/media/wonjun/T7/converted_data/run_%d.root", run));

    UShort_t mtot, fQua[100], fTel[100];
    float sQ3max[100], sQ2max[100];
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("Mtot", 1); chain.SetBranchAddress("Mtot", &mtot);
    chain.SetBranchStatus("fQua", 1); chain.SetBranchAddress("fQua", fQua);
    chain.SetBranchStatus("fTel", 1); chain.SetBranchAddress("fTel", fTel);
    chain.SetBranchStatus("sQ3max", 1); chain.SetBranchAddress("sQ3max", sQ3max);
    chain.SetBranchStatus("sQ2max", 1); chain.SetBranchAddress("sQ2max", sQ2max);

    // histCorr and histCorrWithCut
    TH2D* histCorr[4][4];
    TH2D* histCorr_cut[4][4];
    
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            TString title = Form("q%d_t%d", i+1, j+1);
            histCorr[i][j] = new TH2D("histCorr_"+title, "Corr_"+title+";sQ3max;sQ2max", 600, 0, 300, 400, 0, 80);
            histCorr_cut[i][j] = new TH2D("histCorr_cut_"+title, "cut_corr_"+title+";sQ3max;sQ2max", 600, 0, 300, 400, 0, 80);
        }
    }

    Long64_t nEntries = chain.GetEntries();
    cout << "Total Entries: " << nEntries << endl;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        chain.GetEntry(iEntry);
        int nHit = (mtot > 100) ? 100 : (int)mtot;

        for (int i = 0; i < nHit; i++) {
            int q = fQua[i] - 1;
            int t = fTel[i] - 1;
            if (q < 0 || q >= 4 || t < 0 || t >= 4) continue;
            if (!fTG[q][t]) continue;

            float x = sQ3max[i];
            float y = sQ2max[i];
            double dy = y - fTG[q][t]->Eval(x);

            // histCorr without cut
            histCorr[q][t]->Fill(x, y);

            // dy > (mean - 1*sigma)
            if (dy > cut_limit_down[q][t] && dy < cut_limit_up[q][t]) {
                histCorr_cut[q][t]->Fill(x, y);
                // 컷 통과 데이터에 대해 프로젝션 수행
                histProj_x_cut[q][t]->Fill(x);
                histProj_y_cut[q][t]->Fill(y);
            }
        }
        if(iEntry % 500000 == 0) printf("Processing... %.1f%%\n", 100.0*iEntry/nEntries);
    }

    //TFile* outputFile = new TFile(Form("cut_histCorr_%dMeV.root", Energy), "RECREATE");
    TFile* outputFile = new TFile(Form("d_cut_histCorr_%dMeV.root", Energy), "RECREATE");
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            if(fTG[i][j]) fTG[i][j]->Write();
            if(fGaus[i][j]) fGaus[i][j]->Write();
            if(histCorr[i][j]) histCorr[i][j]->Write();
            if(hDiff_y[i][j]) hDiff_y[i][j]->Write();
            if(histCorr_cut[i][j]) histCorr_cut[i][j]->Write();
            // 프로젝션 히스토그램 저장
            if(histProj_x_cut[i][j]) histProj_x_cut[i][j]->Write();
            if(histProj_y_cut[i][j]) histProj_y_cut[i][j]->Write();
        }
    }

    outputFile->Close();
    fFitResult->Close();
    
    cout << "==========================================" << endl;
    cout << " Output: cut_histCorr_" << Energy << "MeV.root" << endl;
    cout << "==========================================" << endl;
}

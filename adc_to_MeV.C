#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <iostream>

using namespace std;

//Si (sQ2max) 변환 함수: pol2
double Calibrate_Si(double adc, double p0, double p1) {
    //return (p0 * adc * adc) + (p1 * adc);
    return (p0 * adc) + (p1);
}

// CsI (sQ3max) 변환 함수: pol2
double Calibrate_CsI(double adc, double p0, double p1) {
    //return (p0 * adc * adc) + (p1 * adc);
    return (p0 * adc) + (p1);
}

void adc_to_MeV()
{
    int Energy = 60; 
    

    TString inputName = Form("d_cut_histCorr_%dMeV.root", Energy);
    TFile* fIn = new TFile(inputName, "READ");
    if (!fIn || fIn->IsZombie()) {
        cout << "Error: Cannot find " << inputName << endl;
        return;
    }

    // 결과 저장 파일
    TString outputName = Form("hist_%dMeV_MeV.root", Energy);
    TFile* fOut = new TFile(outputName, "RECREATE");

    // =================칼리브레이션 파라미터=================
    double si_p0[4][4], si_p1[4][4];
    double csi_p0[4][4], csi_p1[4][4];

    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            si_p0[i][j] = 0.0; si_p1[i][j] = 0.0;
            csi_p0[i][j] = 0.0; csi_p1[i][j] = 0.0;
        }
    }


    // Q1T2 (i=0, j=1)
    csi_p0[0][1] = 0.38559; csi_p1[0][1] = 0.0;
    si_p0[0][1] = 0.27091; si_p1[0][1] = 0.0;
    // Q4T4 (i=3, j=3)
    csi_p0[3][3] = 0.31677; csi_p1[3][3] = 0.0;
    si_p0[3][3] = 0.26171; si_p1[3][3] = 0.0;
    // Q1T1 (i=0, j=0)
    csi_p0[0][0] = 0.38951; csi_p1[0][0] = 0.0;
    si_p0[0][0] = 0.25190; si_p1[0][0] = 0.0;
    // Q1T3 (i=0, j=2)
    csi_p0[0][2] = 0.35762; csi_p1[0][2] = 0.0;
    si_p0[0][2] = 0.30027; si_p1[0][2] = 0.0;
    // Q3T2 (i=2, j=1)
    csi_p0[2][1] = 0.46902; csi_p1[2][1] = 0.0;
    si_p0[2][1] = 0.28011; si_p1[2][1] = 0.0;
    

    // =================================================================

    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            int q = i + 1;
            int t = j + 1;

            if (si_p0[i][j] == 0.0 && csi_p0[i][j] == 0.0 && si_p1[i][j] == 0.0 && csi_p1[i][j] == 0.0) continue;

            TString histName = Form("histCorr_q%d_t%d", q, t);
            TH2D* h_adc = (TH2D*)fIn->Get(histName);
            if (!h_adc) continue;

            //MeV histogram
            TH2D* h_MeV = new TH2D(histName, 
                                  Form("histCorr_q%d_t%d_MeV;sQ3max (MeV);sQ2max (MeV)", q, t),
                                  600, 0, 100, 400, 0, 30);
            h_MeV->SetOption("COLZ");
            for (int bx = 1; bx <= h_adc->GetNbinsX(); bx++) {
                double adc_x = h_adc->GetXaxis()->GetBinCenter(bx); 
                // CsI channel convert
                double MeV_x = Calibrate_CsI(adc_x, csi_p0[i][j], csi_p1[i][j]);

                for (int by = 1; by <= h_adc->GetNbinsY(); by++) {
                    double content = h_adc->GetBinContent(bx, by);
                    if (content <= 0) continue;

                    double adc_y = h_adc->GetYaxis()->GetBinCenter(by); 
                    // Si channel convert
                    double MeV_y = Calibrate_Si(adc_y, si_p0[i][j], si_p1[i][j]);

                    h_MeV->Fill(MeV_x, MeV_y, content);
                }
            }

            fOut->cd();
            h_MeV->Write();
            cout << ">>> Telescope Q" << q << "T" << t << " converted." << endl;
        }
    }

    fOut->Close();
    fIn->Close();
    cout << "\nJob Done. Output file: " << outputName << endl;
}

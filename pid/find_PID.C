#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>

void find_PID()
{
    gStyle->SetOptStat(0);

    // 1. 설정
    int Energy = 70;
    const char* particles[2] = {"proton", "deutron"};
    double mass_ratios[2] = {1.0, 2.0}; // Z=1(p), Z=1, A=2(d)
    int colors[2] = {kRed+1, kBlue+1};

    // 결과 저장용 파일
    auto file_out = new TFile("ref_PID.root", "RECREATE");

    for (int p = 0; p < 2; p++) {
        // 2. 파일 및 히스토그램 로드
        TString fileName = Form("hist_%s_%dMeV.root", particles[p], Energy);
        TFile* file_input = new TFile(fileName, "READ");
        if (file_input->IsZombie()) {
            cout << "Error: Cannot find " << fileName << endl;
            continue;
        }

        // test.C에서 생성된 전체 합산 히스토그램 "hist_edep" 로드
        TH2D* hist_edep = (TH2D*)file_input->Get("hist_edep");
        if (!hist_edep) continue;

        // 3. 피팅을 위한 데이터 포인트 추출 (Profile 대신 각 X빈의 Peak 지점 탐색)
        TGraph* grPoints = new TGraph();
        int nPoints = 0;
        for (int ix = 1; ix <= hist_edep->GetNbinsX(); ix++) {
            double x_val = hist_edep->GetXaxis()->GetBinCenter(ix);
            if (x_val < 1.2 || x_val > 70) continue;
            // 해당 X 구간의 Y 투영 (Si2)
            TH1D* Proj_y = hist_edep->ProjectionY("Proj_y", ix, ix) ; //ProjectionY("name",firstbin,lastbin)
            if (Proj_y->GetEntries() < 10) { delete Proj_y; continue; }

            double y_max = Proj_y->GetBinCenter(Proj_y->GetMaximumBin()); //각 CsI 에너지 영역별 가장 많은 entry가 있는 si에너지값 사용
            if (y_max > 0.1) {
                grPoints->SetPoint(nPoints++, x_val, y_max);
            }
            delete Proj_y;
        }

        // 4. 피팅 함수 정의 (draw_PID.C의 수식 사용)
        // [0]:g, [1]:mu, [2]:lambda, [3]:Z, [4]:A, [5]:pdy
        TF1* ft = new TF1(Form("ft_%s", particles[p]),
            "pow( pow([0]*x,[1]+1) + pow([2],[1]+1)*pow([3],2)*pow([4],[1]) , 1./([1]+1) ) - [0]*x + [5]",
            0, 80);

        ft->SetParNames("g","mu","lambda","Z","A","pdy");
        ft->SetParameters(1.0, 0.6, 20.0, 1.0, mass_ratios[p], 0.0);
        
        // Z와 A는 고정 (p=1,1 / d=1,2)
        ft->FixParameter(3, 1.0); 
        ft->FixParameter(4, mass_ratios[p]);

        // 피팅 수행
        grPoints->Fit(ft, "RQ");

        // 5. 시각화 및 PDF 저장 (fit_qxtx_CsI 스타일)
        TCanvas* c = new TCanvas(Form("c_%s", particles[p]), particles[p], 800, 650);
        c->SetGrid();
        c->SetLogz();
        
        hist_edep->SetTitle(Form("%s PID Fit (%d MeV)", particles[p], Energy));
        hist_edep->GetXaxis()->SetTitle("Edep_{CsI} (MeV)");
        hist_edep->GetYaxis()->SetTitle("Edep_{Si2} (MeV)");
        hist_edep->Draw("colz");

        ft->SetLineColor(colors[p]);
        ft->SetLineWidth(3);
        ft->Draw("SAME");

        // 파라미터 표시 (fit_qxtx_CsI 스타일)
        double g      = ft->GetParameter(0);
        double mu     = ft->GetParameter(1);
        double lambda = ft->GetParameter(2);
        double pdy    = ft->GetParameter(5);
        double chi2   = ft->GetChisquare();
        int ndf       = ft->GetNDF();

        TLegend* leg = new TLegend(0.45, 0.75, 0.85, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(ft, Form("Fit: %s", particles[p]), "l");
        leg->Draw();

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.035);
        double ty = 0.70;
        tex->DrawLatex(0.45, ty,      Form("g = %.4f, #mu = %.4f", g, mu));
        tex->DrawLatex(0.45, ty-0.05, Form("#lambda = %.4f, pdy = %.4f", lambda, pdy));
        tex->DrawLatex(0.45, ty-0.10, Form("#chi^{2}/ndf = %.2f / %d", chi2, ndf));

        c->SaveAs(Form("fit_ref_%s.pdf", particles[p]));

        // 6. 결과 저장 (ref_PID.root)
        file_out->cd();
        // 히스토그램 이름을 구분하여 저장
        hist_edep->Write(Form("hist_%s", particles[p]));
        ft->Write(); 

        cout << "Finished processing " << particles[p] << endl;
        // file_input->Close(); 
    }

    cout << "\n>>> All results saved to ref_PID.root" << endl;
    // file_out->Close();
}

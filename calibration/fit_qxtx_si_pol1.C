// fit_qxtx_si_pol1.C
void fit_qxtx_si_pol1()
{
    gStyle->SetOptStat(0);

    const int n = 4;
    const int nGraph = 8;

    // Y values (Energy -> MeV)
    double E[n] = {1.241, 1.323, 1.42, 1.537}; //70H, 70L, 60H, 60L

    // Y축(MeV)의 에러(sigma) 값
    double sigma_E[n] = {0.08827, 0.09062, 0.08901, 0.1038};

    // X values (Channel -> adc)
    double Q1T2[n] = {4.589, 4.792, 5.281, 5.716};
    double Q4T4[n] = {4.687, 5.029, 5.447, 5.936};
    
    double Q1T1[n] = {4.956, 5.371, 5.651, 5.908};
    double Q1T3[n] = {4.064, 4.434, 4.674, 5.229};
    double Q2T3[n] = {3.73, 3.638, 3.957, 3.536};
    double Q2T4[n] = {3.568, 3.393, 3.637, 3.391};
    double Q3T2[n] = {4.59, 4.741, 4.929, 5.454};
    double Q4T1[n] = {4.463, 4.369, 4.925, 4.695};

    // X축(adc) 에러(sigma)
    double sigma_Q1T2[n] = {0., 0., 0., 0.}; 
    double sigma_Q4T4[n] = {0., 0., 0., 0.};
    
    double sigma_Q1T1[n] = {0., 0., 0., 0.};
    double sigma_Q1T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T4[n] = {0., 0., 0., 0.};
    double sigma_Q3T2[n] = {0., 0., 0., 0.};
    double sigma_Q4T1[n] = {0., 0., 0., 0.};

    const char* names[nGraph] = 
    {
        "Q1T2", "Q4T4", "Q1T1", "Q1T3",
        "Q2T3", "Q2T4", "Q3T2", "Q4T1"
    };

    double* adc_vals[nGraph] = 
    {
        Q1T2, Q4T4, Q1T1, Q1T3,
        Q2T3, Q2T4, Q3T2, Q4T1
    };

    double* sigma_adc_vals[nGraph] = 
    {
        sigma_Q1T2, sigma_Q4T4, sigma_Q1T1, sigma_Q1T3,
        sigma_Q2T3, sigma_Q2T4, sigma_Q3T2, sigma_Q4T1
    };

    int colors[nGraph] = 
    {
        kBlue+1, kRed+1, kGreen+2, kMagenta+1,
        kOrange+7, kCyan+2, kBlack, kViolet
    };

    for (int i = 0; i < nGraph; i++) {
        TCanvas* c = new TCanvas(Form("c_%s", names[i]), names[i], 800, 650);
        c->SetGrid();
        c->SetTicks();

        TGraphErrors* gr = new TGraphErrors(n, adc_vals[i], E, sigma_adc_vals[i], sigma_E);
        gr->SetName(Form("gr_%s", names[i]));
        gr->SetTitle("");

        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]); 
        gr->SetMarkerSize(1.5);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        gr->GetXaxis()->SetTitle("adc");
        gr->GetYaxis()->SetTitle("MeV");

        double min_adc = TMath::MinElement(n, adc_vals[i]);
        double max_adc = TMath::MaxElement(n, adc_vals[i]);
        
        double fit_min = min_adc - 0.5;
        double fit_max = max_adc + 0.5;

        gr->GetXaxis()->SetLimits(fit_min, fit_max);
        gr->GetYaxis()->SetRangeUser(1.1, 1.7);

        gr->GetXaxis()->SetTitleSize(0.045);
        gr->GetYaxis()->SetTitleSize(0.045);
        gr->GetXaxis()->SetLabelSize(0.035);
        gr->GetYaxis()->SetLabelSize(0.035);

        gr->GetXaxis()->SetTitleOffset(1.0);
        gr->GetYaxis()->SetTitleOffset(1.15);

        // 1차 함수(pol1) 피팅으로 변경
        TF1* f1 = new TF1(Form("f_%s", names[i]), "pol1", fit_min, fit_max);
        f1->SetLineColor(colors[i]);
        f1->SetLineWidth(2);
        f1->FixParameter(0, 0.0);
        gr->Fit(f1, "RMI");

        gr->Draw("P SAME");

        // pol1에서는 p0가 y절편, p1이 기울기
        double p0 = f1->GetParameter(0);
        double p1 = f1->GetParameter(1);
        double chi2   = f1->GetChisquare();
        int ndf       = f1->GetNDF();

        TLegend* leg = new TLegend(0.12, 0.80, 0.35, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->AddEntry(gr, Form("%s", names[i]), "lp");
        leg->Draw();
        
        // 1차 함수 수식으로 텍스트 변경
        TString formulaStr = TString::Format("E = (%.5f)adc + (%.4f)", p1, p0);
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->DrawLatex(0.12, 0.73, formulaStr.Data());
        tex->DrawLatex(0.12, 0.60, TString::Format("#chi^{2}/ndf = %.2f / %d", chi2, ndf).Data());

        c->Update();
        
        // PDF 파일명에 _pol1 추가 (기존에는 fit_%s.pdf 였으므로 fit_si_pol1_%s.pdf로 통일성을 줌)
        c->SaveAs(Form("fit_si_pol1_%s.pdf", names[i]));

        std::cout << names[i]
                  << " : p1 (slope) = " << p1
                  << ", p0 (intercept) = " << p0
                  << ", Fit Range = [" << fit_min << ", " << fit_max << "]"
                  << std::endl;
    }
}

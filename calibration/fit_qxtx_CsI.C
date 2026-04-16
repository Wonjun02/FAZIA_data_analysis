// fit_qxtx_CsI.C
void fit_qxtx_CsI()
{
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1111);

    const int n = 4;
    const int nGraph = 8;

    // Y values (Energy -> MeV)
    double E[n] = {65.98, 60.62, 55.8, 50.82};
    // ★ 수정된 부분: Y축(MeV) 에러 (sigma) ★
    double sigma_E[n] = {0.1964, 1.323, 0.1876, 1.437};

    // X values (Channel -> adc)
    double Q1T2[n] = {172.5, 155.0, 143.3, 127.6};
    double Q4T4[n] = {208.3, 190.0, 176.2, 158.8};
    double Q1T1[n] = {168.3, 150.0, 144.5, 132.7};
    double Q1T3[n] = {187.8, 169.8, 152.4, 136.2};
    double Q2T3[n] = {173.2, 156.7, 148.6, 132.6};
    double Q2T4[n] = {137.4, 120.2, 115.0, 105.0};
    double Q3T2[n] = {140.6, 129.5, 119.1, 105.0};
    double Q4T1[n] = {162.4, 150.3, 142.7, 129.2};

    // X축(adc) 에러 (sigma)
    double sigma_Q1T2[n] = {0., 0., 0., 0.};
    double sigma_Q4T4[n] = {0., 0., 0., 0.};
    double sigma_Q1T1[n] = {0., 0., 0., 0.};
    double sigma_Q1T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T4[n] = {0., 0., 0., 0.};
    double sigma_Q3T2[n] = {0., 0., 0., 0.};
    double sigma_Q4T1[n] = {0., 0., 0., 0.};

    const char* names[nGraph] = {
        "Q1T2", "Q4T4", "Q1T1", "Q1T3",
        "Q2T3", "Q2T4", "Q3T2", "Q4T1"
    };

    double* adc_vals[nGraph] = {
        Q1T2, Q4T4, Q1T1, Q1T3,
        Q2T3, Q2T4, Q3T2, Q4T1
    };

    double* sigma_adc_vals[nGraph] = {
        sigma_Q1T2, sigma_Q4T4, sigma_Q1T1, sigma_Q1T3,
        sigma_Q2T3, sigma_Q2T4, sigma_Q3T2, sigma_Q4T1
    };

    // 색상은 그래프 구분을 위해 유지
    int colors[nGraph] = {
        kBlue+1, kRed+1, kGreen+2, kMagenta+1,
        kOrange+7, kCyan+2, kBlack, kViolet
    };

    // TCanvas* cPdf = new TCanvas("cPdf", "multiPDF", 800, 650);
    // cPdf->Print("fit_all_pages_CsI.pdf[");

    for (int i = 0; i < nGraph; i++) {
        TCanvas* c = new TCanvas(Form("c_%s", names[i]), names[i], 800, 650);
        c->SetGrid();
        c->SetTicks();

        TGraphErrors* gr = new TGraphErrors(n, adc_vals[i], E, sigma_adc_vals[i], sigma_E);
        gr->SetName(Form("gr_%s", names[i]));
        gr->SetTitle("");

        // 마커를 20번(꽉 찬 동그라미)으로 통일
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]);
        gr->SetMarkerSize(1.5);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        gr->GetXaxis()->SetTitle("adc");
        gr->GetYaxis()->SetTitle("MeV");

        // 현재 그래프 데이터의 최소/최대값 찾기
        double min_adc = TMath::MinElement(n, adc_vals[i]);
        double max_adc = TMath::MaxElement(n, adc_vals[i]);
        
        // 피팅 구간: 최솟값-10 ~ 최댓값+10
        double fit_min = min_adc - 10.0;
        double fit_max = max_adc + 10.0;

        // 화면에 보여지는 X축 범위도 피팅 구간에 맞춰서 유동적으로 조절 (양쪽으로 5씩 여유)
        gr->GetXaxis()->SetLimits(fit_min - 5.0, fit_max + 5.0);
        //gr->GetXaxis()->SetLimits(0., 200.);
        // Y축(에너지)은 에러바를 고려해 범위를 조금 더 넓혔습니다.
        gr->GetYaxis()->SetRangeUser(45.0, 70.0);
        //gr->GetYaxis()->SetRangeUser(0., 70.0);

        gr->GetXaxis()->SetTitleSize(0.045);
        gr->GetYaxis()->SetTitleSize(0.045);
        gr->GetXaxis()->SetLabelSize(0.035);
        gr->GetYaxis()->SetLabelSize(0.035);

        gr->GetXaxis()->SetTitleOffset(1.0);
        gr->GetYaxis()->SetTitleOffset(1.15);

        // 피팅 함수
        TF1* f1 = new TF1(Form("f_%s", names[i]), "[0]*x^2+[1]*x", fit_min, fit_max);
        //TF1* f1 = new TF1(Form("f_%s", names[i]), "[0]*x^2+[1]*x", 0., 200.);
        f1->SetLineColor(colors[i]);
        f1->SetLineWidth(2);

        gr->Fit(f1, "RMI");

        gr->Draw("P SAME");

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

        TString formulaStr = TString::Format("E = (%.5f)adc^{2} + (%.4f)adc", p0, p1);
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->DrawLatex(0.12, 0.73, formulaStr.Data()); 
        tex->DrawLatex(0.12, 0.60, TString::Format("#chi^{2}/ndf = %.2f / %d", chi2, ndf).Data());

        c->Update();
        
         c->SaveAs(Form("fit_CsI_%s.pdf", names[i]));

        std::cout << names[i]
                  << " : p0 (x^2) = " << p0
                  << ", p1 (x) = " << p1
                  << ", Fit Range = [" << fit_min << ", " << fit_max << "]"
                  << std::endl;
    }

    // cPdf->Print("fit_all_pages_CsI.pdf]");
}

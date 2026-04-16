// fit_qxtx_si_modified_with_errors.C
void fit_qxtx_si()
{
    gStyle->SetOptStat(0);

    const int n = 4;
    const int nGraph = 8; // 2에서 8로 변경

    // Y values (Energy -> MeV)
    double E[n] = {1.241, 1.323, 1.42, 1.537}; //70H, 70L, 60H, 60L

    // Y축(MeV)의 에러(sigma) 값
    double sigma_E[n] = {0.08827, 0.09062, 0.08901, 0.1038};

    // X values (Channel -> adc)
    double Q1T2[n] = {4.589, 4.792, 5.281, 5.716};
    double Q4T4[n] = {4.687, 5.029, 5.447, 5.936};
    
    // 나머지 텔레스코프 X축(adc) 배열
    double Q1T1[n] = {4.956, 5.371, 5.651, 5.908};
    double Q1T3[n] = {4.064, 4.434, 4.674, 5.229};
    double Q2T3[n] = {3.73, 3.638, 3.957, 3.536};
    double Q2T4[n] = {3.568, 3.393, 3.637, 3.391};
    double Q3T2[n] = {4.59, 4.741, 4.929, 5.454};
    double Q4T1[n] = {4.463, 4.369, 4.925, 4.695};

    // X축(adc) 에러(sigma)
    double sigma_Q1T2[n] = {0., 0., 0., 0.}; 
    double sigma_Q4T4[n] = {0., 0., 0., 0.};
    
    // 나머지 텔레스코프 X축 에러 배열
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

    // 에러 배열들을 포인터로 묶음
    double* sigma_adc_vals[nGraph] = 
    {
        sigma_Q1T2, sigma_Q4T4, sigma_Q1T1, sigma_Q1T3,
        sigma_Q2T3, sigma_Q2T4, sigma_Q3T2, sigma_Q4T1
    };

    // 8개 그래프를 구분하기 위한 색상 추가
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

        // 마커 모양을 20번(꽉 찬 동그라미)으로 통일
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]); 
        gr->SetMarkerSize(1.5);
        gr->SetLineWidth(2);

        // 먼저 데이터를 화면에 그립니다. (피팅 실패와 무관하게 무조건 그려짐)
        gr->Draw("AP");

        gr->GetXaxis()->SetTitle("adc");
        gr->GetYaxis()->SetTitle("MeV");

        // --- 유동적인 X축 범위 설정 ---
        double min_adc = TMath::MinElement(n, adc_vals[i]);
        double max_adc = TMath::MaxElement(n, adc_vals[i]);
        
        // 피팅 및 화면 표시 구간: 최솟값 - 0.5 ~ 최댓값 + 0.5
        double fit_min = min_adc - 0.5;
        double fit_max = max_adc + 0.5;

        // X축 범위를 계산된 값으로 설정하여 데이터가 잘리지 않게 함
        gr->GetXaxis()->SetLimits(fit_min, fit_max);
        
        // Y축은 E 값이 1.241 ~ 1.537 이므로 1.1 ~ 1.7 고정 유지
        gr->GetYaxis()->SetRangeUser(1.1, 1.7);

        gr->GetXaxis()->SetTitleSize(0.045);
        gr->GetYaxis()->SetTitleSize(0.045);
        gr->GetXaxis()->SetLabelSize(0.035);
        gr->GetYaxis()->SetLabelSize(0.035);

        gr->GetXaxis()->SetTitleOffset(1.0);
        gr->GetYaxis()->SetTitleOffset(1.15);

        // 피팅 함수 지정: [0]*x^2 + [1]*x (지정된 범위 적용)
        TF1* f1 = new TF1(Form("f_%s", names[i]), "[0]*x^2+[1]*x", fit_min, fit_max);
        f1->SetLineColor(colors[i]);
        f1->SetLineWidth(2);

        // 피팅 수행. 
        // 만약 특정 데이터 배열의 분포가 너무 이상해서 피팅이 아예 실패하더라도,
        // 이미 위에서 gr->Draw("AP")를 했기 때문에 점들은 무조건 남아 있습니다.
        gr->Fit(f1, "RMI");

        // 선 위로 마커가 올라오도록 다시 그려줌
        gr->Draw("P SAME");

        double p0 = f1->GetParameter(0);
        double p1 = f1->GetParameter(1);
        double chi2   = f1->GetChisquare();
        int ndf       = f1->GetNDF();


        // 텍스트 위치 정렬
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
        
        c->SaveAs(Form("fit_%s.pdf", names[i]));

        std::cout << names[i]
                  << " : p0 (x^2) = " << p0
                  << ", p1 (x) = " << p1
                  << ", Fit Range = [" << fit_min << ", " << fit_max << "]"
                  << std::endl;
    }
}

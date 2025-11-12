// plot_pedestal_comparison.C
{
    // Detector geometry parameters
    const int nlayers = 7;
    const int nrows = 6;
    const int nmods = 6;
    const int nstrips = 32;
    const int ncn = 2;

    const int npulse = 100; //Hypothetically pulsing 100 channels!

    // Threshold parameters (fully consistent with original code)
    const float xthresh = 2.5;  //X-ray ADC threshold
    //const float mthresh = 15; //MIP Threshold?

    // Create histogram (same definition as original code)
    TH2D *hdiff = new TH2D("hdiff", "Pedestal Comparison after Pulse Correction",
                           ncn*nrows, 0, ncn*nrows,
                           nlayers*nmods, 0, nlayers*nmods);

    // Create histogram (same definition as original code)
    TH2D *hmag = new TH2D("hmag", "Magnitude of Improvement",
                           ncn*nrows, 0, ncn*nrows,
                           nlayers*nmods, 0, nlayers*nmods);

    // Create histogram for sigma before correction
    TH2D *hsig_before = new TH2D("hsig_before", "Sigma Map Pre-Correction",
                           nrows*nstrips, 0, nrows*nstrips,
                           nlayers*nmods, 0, nlayers*nmods);

    // Create histogram for sigma before correction
    TH2D *hsig_after = new TH2D("hsig_after", "Sigma Map Pre-Correction",
                           nrows*nstrips, 0, nrows*nstrips,
                           nlayers*nmods, 0, nlayers*nmods);

    // Data storage arrays [layer][row][mod][channel]
    double before_data[7][6][6][32] = {0};
    double after_data[7][6][6][32] = {0};
    int pulsed_ch[7][6][6][ncn] = {0};
    double mag_improv[7][6][6][ncn] = {0};
    double score_improv[7][6][6][ncn] = {0};

    // Read pre-correction data (corresponding to hsig_nocorr)
    std::cout << "Reading output/ped_before_CMN_removal.txt..." << std::endl;
    std::ifstream file_before("output/ped_before_CMN_removal.txt");
    std::string line;
    int count_before = 0;

    while (std::getline(file_before, line)) {
        if (line.empty() || line[0] == '#') continue;

        int l, r, m, ch;
        double mean, value;
        if (sscanf(line.c_str(), "%d %d %d %d %lf %lf", &l, &r, &m, &ch, &mean, &value) == 6) {
            if (l < nlayers && r < nrows && m < nmods && ch < nstrips) {
                before_data[l][r][m][ch] = value;
                int x_bin = r * nstrips + ch;
                int y_bin = l * nmods + m;
                hsig_before->Fill(x_bin + 0.5, y_bin + 0.5, value);  // Fill sigma map

                count_before++;
            }
        }
    }
    file_before.close();
    std::cout << "Read " << count_before << " entries from before file" << std::endl;

    // Read post-correction data (corresponding to hsig_pulsed)
    std::cout << "Reading output/ped_after_pulse_correction.txt..." << std::endl;
    std::ifstream file_after("output/ped_after_pulse_correction.txt");
    int count_after = 0;

    while (std::getline(file_after, line)) {
        if (line.empty() || line[0] == '#') continue;

        int l, r, m, ch;
        double mean, value;
        if (sscanf(line.c_str(), "%d %d %d %d %lf %lf", &l, &r, &m, &ch, &mean, &value) == 6) {
            if (l < nlayers && r < nrows && m < nmods && ch < nstrips) {
                after_data[l][r][m][ch] = value;
                int x_bin = r * nstrips + ch;
                int y_bin = l * nmods + m;
                hsig_after->Fill(x_bin + 0.5, y_bin + 0.5, value);  // Fill sigma map

                count_after++;
            }
        }
    }
    file_after.close();
    std::cout << "Read " << count_after << " entries from after file" << std::endl;

    // Read pulsed channel list
    std::cout << "Reading output/best_pulsed_channel.txt..." << std::endl;
    std::ifstream file_pulch("output/best_pulsed_channel.txt");

    while (std::getline(file_pulch, line)) {
        if (line.empty() || line[0] == '#') continue;

        int l, r, m, ch;
        int pulch0, pulch1;
        if (sscanf(line.c_str(), "%d %d %d %d %d", &l, &r, &m, &pulch0, &pulch1) == 5) {
            if (l < nlayers && r < nrows && m < nmods) {
                pulsed_ch[l][r][m][0] = pulch0;
                pulsed_ch[l][r][m][1] = pulch1;
                //cout << l << r << m << pulch0 << endl;
            }
        }
    }


    // Step 1: Set base background (same logic as original code)
    for (int l = 0; l < nlayers; l++) {
        for (int r = 0; r < nrows; r++) {
            for (int m = 0; m < nmods; m++) {
                for (int ch = 0; ch < nstrips; ch++) {
                    if (before_data[l][r][m][ch] != 0) {  // Channels with data
                        int ncn_val = ch / (nstrips / ncn);
                        int x_bin = r * ncn + ncn_val;
                        int y_bin = l * nmods + m;
                        hdiff->Fill(x_bin + 0.5, y_bin + 0.5, 0.1);  // +0.5 to fill bin center
                    }
                }
            }
        }
    }

    // Step 2: Mark channels with effective pulse correction (identical logic to original code)
    int improvement_count = 0;
    for (int l = 0; l < nlayers; l++) {
        for (int r = 0; r < nrows; r++) {
            for (int m = 0; m < nmods; m++) {
                for (int ch = 0; ch < nstrips; ch++) {
                    double before = before_data[l][r][m][ch];
                    double after = after_data[l][r][m][ch];

                    // Exactly the same conditions as original code
                    if (before > xthresh &&     // Pre-correction > threshold
                        after < xthresh &&      // Post-correction < threshold
                        after > 0) {            // Post-correction > 0

                        int ncn_val = ch / (nstrips / ncn);
                        int x_bin = r * ncn + ncn_val;
                        int y_bin = l * nmods + m;
                        hdiff->Fill(x_bin + 0.5, y_bin + 0.5, 1.0);  // Add weight 1
                        improvement_count++;
                    }
                }
            }
        }
    }

    //Now going to implement a section where we just output in a .txt file the

    std::ofstream myfile;
    myfile.open("output/FullInfo.txt");
    myfile << "X-ray Threshold = " << xthresh << endl;
    myfile << "Layer\tRow\tMod\tStrip\tPre-CMN\tPost-CMN" << endl;

    for (int l = 0; l < nlayers; l++) {
        for (int r = 0; r < nrows; r++) {
            for (int m = 0; m < nmods; m++) { //Going over all modules
                for(int n = 0; n < ncn; n++) {
                    double sum_sigma = 0;
                    int xscore = 0;
                    for(int ch = 0+n*(nstrips/ncn); ch <(n+1)*(nstrips/ncn); ch++ ){
                    if(before_data[l][r][m][ch] > xthresh && after_data[l][r][m][ch] < xthresh && after_data[l][r][m][ch] > 0 && before_data[l][r][m][ch] > 0) xscore++;
                    if(after_data[l][r][m][ch] > 0 && before_data[l][r][m][ch] > 0) sum_sigma = sum_sigma + (before_data[l][r][m][ch] - after_data[l][r][m][ch]);
                    myfile << l << "\t" << r << "\t" << m << "\t" << ch << "\t" << before_data[l][r][m][ch] << "\t" << after_data[l][r][m][ch] << endl;
                    }

                    int x_bin = r * ncn + n;
                    int y_bin = l * nmods + m;
                    hmag->Fill(x_bin + 0.5, y_bin + 0.5, sum_sigma);  // +0.5 to fill bin center
                    mag_improv[l][r][m][n] = sum_sigma;
                    score_improv[l][r][m][n] = xscore;

                    myfile << "Pulsed channel = " << pulsed_ch[l][r][m][n] << endl;
                    myfile << "X-ray score = " << xscore << endl;
                    myfile << "Total sigma improvement of ncn " << n << " = " <<  sum_sigma << endl << endl;
                }

            }
        }
    }

    myfile.close();

    //Ranked pulse channel
    //Iterate over the ncn's. Use score and magnitude to determine which channels to pulse.
    //Start by just making sure you can organize by score and
    std::ofstream rankfile;
    rankfile.open("output/RankedChannels.txt");
    rankfile << "X-ray Threshold = " << xthresh << endl;
    rankfile << "Layer\tRow\tMod\tStrip\tX-score\tMag Improv" << endl;



    rankfile.close();


    std::cout << "Found " << improvement_count << " channels with significant improvement" << std::endl;


    ////////////////////////////////////////////
    /// Plotting section
    ///
    ////////////////////////////////////////////

    // Create grid lines (same as original code)
    std::vector<TLine*> lines;
    for (int i = 0; i <= nrows; i++) {
        lines.push_back(new TLine(i * ncn, 0, i * ncn, nlayers * nmods));
    }
    for (int i = 0; i <= nlayers; i++) {
        lines.push_back(new TLine(0, i * nmods, nrows * ncn, i * nmods));
    }

    // Plotting (same style as original code)
    TCanvas *c1 = new TCanvas("c1", "Pedestal Comparison", 900, 1100);
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);

    hdiff->SetStats(0);
    hdiff->GetXaxis()->SetTitle("row(0-5)*2 + ncn_group(0-1)");
    hdiff->GetXaxis()->SetNdivisions(12, 0, 0, kFALSE);  // 13 ticks

    hdiff->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
    hdiff->GetZaxis()->SetTitle("Improvement Indicator");
    hdiff->GetZaxis()->SetTitleOffset(1.8);
    hdiff->SetMinimum(0);
    hdiff->SetMaximum(nstrips/ncn);  // Maximum value 16

    // Set colors (same as original code)
    gStyle->SetPalette(55);  // kRainBow palette

    hdiff->Draw("COLZ");

    // Draw grid lines
    for (auto line : lines) {
        line->SetLineColor(kBlack);
        line->SetLineWidth(1);
        line->Draw();
    }

    c1->SaveAs("output/pedestal_comparison.png");
    c1->SaveAs("output/pedestal_comparison.pdf");

    std::cout << "Saved pedestal_comparison.png and .pdf" << std::endl;

    std::vector<TLine*> lines2;
    // Create grid lines (same as original code)
    for (int i = 0; i <= nrows; i++) {
        lines2.push_back(new TLine(i * nstrips, 0, i * nstrips, nlayers * nmods));
    }
    for (int i = 0; i <= nlayers; i++) {
        lines2.push_back(new TLine(0, i * nmods, nrows * nstrips, i * nmods));
    }

    // Plotting (same style as original code)
    TCanvas *c2 = new TCanvas("c2", "Sigma Map", 900, 1100);
    c2->SetLeftMargin(0.1);
    c2->SetRightMargin(0.15);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.1);

    hsig_before->SetStats(0);
    hsig_before->GetXaxis()->SetTitle("row(0-5) + strip(0-31)");
    hsig_before->GetXaxis()->SetNdivisions(12, 0, 0, kFALSE);  // 13 ticks

    hsig_before->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
    hsig_before->GetZaxis()->SetTitle("Sigma");
    hsig_before->GetZaxis()->SetTitleOffset(1.8);
    hsig_before->SetMinimum(0);
    hsig_before->SetMaximum(100);

    gStyle->SetPalette();  // kRainBow palette

    hsig_before->Draw("COLZ");

    // Draw grid lines
    for (auto line : lines2) {
        line->SetLineColor(kBlack);
        line->SetLineWidth(1);
        line->Draw();
    }

    c2->SaveAs("output/sigma_map_before.png");

    // Plotting (same style as original code)
    TCanvas *c3 = new TCanvas("c3", "Post-CMN Sigma Map", 900, 1100);
    c3->SetLeftMargin(0.1);
    c3->SetRightMargin(0.15);
    c3->SetTopMargin(0.05);
    c3->SetBottomMargin(0.1);

    hsig_after->SetStats(0);
    hsig_after->GetXaxis()->SetTitle("row(0-5) + strip(0-31)");
    hsig_after->GetXaxis()->SetNdivisions(12, 0, 0, kFALSE);  // 13 ticks

    hsig_after->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
    hsig_after->GetZaxis()->SetTitle("Sigma");
    hsig_after->GetZaxis()->SetTitleOffset(1.8);
    hsig_after->SetMinimum(0);
    hsig_after->SetMaximum(100);

    // Set colors (same as original code)
    //gStyle->SetPalette(55);  // kRainBow palette

    hsig_after->Draw("COLZ");

    // Draw grid lines
    for (auto line : lines2) {
        line->SetLineColor(kBlack);
        line->SetLineWidth(1);
        line->Draw();
    }

    c3->SaveAs("output/sigma_map_after.png");

    //std::cout << "Saved pedestal_comparison.png and .pdf" << std::endl;

    // Plotting (same style as original code)
    TCanvas *c4 = new TCanvas("c4", "Post-CMN Sigma Map", 900, 1100);
    c4->SetLeftMargin(0.1);
    c4->SetRightMargin(0.15);
    c4->SetTopMargin(0.05);
    c4->SetBottomMargin(0.1);

    hmag->SetStats(0);
    hmag->GetXaxis()->SetTitle("row(0-5) + strip(0-31)");
    hmag->GetXaxis()->SetNdivisions(12, 0, 0, kFALSE);  // 13 ticks

    hmag->GetYaxis()->SetTitle("layer(0-6)*6 + mod(0-5)");
    hmag->GetZaxis()->SetTitle("Sigma");
    hmag->GetZaxis()->SetTitleOffset(1.8);
    hmag->SetMinimum(-10);
    hmag->SetMaximum(300);

    //Set colors (same as original code)
    gStyle->SetPalette(55);  // kRainBow palette

    hmag->Draw("COLZ");

    // Draw grid lines
    for (auto line : lines) {
        line->SetLineColor(kBlack);
        line->SetLineWidth(1);
        line->Draw();
    }

    c4->SaveAs("output/sigma_magnitude_improve.png");

    // Clean up memory
    for (auto line : lines) delete line;
}

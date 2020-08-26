float boltz(float x){
    float kt = 155.0;
    return sqrt(x)*expf(-x/kt);
}
float model(float x){
    float v0 = 0,T = 4.7e-7,beta = 1.39,alpha = 0.36,c = 1,m = 1e-6,h = 4.13e-15;
    return powf(x*h/(m*T),alpha)*expf(-powf(((x)*h/(m*T)),beta));
}

int main(void){
    TF1 *f1 = new TF1("f1","model(x)",0,600);
    TF1 *f2 = new TF1("f2","boltz(x)",0,600);

    TH1F *h1 = new TH1F("h1","Histogram",200,0,600);
    h1->FillRandom("f1",1e5);
    h1->SetLineColor(2);
    h1->GetXaxis()->SetTitle("Frequency (1/s)");
    h1->GetYaxis()->SetTitle("Amplitude");
    h1->Scale(1/h1->Integral("width"));
    h1->Draw("hist");

    TH1F *h2 = new TH1F("h2","Histogram2",200,0,600);
    h2->FillRandom("f2",1e5);
    h2->SetLineColor(4);
    h2->Scale(1/h2->Integral("width"));
    h2->Draw("same&hist");

    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry("h1","new model");
    leg->AddEntry("h2","maxwell-boltzman");
    leg->Draw();
}
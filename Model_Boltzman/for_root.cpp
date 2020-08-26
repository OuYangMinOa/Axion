float boltz(float E){
    float normalize = 0.0017;
    return sqrt(E)*expf(-E/(0.013*1.46))/normalize;
    ;}
float model(float E){
    float normalize = 82;
    float v0 = 0,T = 4.7e-7,beta = 1.39,alpha = 0.36,c = 1,m = 1e-6,h = 4.13e-15;
    return powf((E)*h/(m*T),alpha)*expf(-powf(((E)*h/(m*T)),beta))/normalize;
    ;}


TF1 *f1 = new TF1("f1","model(x)",0,600);
f1->SetTitle("Curve");
f1->GetXaxis()->SetTitle("Frequency (1/s)");
f1->GetYaxis()->SetTitle("Amplitude");
f1->Draw();

TF1 *f2 = new TF1("f2","boltz(x)",0,600);
f2->SetLineColor(4);
f2->Draw("same");

TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry("f1","new model");
leg->AddEntry("f2","maxwell-boltzman");
leg->Draw()



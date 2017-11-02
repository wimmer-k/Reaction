#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TKey.h"
#include "TMath.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TColor.h"

#include "CommandLineInterface.hh"
#include "Kinematics.hh"
#include "Nucleus.hh"
#ifndef rad2deg
#define rad2deg                       180./(TMath::Pi())
#endif
#ifndef deg2rad
#define deg2rad                       (TMath::Pi())/180.
#endif

using namespace TMath;
using namespace std;
int main(int argc, char* argv[]){
  vector<int> projectile;
  vector<int> target;
  vector<double> transfer;
  vector<double> coulex;
  int chargeex=0;
  int pickup=0;
  int protont=0;
  int fusion=0;
  double ebeam;
  double steps=0;
  bool rutherford = false;
  bool cm2lab = false;
  bool iso = false;
  bool test = false;
  char* rootFile = NULL;
  char* cmFile = NULL;
  char* psFile = NULL;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-p", "projectile N | Z", &projectile);  
  interface->Add("-t", "target N | Z", &target);  
  interface->Add("-tr", "neutrons to be transferred (from target to beam) | excitation energy [MeV] | nr of steps", &transfer);    
  interface->Add("-pt", "protons to be transferred (from target to beam)", &protont);    
  interface->Add("-pi", "neutrons to be transferred (from beam to target)", &pickup);    
  interface->Add("-ce", "charge exchange (from target to beam)", &chargeex);    
  interface->Add("-co", "coulex states", &coulex);    
  interface->Add("-f", "fusion", &fusion);  
  interface->Add("-e", "energy [MeV]", &ebeam);  
  interface->Add("-s", "step size (opt)", &steps);
  interface->Add("-ru", "Rutherford cross sections", &rutherford);
  interface->Add("-is", "generate cm distribution for isotropic lab", &iso);
  interface->Add("-cf", "file with splines of cm cross section", &cmFile);
  interface->Add("-o", "output root file", &rootFile);
  interface->Add("-te", "test", &test);
  interface->Add("-ps", "ps file for kinematics", &psFile);

  interface->CheckFlags(argc, argv);
  TFile* outfile = NULL;
  if(test == true){
    cout<<"testing stuff"<< endl;
  }
  else if(!rutherford && cmFile==NULL && iso == false)
    cout<<"just kinematics no cross section"<< endl;
  else if(rootFile == NULL){
    cerr<<"You have to provide a output file!"<<endl;
      exit(1);
  }
  if(psFile == NULL)
    psFile = (char*)"kinematics.ps";
  if(rootFile != NULL){
    outfile = new TFile(rootFile,"recreate");
    if(outfile->IsZombie())
      return 4;
  }
  if(cmFile!=NULL)
    cm2lab=true;


  TColor *orange = gROOT->GetColor(5);
  orange->SetRGB(1.0, 0.612, 0.002); 
  
  TColor *green = gROOT->GetColor(3);
  green->SetRGB(0.15, 0.7, 0.15);

  if(steps==0)
    steps=1;

  Nucleus *proj;
  Nucleus *targ;
  Nucleus *reco;
  Nucleus *ejec;

  int numb =0;

  if(projectile.size() == 2){
    cout << "projectile Z " << projectile[1] << " N " << projectile[0] << endl;
    proj = new Nucleus(projectile[1],projectile[0], (char*)MASSFILE);
  }
  else{
    cerr<<"no or incorrect Projectile provided!";
    for(UShort_t i=0; i<projectile.size(); i++){
      cerr<<projectile[i]<<" ";
    }
    cerr<<endl;
    exit(1);
  }

  if(target.size() == 2){
    cout << "target Z " << target[1] << " N " << target[0] << endl;
    targ = new Nucleus(target[1],target[0], (char*)MASSFILE);
  }
  else{
    cerr<<"flag -p provided but not two arguments following: ";
    for(UShort_t i=0; i<target.size(); i++){
      cerr<<target[i]<<" ";
    }
    exit(1);
  }

  int trexmin = 7;
  int trexmax = 74;
  int oldmin = 16;
  int oldmax = 54;
  
  Kinematics *boom;
  Kinematics *boomar[17];
  int trst =1;
  cout << "beam energy " << ebeam << endl;
  cout << "proj " << proj->GetMass() << " targ " << targ->GetMass() << endl;
  if(transfer.size() > 0&&transfer.size() <4){
    cout << "transfer reaction deltaN = " << transfer[0] <<endl;
    cout << "----------------------------"<<endl;
    reco = new Nucleus(target[1],target[0]-(int)transfer[0], (char*)MASSFILE);
    cout << "recoil Z " << target[1] << "\tN " << target[0]-transfer[0] <<"\t" << reco->GetMass()<< endl;
    ejec = new Nucleus(projectile[1],projectile[0]+(int)transfer[0], (char*)MASSFILE);
    cout << "ejectile Z " << projectile[1] << "\tN " << projectile[0]+transfer[0] <<"\t" << ejec->GetMass()<< endl;
    if(transfer.size() == 3){
      cout << transfer[2]<< " steps to state with excitation energy = " << transfer[1] <<endl;
      for(int k=0;k<transfer[2]&&transfer[2]<17;k++){
	boomar[k] = new Kinematics(proj, targ, reco, ejec, ebeam, transfer[1]*k/(transfer[2]-1));
      }
      trst =(int)transfer[2];
      numb = trst;
    }
    if(transfer.size() == 2){
      cout << "in state with excitation energy = " << transfer[1] <<endl;
      boom = new Kinematics(proj, targ, reco, ejec, ebeam, transfer[1]);
      boomar[0] = new Kinematics(proj, targ, reco, ejec, ebeam, transfer[1]);
      trst =1;
    }
    else
      boom = new Kinematics(proj, targ, reco, ejec, ebeam, 0);

    //cout << "ejec " << ejec->GetMass() << " reco " << reco->GetMass() << endl;

  }
  else if(pickup!=0){
    cout << "transfer reaction deltaN = " << -pickup <<endl;
    cout << "----------------------------"<<endl;
    reco = new Nucleus(target[1],target[0]+(int)pickup, (char*)MASSFILE);
    cout << "recoil Z " << target[1] << " N " << target[0]+pickup << endl;
    ejec = new Nucleus(projectile[1],projectile[0]-(int)pickup, (char*)MASSFILE);
    cout << "ejectile Z " << projectile[1] << " N " << projectile[0]-pickup << endl;
    cout << "ejec Z " << ejec->GetZ() << " N " << ejec->GetN() << " reco Z " << reco->GetZ() << " N " << reco->GetZ() << endl;

    boom = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
  }
  else if(protont!=0){
    cout << "proton transfer reaction deltaP = " << protont <<endl;
    cout << "----------------------------"<<endl;
    reco = new Nucleus(target[1]-protont,target[0], (char*)MASSFILE);
    cout << "recoil Z " << target[1]-protont << " N " << target[0] << endl;
    ejec = new Nucleus(projectile[1]+protont,projectile[0], (char*)MASSFILE);
    cout << "ejectile Z " << projectile[1]+protont << " N " << projectile[0] << endl;
    //cout << "ejec " << ejec->GetMass() << " reco " << reco->GetMass() << endl;
 
    boom = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
  }
  else if(chargeex!=0){
    cout << "charge exchange deltaZ = " << chargeex <<endl;
    cout << "----------------------------"<<endl;
    reco = new Nucleus(target[1]+chargeex,target[0]-chargeex, (char*)MASSFILE);
    cout << "recoil Z " << target[1]+chargeex << " N " << target[0]-chargeex << endl;
    ejec = new Nucleus(projectile[1]-chargeex,projectile[0]+chargeex, (char*)MASSFILE);
    cout << "ejectile Z " << projectile[1] -chargeex<< " N " << projectile[0]+chargeex << endl;
    //cout << "ejec " << ejec->GetMass() << " reco " << reco->GetMass() << endl;

    boom = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
    

  }
  else if(fusion!=0){
    cout << "fusion reaction " <<endl;
    cout << "----------------------------"<<endl;
    reco = new Nucleus(target[1]+projectile[1],target[0]+projectile[0], (char*)MASSFILE);
    cout << "recoil Z " << target[1]+projectile[1] << " N " << target[0]+projectile[0] << endl;
    ejec = new Nucleus(0,0, (char*)MASSFILE);
    cout << "ejectile Z " << 0 << " N " << 0 << endl;
    cout << "ejec " << ejec->GetMass() << " reco " << reco->GetMass() << endl;

    boom = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
  }
  else if(coulex.size() > 0){
    cout << "coulex " <<endl;
    cout << "----------------------------"<<endl;
    ejec = new Nucleus(projectile[1],projectile[0], (char*)MASSFILE);
    reco = new Nucleus(target[1],target[0], (char*)MASSFILE);
    trst = coulex.size();
    for(UShort_t k=0;k<coulex.size()&&coulex.size()<17;k++){
      boomar[k] = new Kinematics(proj, targ, targ, proj, ebeam, coulex[k]);
    }
    boom = new Kinematics(proj, targ, ebeam);
    numb =2;

  }
  else{
    cout << "elastic scattering " <<endl;
    cout << "----------------------------"<<endl;
    ejec = new Nucleus(projectile[1],projectile[0], (char*)MASSFILE);
    reco = new Nucleus(target[1],target[0], (char*)MASSFILE);
	
    boom = new Kinematics(proj, targ, ebeam);
  }
  cout << "Qvalue " << boom->GetQValue() << endl;
  cout << "----------------------------"<<endl;
  cout << "norm kin beam energy " <<endl;
  cout << boom->NormalkinEnergy() << endl;
  cout << "----------------------------"<<endl;
  if(fusion){
    cout << "barrier " <<endl;
    cout << proj->GetZ()*targ->GetZ()/137.*197./(proj->GetRadius()+targ->GetRadius()+3.54) << endl;
    //cout << proj->GetZ()*targ->GetZ()/137.*197./(1.2*pow(proj->GetA(),1./3.)+1.2*pow(targ->GetA(),1./3.)) << endl;
    //cout << "proj " << proj->GetRadius() << " my " << 1.2*pow(proj->GetA(),1./3.) << endl;
    //cout << "targ " << targ->GetRadius() << " my " << 1.2*pow(targ->GetA(),1./3.) << endl;
    //cout << 1./(51.58/(proj->GetZ()*targ->GetZ()/137.*197.)) << endl;
    cout << "----------------------------"<<endl;
  }
  cout << "cm energy " <<endl;
  cout << boom->GetCmEnergy() << endl;
  cout << "----------------------------"<<endl;
  cout << "maximum scattering angle: " << boom->GetMaxAngle(3)*180./PI << endl;
  cout << "lab energy there:         " << boom->ELab(boom->GetMaxAngle(3), 3) << endl;
  cout << "beta " <<endl;
  cout << boom->GetBetacm() << endl;
  for(int i=0;i<4;i++)
    cout << i << "\t"<< boom->GetBetacm(i)<< "\t"<< boom->GetV(i) << "\t" << boom->GetEcm(i) << "\t" << boom->GetTcm(i) << endl;

  //cout << "closest approach = " << proj->GetZ()<<"*"<<targ->GetZ()<<"*1.44/2./"<<boom->GetTlab(0) << " = " << proj->GetZ()*targ->GetZ()*1.44/2./boom->GetTlab(0) << endl;
  if(!rutherford && !cm2lab && !iso &&!test){
    cout << "kinematics " << endl;

    double maxelab_t =0.;

    vector<TSpline3*> EvsTh_lab_p;
    vector<TSpline3*> EvsTh_lab_t;
    vector<TSpline3*> cmvslab_t;
    vector<TSpline3*> cmvslab_p;
    vector<int> n_p;
    vector<int> n_t;
 
    EvsTh_lab_p.resize(trst);
    EvsTh_lab_t.resize(trst);
    cmvslab_t.resize(trst);
    cmvslab_p.resize(trst);
    n_p.resize(trst);
    n_t.resize(trst);

    vector<vector<double> > thetalab_p;
    vector<vector<double> > energylab_p;
    vector<vector<double> > thetalab_t;
    vector<vector<double> > energylab_t;
    vector<vector<double> > cm_p;
    vector<vector<double> > lab_p;
    vector<vector<double> > cm_t;
    vector<vector<double> > lab_t;

    thetalab_p.resize(trst);
    energylab_p.resize(trst);
    thetalab_t.resize(trst);
    energylab_t.resize(trst);
    cm_p.resize(trst);
    lab_p.resize(trst);
    cm_t.resize(trst);
    lab_t.resize(trst);
    for(UShort_t s=0;s<trst;s++){
      thetalab_p[s].resize((int)(180./steps));
      energylab_p[s].resize((int)(180./steps));
      thetalab_t[s].resize((int)(180./steps));
      energylab_t[s].resize((int)(180./steps));
      cm_p[s].resize((int)(180./steps));
      lab_p[s].resize((int)(180./steps));
      cm_t[s].resize((int)(180./steps));
      lab_t[s].resize((int)(180./steps));
    }

    if(trst==1)
      boomar[0]=boom;
    for(int k=0;k<trst&&trst<17;k++){
      EvsTh_lab_p[k] = boomar[k]->Evslab(1, 179, steps,3);
      EvsTh_lab_t[k] = boomar[k]->Evslab(1, 179, steps,2);
      n_p[k]=0;
      n_t[k]=0;
      if(rootFile != NULL){
	//cout << "writing " << Form("EvsTh_lab_p_%d",k) << endl;
	EvsTh_lab_p[k]->Write(Form("EvsTh_lab_p_%d",k),TObject::kOverwrite);
	//cout << "writing " << Form("EvsTh_lab_t_%d",k) << endl;
	EvsTh_lab_t[k]->Write(Form("EvsTh_lab_t_%d",k),TObject::kOverwrite);
      }
      for(int i=0;i<EvsTh_lab_p[k]->GetNp();i++){
	EvsTh_lab_p[k]->GetKnot(i, thetalab_p[k][i], energylab_p[k][i]); 
	energylab_p[k][i]/=1000;
	if((thetalab_p[k][i]>90.1||thetalab_p[k][i]<89.9)&&thetalab_p[k][i]>0&&energylab_p[k][i]>0){
	  thetalab_p[k][n_p[k]]=thetalab_p[k][i];
	  energylab_p[k][n_p[k]]=energylab_p[k][i];
	  n_p[k]++;
	}
	
	//energylab_p[k][i]/=1000.;
	//cout << " theta " << thetalab_p[k][i] << " cs " << energylab_p[k][i] << endl;      
      }
      for(int i=0;i<EvsTh_lab_t[k]->GetNp();i++){
	EvsTh_lab_t[k]->GetKnot(i, thetalab_t[k][i], energylab_t[k][i]);      
	//energylab_t[k][i]/=1000.;
	//cout << " theta " << thetalab_t[k][i] << " cs " << energylab_t[k][i] << endl;      

	energylab_t[k][i]/=1000;
	if(thetalab_t[k][i]>0&&energylab_t[k][i]>0){
	  thetalab_t[k][n_t[k]]=thetalab_t[k][i];
	  energylab_t[k][n_t[k]]=energylab_t[k][i];
	  n_t[k]++;
	}
	if(energylab_t[k][i]>maxelab_t)
	  maxelab_t=energylab_t[k][i];
      }
      cmvslab_t[k] = boomar[k]->cmvslab(1, 179, steps,2);
      cmvslab_p[k] = boomar[k]->cmvslab(1, 179, steps,3);
      for(int i=0;i<cmvslab_p[k]->GetNp();i++){
	cmvslab_p[k]->GetKnot(i, lab_p[k][i], cm_p[k][i]);  
	//cm_p[k][i] = 180-cm_p[k][i];
      }
      for(int i=0;i<cmvslab_t[k]->GetNp();i++){
	cmvslab_t[k]->GetKnot(i, lab_t[k][i], cm_t[k][i]); 
	cm_t[k][i] = 180-cm_t[k][i];
      }
    }
    TSpline3* EvsTh_cm_p;
    TSpline3* EvsTh_cm_t;
    EvsTh_cm_p = boom->Evscm(1, 179, steps,3);
    EvsTh_cm_t = boom->Evscm(1, 179, steps,2);
    double maxecm_t =0.;
    double maxecm_p =0.;
    vector<double> thetacm_p;
    vector<double> energycm_p;
    thetacm_p.resize(EvsTh_cm_p->GetNp());
    energycm_p.resize(EvsTh_cm_p->GetNp());
    int ncm_p=0;
    int ncm_t=0;
    for(int i=0;i<EvsTh_cm_p->GetNp();i++){
      EvsTh_cm_p->GetKnot(i, thetacm_p[i], energycm_p[i]);      
      energycm_p[i]/=1000.;
      if(thetacm_p[i]<180&&thetacm_p[i]>0&&energycm_p[i]>0){
	thetacm_p[ncm_p]=thetacm_p[i];
	energycm_p[ncm_p]=energycm_p[i];
	ncm_p++;
      }
      if(energycm_p[i]>maxecm_p)
	maxecm_p=energycm_p[i];
      //cout << " theta " << x[i] << " cs " << y[i] << endl;      
    }
    vector<double> thetacm_t;
    vector<double> energycm_t;
    thetacm_t.resize(EvsTh_cm_t->GetNp());
    energycm_t.resize(EvsTh_cm_t->GetNp());
    for(int i=0;i<EvsTh_cm_t->GetNp();i++){
      EvsTh_cm_t->GetKnot(i, thetacm_t[i], energycm_t[i]);      
      energycm_t[i]/=1000.;
      if(thetacm_t[i]<180&&thetacm_t[i]>0&&energycm_t[i]>0){
	thetacm_t[ncm_t]=thetacm_t[i];
	energycm_t[ncm_t]=energycm_t[i];
	ncm_t++;
      }
      if(energycm_t[i]>maxecm_t)
	maxecm_t=energycm_t[i];
      //cout << " theta " << x[i] << " cs " << y[i] << endl;      
    }

    

    TCanvas* div4 = new TCanvas("div4","div4", 800, 600);
    TPaveLabel *pl = new TPaveLabel(0.02,0.95,0.98,0.99,Form("%s(%s,%s)%s at %.2f MeV",targ->GetSymbol(),proj->GetSymbol(),ejec->GetSymbol(),reco->GetSymbol(),ebeam),"br");
    pl->SetShadowColor(0);
    pl->SetFillColor(0);
    pl->SetBorderSize(0);
    pl->Draw();

    TPad *p[4];
    p[0] = new TPad(Form("p%d",0),Form("p%d",0),0.02,0.47,0.49,0.94);
    p[1] = new TPad(Form("p%d",1),Form("p%d",1),0.51,0.47,0.98,0.94);
    p[2] = new TPad(Form("p%d",2),Form("p%d",2),0.02,0.00,0.49,0.47);
    p[3] = new TPad(Form("p%d",3),Form("p%d",3),0.51,0.00,0.98,0.47);
    
    for(int i=0;i<4;i++){
      p[i]->Draw();
      p[i]->SetLeftMargin(0.12);
      p[i]->SetRightMargin(0.05);
      p[i]->SetBottomMargin(0.11);
      p[i]->SetTopMargin(0.07);
    }
    TAxis *xa, *ya;  
    p[0]->cd();

    double dx[2] = {0,180};
    double dy[2] = {0,ebeam*1200};
    TGraph *dg = new TGraph(2,dx,dy);
    dg->SetTitle("");
    xa = dg->GetXaxis();
    xa->SetRangeUser(0.,180.);
    xa->SetTitle("#vartheta_{lab} [#circ]");
    xa->SetTitleOffset(1.0);
    xa->SetTitleSize(0.05);
    xa->SetLabelOffset(0.008);
    xa->SetLabelColor(1);
    xa->SetLabelSize(0.05);
    xa->SetTickLength(0.03);
    ya = dg->GetYaxis();
    ya->SetRangeUser(0., ebeam*1.2);
    ya->SetTitle("Energy [MeV]");
    ya->SetTitleOffset(1.2);
    ya->SetTitleSize(0.05);
    ya->SetLabelOffset(0.005);
    ya->SetLabelColor(1);
    ya->SetLabelSize(0.05);
    ya->SetTickLength(0.03);
    dg->SetMarkerSize(0);
    dg->Draw("AP");
    TGraph *recog = new TGraph(n_t[0],&thetalab_t[0][0],&energylab_t[0][0]);
    recog->SetLineColor(2);
    recog->Draw("sameL");
    TGraph *ejecg = new TGraph(n_p[0],&thetalab_p[0][0],&energylab_p[0][0]);
    ejecg->SetLineColor(3);
    ejecg->Draw("sameL");
    TLegend *leg_tr = new TLegend(0.5, 0.75, 0.8, 0.9);

    leg_tr->AddEntry(recog, Form(" recoil %s",reco->GetSymbol()), "L");
    leg_tr->AddEntry(ejecg, Form(" ejectile %s",ejec->GetSymbol()), "L");
    leg_tr->SetBorderSize(0);
    leg_tr->SetFillColor(0);
    leg_tr->SetFillStyle(0);
    leg_tr->SetTextSize(0.035);
    leg_tr->Draw();
    if(rootFile != NULL){
      recog->Write("recoil",TObject::kOverwrite);
      ejecg->Write("ejectile",TObject::kOverwrite);
    }

    p[1]->cd();
    vector <TGraph*> recog2;
    recog2.resize(trst);
    recog2[0] = new TGraph(n_t[0],&thetalab_t[0][0],&energylab_t[0][0]);
    recog2[0]->SetTitle("");
    xa = recog2[0]->GetXaxis();
    xa->SetRangeUser(0.,180.);
    xa->SetTitle("#vartheta_{lab} [#circ]");
    xa->SetTitleOffset(1.0);
    xa->SetTitleSize(0.05);
    xa->SetLabelOffset(0.008);
    xa->SetLabelColor(1);
    xa->SetLabelSize(0.05);
    xa->SetTickLength(0.03);
    ya = recog2[0]->GetYaxis();
    ya->SetRangeUser(0., maxelab_t*1.2);
    ya->SetTitle("Energy [MeV]");
    ya->SetTitleOffset(1.2);
    ya->SetTitleSize(0.05);
    ya->SetLabelOffset(0.005);
    ya->SetLabelColor(1);
    ya->SetLabelSize(0.05);
    ya->SetTickLength(0.03);
    recog2[0]->SetLineColor(2);
    recog2[0]->Draw("AL");
    if(rootFile != NULL)
      recog2[0]->Write("EvsTh_lab_0",TObject::kOverwrite);
    for(int k=1; k<trst;k++){
      recog2[k] = new TGraph(n_t[k],&thetalab_t[k][0],&energylab_t[k][0]);
      if(rootFile != NULL)
	recog2[k]->Write(Form("EvsTh_lab_%d",k),TObject::kOverwrite);
      recog2[k]->SetLineColor(2);
      recog2[k]->Draw("L");
    }
    leg_tr = new TLegend(0.5, 0.75, 0.8, 0.9);

    if(trst>1)
      leg_tr->AddEntry(recog2[0], Form("E_{eexc} %.3f MeV",0.0), "");
    else if(transfer.size()==0)
      leg_tr->AddEntry(recog2[0],"elastic", "");
    else if(transfer.size()==1)
      leg_tr->AddEntry(recog2[0], Form("E_{eexc} %.3f MeV",transfer[1]), "");
    for(int k=1; k<trst;k++){
      if(transfer.size()>0)
	leg_tr->AddEntry(recog2[k], Form("E_{eexc} %.3f MeV",transfer[1]*k/(transfer[2]-1)), "");
      else if(coulex.size()>0)
	leg_tr->AddEntry(recog2[k], Form("E_{eexc} %.3f MeV",coulex[k]), "");
    }
    leg_tr->SetBorderSize(0);
    leg_tr->SetFillColor(0);
    leg_tr->SetFillStyle(0);
    leg_tr->SetTextSize(0.035);
    leg_tr->Draw();



    p[2]->cd();
    vector<TGraph*> rec_labcm;
    vector<TGraph*> eje_labcm;
    rec_labcm.resize(trst);
    eje_labcm.resize(trst);
    rec_labcm[0] = new TGraph(cmvslab_t[0]->GetNp(),&cm_t[0][0],&lab_t[0][0]);
    rec_labcm[0]->SetTitle("");
    xa = rec_labcm[0]->GetXaxis();
    xa->SetRangeUser(0.,180.);
    xa->SetTitle("#vartheta_{cm} [#circ]");
    xa->SetTitleOffset(1.0);
    xa->SetTitleSize(0.05);
    xa->SetLabelOffset(0.008);
    xa->SetLabelColor(1);
    xa->SetLabelSize(0.05);
    xa->SetTickLength(0.03);
    ya = rec_labcm[0]->GetYaxis();
    ya->SetRangeUser(0., 180.);
    ya->SetTitle("#vartheta_{lab} [#circ]");
    ya->SetTitleOffset(1.2);
    ya->SetTitleSize(0.05);
    ya->SetLabelOffset(0.005);
    ya->SetLabelColor(1);
    ya->SetLabelSize(0.05);
    ya->SetTickLength(0.03);
    rec_labcm[0]->SetLineColor(2);
    rec_labcm[0]->Draw("AL");
    eje_labcm[0] = new TGraph(cmvslab_p[0]->GetNp(),&cm_p[0][0],&lab_p[0][0]);
    eje_labcm[0]->SetLineColor(3);
    eje_labcm[0]->Draw("sameL");

    for(int k=1; k<trst;k++){
      rec_labcm[k] = new TGraph(cmvslab_t[k]->GetNp(),&cm_t[k][0],&lab_t[k][0]);
      rec_labcm[k]->SetLineColor(2);
      rec_labcm[k]->Draw("L");
      eje_labcm[k] = new TGraph(cmvslab_p[k]->GetNp(),&cm_t[k][0],&lab_p[k][0]);
      eje_labcm[k]->SetLineColor(3);
      eje_labcm[k]->Draw("L");
    }

    leg_tr = new TLegend(0.5, 0.75, 0.8, 0.9);
    leg_tr->AddEntry(rec_labcm[0], Form("recoil %s",reco->GetSymbol()), "L");
    leg_tr->AddEntry(eje_labcm[0], Form("ejectile %s",ejec->GetSymbol()), "L");
    leg_tr->SetBorderSize(0);
    leg_tr->SetFillColor(0);
    leg_tr->SetFillStyle(0);
    leg_tr->SetTextSize(0.035);
    leg_tr->Draw();

    p[3]->cd();
    recog = new TGraph(ncm_t,&thetacm_p[0],&energycm_t[0]);
    recog->SetTitle("");
    xa = recog->GetXaxis();
    xa->SetRangeUser(0.,180.);
    xa->SetTitle("#vartheta_{cm} [#circ]");
    xa->SetTitleOffset(1.0);
    xa->SetTitleSize(0.05);
    xa->SetLabelOffset(0.008);
    xa->SetLabelColor(1);
    xa->SetLabelSize(0.05);
    xa->SetTickLength(0.03);
    ya = recog->GetYaxis();
    ya->SetRangeUser(0., maxecm_t);
    ya->SetTitle("Energy [MeV]");
    ya->SetTitleOffset(1.2);
    ya->SetTitleSize(0.05);
    ya->SetLabelOffset(0.005);
    ya->SetLabelColor(1);
    ya->SetLabelSize(0.05);
    ya->SetTickLength(0.03);
    recog->SetLineColor(2);
    recog->Draw("AL");
    ejecg = new TGraph(ncm_p,&thetacm_p[0],&energycm_p[0]);
    ejecg->SetLineColor(3);
    ejecg->Draw("sameL");
    leg_tr = new TLegend(0.5, 0.75, 0.8, 0.9);

    leg_tr->AddEntry(recog, Form("recoil %s",reco->GetSymbol()), "L");
    leg_tr->AddEntry(ejecg, Form("ejectile %s",ejec->GetSymbol()), "L");
    leg_tr->SetBorderSize(0);
    leg_tr->SetFillColor(0);
    leg_tr->SetFillStyle(0);
    leg_tr->SetTextSize(0.035);
    leg_tr->Draw();
      
    div4->Print(Form("%s",psFile));
    
  }
  else if(rutherford){
    cout << "rutherford " << endl;
    TSpline3* splinecm;
    TSpline3* splinelabt;
    TSpline3* splinelabp;
    cout <<"creating splinecm" << endl;
    splinecm = boom->Ruthvscm(1, 179, steps); 
    cout << splinecm->GetNp() << endl;
    splinecm->SetLineWidth(1);
    splinecm->SetName("cm");
    splinecm->SetLineWidth(1);
    splinecm->Write("",TObject::kOverwrite);
    ofstream txtfile("lab_ruther.dat");    
    cout <<"creating splinelab projectile" << endl;
    splinelabp = boom->Ruthvslab(1, 179, steps,3); 
    if(splinelabp->GetXmax() < splinelabp->GetXmin() ){
      cout << " inverting spline xmax < xmin!" << endl;
      cout << splinelabp->GetNp() << endl;
      vector<double> x;
      vector<double> y;
      x.resize(splinelabp->GetNp());
      y.resize(splinelabp->GetNp());
      for(int i=0;i<splinelabp->GetNp();i++){
	splinelabp->GetKnot(splinelabp->GetNp()-i-1, x[i], y[i]);      
	//cout << " theta " << x[i] << " cs " << y[i] << endl;      
      }
      TGraph* graph = new TGraph(splinelabp->GetNp(), &x[0], &y[0]);
      TSpline3* spline = new TSpline3("pro",graph);
      spline->SetLineColor(3);
      spline->SetName("pro");
      spline->SetLineWidth(1);
      spline->Write("",TObject::kOverwrite);
      delete spline;
    }
    else{
      cout << " not inverting spline xmax > xmin!" << endl;
      splinelabp->SetLineColor(3);
      splinelabp->SetName("pro");
      splinelabp->SetLineWidth(1);
      splinelabp->Write("",TObject::kOverwrite);
    }
    cout <<"creating splinelab target" << endl;
    splinelabt = boom->Ruthvslab(1, 179, steps,2); 
    if(splinelabt->GetXmax() < splinelabt->GetXmin() ){
      cout << " inverting spline xmax < xmin!" << endl;
      cout << splinelabt->GetNp() << endl;
      vector<double> x;
      vector<double> y;
      x.resize(splinelabt->GetNp());
      y.resize(splinelabt->GetNp());
      for(int i=0;i<splinelabt->GetNp();i++){
	splinelabt->GetKnot(splinelabt->GetNp()-i-1, x[i], y[i]);      
	txtfile << x[i] << "\t" << y[i] << endl;
      }
      TGraph* graph = new TGraph(splinelabt->GetNp(), &x[0], &y[0]);
      TSpline3* spline = new TSpline3("tar",graph);
      spline->SetLineColor(3);
      spline->SetName("tar");
      spline->SetLineWidth(1);
      spline->Write("",TObject::kOverwrite);
      delete spline;
    }
    else{
      cout << " not inverting spline xmax > xmin!" << endl;
      vector<double> x;
      vector<double> y;
      x.resize(splinelabt->GetNp());
      y.resize(splinelabt->GetNp());

      for(int i=0;i<splinelabt->GetNp();i++){
	splinelabt->GetKnot(i, x[i], y[i]);      
	txtfile << x[i] << "\t" << y[i] << endl;
      }
      splinelabt->SetLineColor(3);
      splinelabt->SetName("tar");
      splinelabt->SetLineWidth(1);
      splinelabt->Write("",TObject::kOverwrite);
    }
    if(rootFile != NULL){
      outfile->Close();
      cout <<"file closed" << endl;
    }
    TFile *fsharc = new TFile("sharc.root","RECREATE");
    splinecm->Write();
    fsharc->Close();
    delete boom;
    delete splinecm;
    delete splinelabp;
  }
  else if(cm2lab){
    cout << "cm2lab " << endl;
    //adjust here!!
    //Kinematics *boomboom[5];
    //boomboom[0] = new Kinematics(proj, targ, ebeam);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.886);
    //boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.100);
    //boomboom[4] = new Kinematics(proj, targ, reco, ejec, ebeam, 2.870);
    //Kinematics *boomboom[4];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.886);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.155);
    //boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 2.870);
    //Kinematics *boomboom[3];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
    //boomboom[1] = new Kinematics(proj, targ, ebeam);
    //boomboom[2] = new Kinematics(proj, targ, ebeam);
    //Kinematics *boomboom[4];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.570);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 2.710);
    //boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 3.892);
    
    //Te for KAKENHI
    Kinematics *boomboom[7];
    for(int i=0;i<7;i++)
      boomboom[i] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.0);

    //Ti for DOE proposal
    //Kinematics *boomboom[4];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.863);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.237);
    //boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.576);
    //Fe for DOE proposal
    //Kinematics *boomboom[2];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);

    //Xe for triumf 1508
    //Kinematics *boomboom[6];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[4] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[5] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);

    //Sr for triumf 1389
    //Kinematics *boomboom[3];
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 0.815);
    //boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.793);

    //for(int i=0;i<5;i++){
    //  cout << "boomboom["<<i<<"]->GetVcm(2) " << boomboom[i]->GetVcm(2) << " boomboom["<<i<<"]->GetBetacm(2) " << boomboom[i]->GetBetacm(2) << endl;
    //}
    cout << "ticm = " << boomboom[1]->GetCmEnergy(ebeam)-proj->GetMass()-targ->GetMass() << endl;
    TCanvas* div2 = new TCanvas("div2","div2", 800, 600);
    div2->Divide(2,2);
    TPad *pad = new TPad("pad","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    TPad *pad3 = new TPad("pad3","",0,0,1,1);
    TPad *pad4 = new TPad("pad4","",0,0,1,1);
   
    TAxis *xa, *ya;  
    char* Name = NULL;
    TFile *infile = new TFile(cmFile);    
    TIter next(infile->GetListOfKeys());
    TKey *key;
    TSpline3* splinecm;
    int j=0;
    while ((key=(TKey*)next()))
      if(strcmp(key->GetClassName(),"TSpline3") == 0)
	j++;
    TIter anext(infile->GetListOfKeys());
    vector<TGraph*> graph;
    vector<TGraph*> graphlab;
    vector<TGraph*> graphlab2;
    vector<TGraph*> graphlabcs;
    vector<TGraph*> graphlabcm;
    vector<TH1F*> histcm;
    vector<TH1F*> histlab;

    graph.resize(j);
    graphlab.resize(j);
    graphlab2.resize(j);
    graphlabcs.resize(j);
    graphlabcm.resize(j);
    histcm.resize(j);
    histlab.resize(j);
    
    j=0;
    bool inverse = false;
    while ((key=(TKey*)anext())){
      if(strcmp(key->GetClassName(),"TSpline3") == 0){
	Name = (char*)key->GetName();
	cout << Name << endl;
	infile->cd();
	splinecm = (TSpline3*)infile->Get(Name);
	outfile->cd();
	ofstream txtfile(Form("lab_%s.dat",Name));
	cout << splinecm->GetNp() << endl;
	vector<double> x;
	vector<double> y;
	vector<double> xlab;
	vector<double> ylab;
	vector<double> xlabr;
	vector<double> ylabr;
	vector<double> ylabcs;

	x.resize(splinecm->GetNp());
	y.resize(splinecm->GetNp());
	xlab.resize(splinecm->GetNp());
	ylab.resize(splinecm->GetNp());
	xlabr.resize(splinecm->GetNp());
	ylabr.resize(splinecm->GetNp());
	ylabcs.resize(splinecm->GetNp());

	double crosscm[3]={0.,0.,0.};
	double crosslab[3]={0.,0.,0.};
	inverse = false;
	histcm[j] = new TH1F(Form("hcm_%d",j),Form("hcm_%d",j),180,0,180);
	histlab[j] = new TH1F(Form("hlab_%d",j),Form("hlab_%d",j),180,0,180);
	if(splinecm->GetXmax() < splinecm->GetXmin() ){
	  crosslab[0]=0.;
	  crosslab[1]=0.;
	  crosslab[2]=0.;
	  crosscm[0]=0.;
	  crosscm[1]=0.;
	  crosscm[2]=0.;
	  cout << " inverting spline xmax < xmin!" << endl;
	  for(int i=0;i<splinecm->GetNp();i++){
	    splinecm->GetKnot(splinecm->GetNp()-i-1, x[i], y[i]);
	    if(i>0){
	      crosscm[0] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
	      if(x[i]>29&&x[i]<165)
		crosscm[1] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
	      if(x[i]>72&&x[i]<149){
		crosscm[2] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
		//cout << x[i] << "\t" << crosscm[2] << endl;
	      }
	    }
	    //cout << " theta " << x[i] << " cs " << y[i] << endl;  
	    xlab[i] = boomboom[j]->Angle_cm2lab(boomboom[j]->GetVcm(2),(180-x[i])*deg2rad)*rad2deg;
	    ylab[i] = boomboom[j]->Sigma_cm2lab(x[i]*deg2rad,y[i]);
	    if(i>0){
	      crosslab[0] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	      if(xlab[i]>trexmin&&xlab[i]<trexmax)
		crosslab[1] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	      if(xlab[i]>oldmin&&xlab[i]<oldmax)
		crosslab[2] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	    }
	    txtfile << xlab[i] << "\t" << ylab[i] << endl;
	    //cout << x[i] << "\t" << y[i] << "\t" << xlab[i] << "\t" << ylab[i] << endl;
	    histcm[j]->Fill(x[i],y[i]);
	    histlab[j]->Fill(xlab[i],ylab[i]);
	  }
	  graph[j] = new TGraph(splinecm->GetNp(), &x[0], &y[0]);
	  if(xlab[0] > xlab[splinecm->GetNp()]){
	    crosslab[0]=0.;
	    crosslab[1]=0.;
	    crosslab[2]=0.;
	    inverse = true;
	    cout << " inverting lab spline xmax < xmin!" << endl;
	    for(int i=0;i<splinecm->GetNp();i++){
	      xlabr[i]=xlab[splinecm->GetNp()-i-1];
	      ylabr[i]=ylab[splinecm->GetNp()-i-1];
	      if(i>0){
		crosslab[0] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
		if(xlabr[i]>trexmin&&xlabr[i]<trexmax)
		  crosslab[1] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
		if(xlabr[i]>oldmin&&xlabr[i]<oldmax)
		  crosslab[2] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
	      }
	    }
	  }
	  if(!inverse)
	    graphlab[j] = new TGraph(splinecm->GetNp(), &xlab[0], &ylab[0]);
	  else
	    graphlab[j] = new TGraph(splinecm->GetNp(), &xlabr[0], &ylabr[0]);
	  
	  graphlab2[j] = new TGraph(splinecm->GetNp(), &xlab[0], &ylab[0]);
	  graphlabcm[j] = new TGraph(splinecm->GetNp(), &x[0], &xlab[0]);
	  TSpline3* spline = new TSpline3(Form("cm_%d",j),graph[j]);
	  spline->SetLineColor(3);
	  spline->SetName(Form("cm_%d",j));
	  spline->SetLineWidth(1);
	  spline->Write("",TObject::kOverwrite);
	  delete spline;
	  TSpline3* splinel = new TSpline3(Form("lab_%d",j),graphlab[j]);
	  splinel->SetLineColor(2);
	  splinel->SetName(Form("lab_%d",j));
	  splinel->SetLineWidth(1);
	  splinel->Write("",TObject::kOverwrite);
	  delete spline;
	  graph[j]->Write(Form("dsdO_cm_%d",j),TObject::kOverwrite);
	  graphlab[j]->Write(Form("dsdO_lab_%d",j),TObject::kOverwrite);
	  graphlabcm[j]->Write(Form("lab_cm_%d",j),TObject::kOverwrite);
	}
	else{
	  cout << " not inverting spline xmax > xmin!" << endl;
	  crosslab[0]=0.;
	  crosslab[1]=0.;
	  crosslab[2]=0.;
	  crosscm[0]=0.;
	  crosscm[1]=0.;
	  crosscm[2]=0.;
	  for(int i=0;i<splinecm->GetNp();i++){
	    splinecm->GetKnot(i, x[i], y[i]);  
	    //if(j==1)
	    //  cout << " cm theta " << x[i] << " cs " << y[i] << endl;      
	    if(i>0){
	      crosscm[0] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
	      if(x[i]>29&&x[i]<165)
		crosscm[1] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
		if(x[i]>72&&x[i]<149){
		crosscm[2] += (y[i]*sin(x[i]*deg2rad)+y[i-1]*sin(x[i-1]*deg2rad))/2.*(x[i]-x[i-1])*deg2rad;
		//cout << x[i] << "\t" << crosscm[2] << endl;
	      }
	    }
	    xlab[i] = boomboom[j]->Angle_cm2lab(boomboom[j]->GetVcm(2),(180-x[i])*deg2rad)*rad2deg;
	    ylab[i] = boomboom[j]->Sigma_cm2lab(x[i]*deg2rad,y[i]);
	    if(i>0){
	      crosslab[0] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	      if(xlab[i]>trexmin&&xlab[i]<trexmax)
		crosslab[1] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	      if(xlab[i]>oldmin&&xlab[i]<oldmax)
		crosslab[2] += (ylab[i]*sin(xlab[i]*deg2rad)+ylab[i-1]*sin(xlab[i-1]*deg2rad))/2.*(xlab[i]-xlab[i-1])*deg2rad;
	    }
	    //cout << x[i] << "\t" << y[i] << "\t" << xlab[i] << "\t" << ylab[i] << endl;
	    txtfile << xlab[i] << "\t" << ylab[i] << endl;
	    histcm[j]->Fill(x[i],y[i]*sin(x[i]*deg2rad));
	    histlab[j]->Fill(xlab[i],ylab[i]*sin(xlab[i]*deg2rad));
	  }
	  graph[j] = new TGraph(splinecm->GetNp(), &x[0], &y[0]);
	  if(xlab[0] > xlab[splinecm->GetNp()]){
	    crosslab[0]=0.;
	    crosslab[1]=0.;
	    crosslab[2]=0.;
	    inverse = true;
	    cout << " inverting lab spline xmax < xmin!" << endl;
	    for(int i=0;i<splinecm->GetNp();i++){
	      xlabr[i]=xlab[splinecm->GetNp()-i-1];
	      ylabr[i]=ylab[splinecm->GetNp()-i-1];
	      if(i>0){
		crosslab[0] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
		if(xlabr[i]>trexmin&&xlabr[i]<trexmax)
		  crosslab[1] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
		if(xlabr[i]>oldmin&&xlabr[i]<oldmax)
		  crosslab[2] += (ylabr[i]*sin(xlabr[i]*deg2rad)+ylabr[i-1]*sin(xlabr[i-1]*deg2rad))/2.*(xlabr[i]-xlabr[i-1])*deg2rad;
	      }
	    }
	  }
	  if(!inverse)
	    graphlab[j] = new TGraph(splinecm->GetNp(), &xlab[0], &ylab[0]);
	  else
	    graphlab[j] = new TGraph(splinecm->GetNp(), &xlabr[0], &ylabr[0]);
	  
	  graphlab2[j] = new TGraph(splinecm->GetNp(), &xlab[0], &ylab[0]);
	  graphlabcm[j] = new TGraph(splinecm->GetNp(), &x[0], &xlab[0]);
	  TSpline3* spline = new TSpline3(Form("cm_%d",j),graph[j]);
	  spline->SetLineColor(3);
	  spline->SetName(Form("cm_%d",j));
	  spline->SetLineWidth(1);
	  spline->Write("",TObject::kOverwrite);
	  delete spline;
	  TSpline3* splinel = new TSpline3(Form("lab_%d",j),graphlab[j]);
	  splinel->SetLineColor(2);
	  splinel->SetName(Form("lab_%d",j));
	  splinel->SetLineWidth(1);
	  splinel->Write("",TObject::kOverwrite);
	  delete splinel;
	  graph[j]->Write(Form("dsdO_cm_%d",j),TObject::kOverwrite);
	  graphlab[j]->Write(Form("dsdO_lab_%d",j),TObject::kOverwrite);
	  graphlabcm[j]->Write(Form("lab_cm_%d",j),TObject::kOverwrite);
	}
	//cout << j << " cm " << histcm[j]->Integral() << endl;
	//cout << j << " lab " << histlab[j]->Integral() << endl;
	cout << j << " lab all " << crosslab[0]*2*(TMath::Pi()) << "\tcm all " << crosscm[0]*2*(TMath::Pi()) << endl;
	cout << j << " lab trex " << crosslab[1]*2*(TMath::Pi()) << "\tcm trex " << crosscm[1]*2*(TMath::Pi()) << endl;
	cout << j << " lab old " << crosslab[2]*2*(TMath::Pi()) << "\tcm old " << crosscm[2]*2*(TMath::Pi()) << endl;

	div2->cd(1);
	if(j==0){
	  pad->Draw();
	  pad->cd();
	  //pad->SetLogy();
	  xa = graph[0]->GetXaxis();
	  xa->SetLimits(0, 180);
	  xa->SetTitle("#theta_{cm} [#circ]");
	  xa->SetTitleOffset(0.9);
	  xa->SetTitleSize(0.05);
	  xa->SetLabelOffset(0.008);
	  xa->SetLabelSize(0.04);
	  xa->SetTickLength(0.03);
	  
	  ya = graph[0]->GetYaxis();
	  //ya->SetRangeUser(1e-4,1e1);
	  ya->SetTitle("cross section [mb/sr]");
	  ya->SetTitleOffset(0.9);
	  ya->SetTitleSize(0.05);
	  ya->SetLabelOffset(0.005);
	  ya->SetLabelSize(0.04);
	  ya->SetTickLength(0.03);

	  graph[0]->Draw("AL");
	}
	else{
	  pad->cd();
	  graph[j]->SetLineColor(1+j);
	  graph[j]->Draw("L");
	}
	div2->cd(2);
	if(j==0){
	  pad2->Draw();
	  pad2->cd();
	  //pad2->SetLogy();
	  xa = graphlab[0]->GetXaxis();
	  xa->SetLimits(0, 180);
	  xa->SetTitle("#theta_{lab} [#circ]");
	  xa->SetTitleOffset(0.9);
	  xa->SetTitleSize(0.05);
	  xa->SetLabelOffset(0.008);
	  xa->SetLabelSize(0.04);
	  xa->SetTickLength(0.03);
	  
	  ya = graphlab[0]->GetYaxis();
	  ya->SetRangeUser(1e-4,1e1);
	  ya->SetTitle("cross section [mb/sr]");
	  ya->SetTitleOffset(0.9);
	  ya->SetTitleSize(0.05);
	  ya->SetLabelOffset(0.005);
	  ya->SetLabelSize(0.04);
	  ya->SetTickLength(0.03);
	  
	  graphlab[0]->Draw("AL");
	}
	else{
	  pad2->cd();
	  graphlab[j]->SetLineColor(1+j);
	  graphlab[j]->Draw("L");
	  
	}
	div2->cd(3);
	if(j==0){
	  pad3->Draw();
	  pad3->cd();
	  xa = graphlab2[0]->GetXaxis();
	  xa->SetLimits(0, 180);
	  xa->SetTitle("#theta_{lab} [#circ]");
	  xa->SetTitleOffset(0.9);
	  xa->SetTitleSize(0.05);
	  xa->SetLabelOffset(0.008);
	  xa->SetLabelSize(0.04);
	  xa->SetTickLength(0.03);
	  
	  ya = graphlab2[0]->GetYaxis();
	  ya->SetRangeUser(0,0.4);
	  ya->SetTitle("cross section [mb/sr]");
	  ya->SetTitleOffset(0.9);
	  ya->SetTitleSize(0.05);
	  ya->SetLabelOffset(0.005);
	  ya->SetLabelSize(0.04);
	  ya->SetTickLength(0.03);
	  
	  graphlab2[0]->Draw("AL");
	}
	else{
	  pad3->cd();
	  graphlab2[j]->SetLineColor(1+j);
	  graphlab2[j]->Draw("L");
	}
	div2->cd(4);
	if(j==0){
	  pad4->Draw();
	  pad4->cd();
	  xa = graphlabcm[0]->GetXaxis();
	  xa->SetLimits(0, 180);
	  xa->SetTitle("#theta_{cm} [#circ]");
	  xa->SetTitleOffset(0.9);
	  xa->SetTitleSize(0.05);
	  xa->SetLabelOffset(0.008);
	  xa->SetLabelSize(0.04);
	  xa->SetTickLength(0.03);
	  
	  ya = graphlabcm[0]->GetYaxis();
	  ya->SetRangeUser(0,180);
	  ya->SetTitle("#theta_{lab} [#circ]");
	  ya->SetTitleOffset(0.9);
	  ya->SetTitleSize(0.05);
	  ya->SetLabelOffset(0.005);
	  ya->SetLabelSize(0.04);
	  ya->SetTickLength(0.03);
	  
	  graphlabcm[0]->Draw("AL");
	}
	else{
	  pad4->cd();
	  graphlabcm[j]->SetLineColor(1+j);
	  graphlabcm[j]->Draw("L");
	}

	j++;

      }
    }
    infile->Close();
    outfile->Close();
    div2->Print(Form("%s",psFile));
  }
  else if(iso){
    cout << "isotropic distribution in lab system transformed to cm for simulation " << endl;
    outfile->cd();
    vector<Kinematics*> boomboom;   
    boomboom.resize(numb);
    //adjust here
    //boomboom[0] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.017);
    //boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 3.221);
    //boomboom[0] = new Kinematics(proj, targ, targ, proj, ebeam, 0.000);
    //boomboom[1] = new Kinematics(proj, targ, targ, proj, ebeam, 0.000);
    if(coulex.size() > 0){
      boomboom[0] = new Kinematics(proj, targ, ebeam);
      boomboom[1] = new Kinematics(proj, targ, ebeam);
    }
    else{
      boomboom[1] = new Kinematics(proj, targ, reco, ejec, ebeam, 1.570);
      boomboom[2] = new Kinematics(proj, targ, reco, ejec, ebeam, 2.710);
      boomboom[3] = new Kinematics(proj, targ, reco, ejec, ebeam, 3.892);
    }
    vector<TGraph*> graph;
    vector<TGraph*> graphlab;
    vector<TSpline3*> spline;
    vector<TSpline3*> splinelab;
    vector<int> ctr;
    graph.resize(numb);
    graphlab.resize(numb);
    spline.resize(numb);
    splinelab.resize(numb);
    ctr.resize(numb);

    vector<vector<double> > x;
    vector<vector<double> > y;
    vector<vector<double> > xr;
    vector<vector<double> > yr;
    vector<vector<double> > xlab;
    vector<vector<double> > ylab;

    x.resize(numb);
    y.resize(numb);
    xr.resize(numb);
    yr.resize(numb);
    xlab.resize(numb);
    ylab.resize(numb);
    for(UShort_t s=0;s<numb;s++){
      x[s].resize(180);
      y[s].resize(180);
      xr[s].resize(180);
      yr[s].resize(180);
      xlab[s].resize(180);
      ylab[s].resize(180);
    }


    bool inverse = false;
    for(int j=0;j<numb;j++){
      for(int i=0;i<179;i++){
	xlab[j][i] = i+1;
	ylab[j][i] = 1;
	x[j][i] = boomboom[j]->Angle_lab2cm(boomboom[j]->GetVcm(2),xlab[j][i]*deg2rad)*rad2deg;
	y[j][i] = boomboom[j]->Sigma_lab2cm(x[j][i]*deg2rad,ylab[j][i]); //was xlab!!
	if(coulex.size() > 0){
	  x[0][i] = boomboom[j]->Angle_lab2cm(boomboom[j]->GetVcm(3),xlab[j][i]*deg2rad)*rad2deg;
	  x[1][i] = boomboom[j]->Angle_lab2cm(boomboom[j]->GetVcm(2),xlab[j][i]*deg2rad)*rad2deg;
	  y[0][i] = boomboom[j]->Sigma_lab2cm((180.-x[0][i])*deg2rad,ylab[j][i]);
	  y[1][i] = boomboom[j]->Sigma_lab2cm(x[1][i]*deg2rad,ylab[j][i]);

	}
	if(x[j][i]>0 && x[j][i]<180 && y[j][i]>0){
	  ctr[j]++;
	//cout << "xlab[j][i] " << xlab[j][i] << " ylab[j][i] " << ylab[j][i]; 
	  //cout << " x["<<j<<"]["<<i<<"] " << x[j][i] << "  y["<<j<<"]["<<i<<"] " << y[j][i] <<endl; 
	}
	//cout << " re xlab" << boomboom[j]->Sigma_cm2lab(x[j][i]*deg2rad,y[j][i]) << " re ylab " <<  boomboom[j]->Angle_cm2lab(boomboom[j]->GetVcm(2),x[j][i]*deg2rad)*rad2deg << endl;
      }
      //cout << "x["<<j<<"][0] " << x[j][0] << " x["<<j<<"]["<<ctr[j-1]<<"] " << x[j][ctr[j]-1] << endl ;
      if(x[j][0]>x[j][ctr[j]-1]){
	inverse =true;
	cout << " inverting lab spline xmax < xmin! " << j  << endl;
	for(int i=0;i<ctr[j];i++){
	  xr[j][ctr[j]-i-1]=x[j][i];
	  yr[j][ctr[j]-i-1]=y[j][i];
	  
	}
	
      }
      graph[j] = new TGraph(180, &x[j][0], &y[j][0]);
      graphlab[j] = new TGraph(180, &xlab[j][0], &ylab[j][0]);
      if(!inverse)
	spline[j] = new TSpline3(Form("cm_%d",j),&x[j][0],&y[j][0],ctr[j]);
      else
	spline[j] = new TSpline3(Form("cm_%d",j),&xr[j][0],&yr[j][0],ctr[j]);
	
      spline[j]->SetLineColor(3);
      spline[j]->SetName(Form("cm_%d",j));
      spline[j]->SetLineWidth(1);
      spline[j]->Write("",TObject::kOverwrite);
      splinelab[j] = new TSpline3(Form("lab_%d",j),&xlab[j][0],&ylab[j][0],180);
      splinelab[j]->SetLineColor(2);
      splinelab[j]->SetName(Form("lab_%d",j));
      splinelab[j]->SetLineWidth(1);
      splinelab[j]->Write("",TObject::kOverwrite);
      graph[j]->Write(Form("dsdO_cm_%d",j),TObject::kOverwrite);
      graphlab[j]->Write(Form("dsdO_lab_%d",j),TObject::kOverwrite);
      
    }
    for(int i=0;i<180;i++){
      for(int j=0;j<numb;j++){
	cout << x[j][i] << "\t" << y[j][i] << "\t";
      }
      cout << endl;
    }
    outfile->Close();
  }
  else if(test){
    Kinematics *coulex = new Kinematics(proj, targ, ebeam);
    double ejectileangle[90];
    double ejectileenergy[90];
    double recoilangle[90];
    double recoilenergy[90];
    double calrecoilangle[90];
    double calrecoilenergy[90];
    double maxen = coulex->ELab(coulex->GetMaxAngle(3)-0.000001,3);
    cout << "energy at max angle " << maxen << endl;
    for(int i=0;i<90;i++){
      coulex->Final((i+0.1)*TMath::Pi()/180.,2);
      ejectileangle[i]=coulex->GetThetalab(3)*180./TMath::Pi();
      recoilangle[i]=i+0.1;
      ejectileenergy[i]=coulex->GetTlab(3);//(coulex->ELab(coulex->GetThetalab(3),3);
      recoilenergy[i]=coulex->GetTlab(2);//coulex->ELab(i*TMath::Pi()/180.,2);
      //cout << " recoil " << recoilangle[i] << " energy " << recoilenergy[i] << " ejectile " << ejectileangle[i] << " energy " << ejectileenergy[i] << endl;
      coulex->SetAngles(ejectileangle[i]*TMath::Pi()/180., 3,ejectileenergy[i]>maxen);
      calrecoilangle[i]=coulex->GetThetalab(2)*180./TMath::Pi();
      calrecoilenergy[i]=coulex->GetTlab(2);
      //cout << " recoil " << calrecoilangle[i] << " energy " << calrecoilenergy[i] << endl;
    }
    TCanvas* div2 = new TCanvas("div2","div2", 800, 600);
    div2->cd();
    TGraph *recoilgr = new TGraph(90,recoilangle,recoilenergy);
    recoilgr->SetLineColor(3);
    recoilgr->Draw("AL");
    TGraph *ejectilegr = new TGraph(90,ejectileangle,ejectileenergy);
    ejectilegr->SetLineColor(2);
    ejectilegr->Draw("L");
    TGraph *calrecoilgr = new TGraph(90,calrecoilangle,calrecoilenergy);
    calrecoilgr->SetLineColor(4);
    calrecoilgr->Draw("L");
    cout << "max angle " << coulex->GetMaxAngle(3)*180./TMath::Pi() << endl;
    div2->Print(Form("%s",psFile));
    /*
    cout << "---------------------------------" << endl;
    for(int i=0;i<90;i++){
      cout << " ejectile " << ejectileangle[i] << " energy " << coulex->ELab(ejectileangle[i],2);
      cout << " recoil " << coulex->GetThetalab(2) << " energy " << coulex->ELab(coulex->GetThetalab(2),2)<<endl;
    }
    */
    /*
    Kinematics *boomboom = new Kinematics(proj, targ, reco, ejec, ebeam, 0.000);
    double angle;
    double sigma;
    


    for(int i=0; i<180; i++){
      angle = i*deg2rad;
      sigma =1;
      boomboom->Transform2cm(angle, sigma);
      cout << angle*rad2deg << "\t" << sigma << endl;
    }
    */
  }
  if(rootFile != NULL){
     if(outfile != NULL)
       outfile->Close();
  }
  cout << "end of prog" << endl;
}

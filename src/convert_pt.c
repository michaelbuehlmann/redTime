#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_spline.h>

double get_hubble(int model_no, char* params_file);
double get_f_cb(int model_no, char* red_file);
int process_PT_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *D, double *Pk_nu, char* red_file);

int main(int argc, char**argv) {
  if(argc != 4) {
    fprintf(stderr, "Invalid arguments\n");
    return -1;
  }
  double kmin = 1e-3;
  double kmax = 5.;
  //  int nk = 351;
  int nk = 3000;
  double *k = (double *)malloc(nk*sizeof(double));

  int model_no = 1;
  int step_no = atoi(argv[1]);

  char* params_file = argv[2]; // direct path to params file
  char* red_file = argv[3]; // path to redTime output folder

  int i;
  int nk1 = 50;
  int nk2 = 200;

  // this sets up the k spacing used for the MT emulator
  for (i = 0; i < nk1; i++) {
    k[i] = log10(kmin)+(float)i*(log10(0.04)-log10(kmin))/nk1;
    k[i] = pow(10.,k[i]);
  }

  for (i = nk1; i < nk2; i++) {
    k[i] = 0.04 + (float)(i-nk1)*(0.200 - 0.04)/(nk2-nk1-1);
  }

  for(i = nk2; i < nk; i++) {
    k[i] = log10(0.201)+(float)(i-nk2)*(log10(kmax)-log10(0.201))/((nk-nk2)-1);
    k[i] = pow(10.,k[i]);
  }

  for (model_no = 1; model_no < 30; model_no++) {
    double *Pk_pt = (double *)malloc(nk*sizeof(double));
    double *D = (double *)malloc(nk*sizeof(double));
    double *Pk_nu = (double *)malloc(nk*sizeof(double));

    double *k_pt = (double *)malloc(nk*sizeof(double));

    double h = get_hubble(model_no, params_file);   // get h for each model from the design
    double f_cb = get_f_cb(model_no, params_file);  // calculate f_cb for each model


    int nk_pt = process_PT_runs(h, model_no, step_no, nk, k_pt, Pk_pt, D, Pk_nu, red_file); //read D+, pk_pt, pk_nu


    int kk;
    char k_filename[256], pk_filename[256];

    // these will be the output files, come up with a naming convention or input one
    // note to self, need to create sub-directories for redshifts
    sprintf(k_filename, "%s/STEP%d/k_M%03d_no_interp_test.dat",red_file, step_no, model_no);
    sprintf(pk_filename, "%s/STEP%d/pk_M%03d_no_interp_test.dat",red_file, step_no, model_no);

    FILE *fp_k = fopen(k_filename, "w");
    FILE *fp_pk = fopen(pk_filename, "w");

    gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,nk_pt);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp, k_pt, D, nk_pt);


    for (kk = 0; kk < nk; kk++) {
      //Correct PT by a factor of f_cb due to differing definitions of delta in Amol's PT
      //and the power spectrum output from HACC
      fprintf(fp_pk, "%lf ", Pk_pt[kk]*f_cb*f_cb);
      fprintf(fp_k, "%lf ", k_pt[kk]);
    }
    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);

    free(Pk_pt);  free(D); free(Pk_nu);
    free(k_pt);
    fclose(fp_k); fclose(fp_pk);
  }
  free(k);
  return(0);
}

double get_hubble(int model_no, char* params_file) {
  double h;
  double Om, Omb, Omnu, w0, wa, s8, n_s;
  char name[5];
  FILE *fp = fopen(params_file,"r"); // this is just parameters file

  if (fp == NULL) {
    fprintf(stderr, "Couldn't read: %s\n", params_file);
    return(-1);
  }
  //      M000   0.1335   0.02258    0.8   0.71   0.963   -1.0   0.0   0.0
  char tmp[300];
  for (int i=0; i < 5; i++){
          fgets(tmp,300,fp);
  }

  for (int i=0; i < model_no; i++) {
    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf", name, &Om, &Omb, &s8, &h, &n_s,  &w0, &wa, &Omnu);
  }

  fclose(fp);
  return(h);
}

double get_f_cb(int model_no, char* params_file) {
  double f_cb;
  double h, Om, Omb, Omnu, w0, wa, s8, n_s;

  char name[5];
  FILE *fp = fopen(params_file,"r");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't read: %s\n", params_file);
    return(-1);
  }
  //      M000   0.1335   0.02258    0.8   0.71   0.963   -1.0   0.0   0.0
  char tmp[300];
  for (int i=0; i < 5; i++){
	  fgets(tmp,300,fp);
  }
  for (int i=0; i < model_no; i++) {
    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf", name, &Om, &Omb, &s8, &h, &n_s,  &w0, &wa, &Omnu);
  }

  f_cb = (Om-Omnu)/Om;
  fclose(fp);
  return(f_cb);
}


int process_PT_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *D, double *Pk_nu, char* red_file) {
  char sed_call[256];

  sprintf(sed_call, "sed '/^#/ d' %s/redTime_M%03d.dat > junk.dat", red_file, model_no);

  system(sed_call);
  FILE *fp = fopen("junk.dat", "r");

  // need to confirm this
  int nread = 1, nk_pt=128, nz=33;//nk_pt=202, nz=26;

  double any[15];
  size_t nchar = 1000;
  char *line  = (char *)malloc((nchar+1)*sizeof(char));

  double *k_h = (double *) malloc(nk_pt*sizeof(double));
  double *Pk_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *Pk_nu_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *D_h = (double *) malloc(nz*nk_pt*sizeof(double));

  int z_no;
  int steps[8] = {163, 189, 247, 300, 347, 401, 453, 499};
  int output_z[8] = {9,11,14,18, 24, 28, 31, 32};

  for (int i = 0; i < 8; i++)
    if (step_no == steps[i])
      z_no = i;


  int i = 0, j = 0;
  while (nread!=EOF) {
    // save values
    nread = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &(k_h[j]), &(D_h[nk_pt*i+j]), &(any[1]), &(any[2]), &(any[3]), &(any[4]), &(Pk_nu_h[nk_pt*i+j]), &(Pk_h[nk_pt*i+j]), &(any[6]), &(any[7]), &(any[8]), &(any[9]), &(any[10]), &(any[11]), &(any[12]), &(any[13]), &(any[14]));
    if (nread!=EOF) {
      k_h[j] *= h;
      Pk_h[nk_pt*i+j] /= pow(h,3.);
      Pk_nu_h[nk_pt*i+j] /=pow(h,3.);
    }

    if (j == nk_pt-1) {
      // next input will be new redshift, start over
      i++;
      j = 0;

    } else {
      j++;
    }
  }
  fclose(fp);

  double D0 = D_h[nk_pt*(output_z[z_no]+1)-1];

  for (int kk=0; kk < nk; kk++) {
    if (kk < nk_pt) {
      k[kk] = k_h[kk];
      Pk[kk] = Pk_h[nk_pt*output_z[z_no]+kk];
      D[kk] = D_h[nk_pt*output_z[z_no]+kk]/D0;
      Pk_nu[kk] = Pk_nu_h[nk_pt*output_z[z_no]+kk];
    } else {
      k[kk] = 0;
      Pk[kk] = 0;
      D[kk] = 1.;
      Pk_nu[kk] = 0;
    }
  }

  free(k_h); free(Pk_h); free(D_h); free(Pk_nu_h);
  free(line);
  return(nk_pt);
}

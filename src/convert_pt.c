#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_spline.h>

double get_hubble(int model_no, char* params_file);
double get_f_cb(int model_no, char* red_file);
void process_PT_runs(double h, int model_no, int step_no, double *k_pt, double *Pk, double *D, double *Pk_nu, char* red_file, int nk_pt);

int main(int argc, char**argv) {
  if(argc != 6) {
    fprintf(stderr, "Invalid arguments\n");
    return -1;
  }
  int model_no = atoi(argv[1]); // total number of models to evaluate
  int step_no = atoi(argv[2]); // analysis step number to run
  int nk_pt = atoi(argv[3]); // number of k points per perturbation theory output, default to 128
  char* params_file = argv[4]; // direct path to params file
  char* red_file = argv[5]; // path to redTime output folder

  // make output directory
  char mkdir_call[256];
  sprintf(mkdir_call, "mkdir -p %s/STEP%d", red_file, step_no);
  {
    int res = system(mkdir_call);
    if(res)
      abort();
  }

  for (int mn = 1; mn < (model_no+1); mn++) {
    double *Pk_pt = (double *)malloc(nk_pt*sizeof(double));
    double *D = (double *)malloc(nk_pt*sizeof(double));
    double *Pk_nu = (double *)malloc(nk_pt*sizeof(double));
    double *k_pt = (double *)malloc(nk_pt*sizeof(double));
    double h = get_hubble(mn, params_file);   // get h for each model from the design
    double f_cb = get_f_cb(mn, params_file);  // calculate f_cb for each model
    process_PT_runs(h, mn, step_no, k_pt, Pk_pt, D, Pk_nu, red_file, nk_pt); //read D+, pk_pt, pk_nu
    int kk;

    // create output files
    char k_filename[256], pk_filename[256];
    sprintf(k_filename, "%s/STEP%d/k_M%03d_no_interp_test.dat",red_file, step_no, mn);
    sprintf(pk_filename, "%s/STEP%d/pk_M%03d_no_interp_test.dat",red_file, step_no, mn);
    FILE *fp_k = fopen(k_filename, "w");
    FILE *fp_pk = fopen(pk_filename, "w");
    // interpolate PT power spectra
    gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,nk_pt);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    gsl_spline_init(sp, k_pt, D, nk_pt);

    for (kk = 0; kk < nk_pt; kk++) {
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

  return(0);
}

double get_hubble(int model_no, char* params_file) {
  double h;
  double Om, Omb, Omnu, w0, wa, s8, n_s;
  char name[5];
  FILE *fp = fopen(params_file,"r");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't read: %s\n", params_file);
    return(-1);
  }
  // skip comment lines
  char tmp[300];
  for (int i=0; i < 5; i++){
    char* t = fgets(tmp,300,fp);
    if(t != tmp)
      abort();
  }
  // read cosmological parameters
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
  // skip comment lines
  char tmp[300];
  for (int i=0; i < 5; i++){
      char* t = fgets(tmp,300,fp);
      if(t != tmp)
        abort();
  }
  // read cosmological parameters
  for (int i=0; i < model_no; i++) {
    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf", name, &Om, &Omb, &s8, &h, &n_s,  &w0, &wa, &Omnu);
  }

  f_cb = (Om-Omnu)/Om;
  fclose(fp);
  return(f_cb);
}

void process_PT_runs(double h, int model_no, int step_no, double *k_pt, double *Pk, double *D, double *Pk_nu, char* red_file, int nk_pt) {
  char sed_call[256];
  sprintf(sed_call, "sed '/^#/ d' %s/redTime_M%03d.dat > junk.dat", red_file, model_no);
  {
    int res = system(sed_call);
    if(res)
      abort();
  }

  FILE *fp = fopen("junk.dat", "r");
  int nread = 1, nz=33;
  double any[15];
  size_t nchar = 1000;
  char *line  = (char *)malloc((nchar+1)*sizeof(char));

  double *k_h = (double *) malloc(nk_pt*sizeof(double));
  double *Pk_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *Pk_nu_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *D_h = (double *) malloc(nz*nk_pt*sizeof(double));


  int steps[8] = {163, 189, 247, 300, 347, 401, 453, 499};
  int output_z[8] = {9,11,14,18, 24, 28, 31, 32};

  int z_no = 0;
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

  for (int kk=0; kk < nk_pt; kk++) {
    if (kk < nk_pt) {
      k_pt[kk] = k_h[kk];
      Pk[kk] = Pk_h[nk_pt*output_z[z_no]+kk];
      D[kk] = D_h[nk_pt*output_z[z_no]+kk]/D0;
      Pk_nu[kk] = Pk_nu_h[nk_pt*output_z[z_no]+kk];
    }
  }
  free(k_h); free(Pk_h); free(D_h); free(Pk_nu_h);
  free(line);
  return;
}

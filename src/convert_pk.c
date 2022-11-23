#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<gsl/gsl_spline.h>

double get_hubble(int model_no);
double get_f_cb(int model_no);
int process_PT_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *D, double *Pk_nu);
int process_PM_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *Pk_err);
int process_HACC_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *Pk_err);

int main(int argc, char**argv) {
  if(argc != 2) {
    fprintf(stderr, "Invalid arguments\n");
    return -1;
  }
  double kmin = 1e-3;
  double kmax = 5.;
  //  int nk = 351;
  int nk = 3000;
  double *k = (double *)malloc(nk*sizeof(double));

  int model_no = 0;
  int step_no = atoi(argv[1]);

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

  for (model_no = 0; model_no < 112; model_no++) {
    double *Pk_pt = (double *)malloc(nk*sizeof(double));
    double *Pk_pm = (double *)malloc(16*nk*sizeof(double));
    double *Pk_pm_err = (double *)malloc(16*nk*sizeof(double));
    double *Pk_hacc = (double *)malloc(nk*sizeof(double));
    double *Pk_hacc_err = (double *)malloc(nk*sizeof(double));
    double *D = (double *)malloc(nk*sizeof(double));
    double *Pk_nu = (double *)malloc(nk*sizeof(double));

    double *k_pt = (double *)malloc(nk*sizeof(double));
    double *k_pm = (double *)malloc(16*nk*sizeof(double));
    double *k_hacc = (double *)malloc(nk*sizeof(double));

    double h = get_hubble(model_no);   // get h for each model from the design
    double f_cb = get_f_cb(model_no);  // calculate f_cb for each model

    int nk_pt = process_PT_runs(h, model_no, step_no, nk, k_pt, Pk_pt, D, Pk_nu); //read D+, pk_pt, pk_nu
    process_PM_runs(h, model_no, step_no, nk, k_pm, Pk_pm, Pk_pm_err); //read pk_PM
    process_HACC_runs(h, model_no, step_no, nk, k_hacc, Pk_hacc, Pk_hacc_err);// read pk_hacc

    int kk, pm_no;
    char k_filename[256], pk_filename[256], err_filename[256];

    // sprintf(filename, "/Users/jkwan/emulators/PT/linear_nu/STEP%d/pk_linear_nu_M%03d_test.dat", step_no, model_no);

    sprintf(k_filename, "/Users/jkwan/emulators/STEP%d/k_M%03d_no_interp_test.dat", step_no,model_no);
    sprintf(pk_filename, "/Users/jkwan/emulators/STEP%d/pk_M%03d_no_interp_test.dat", step_no,model_no);
    sprintf(err_filename, "/Users/jkwan/emulators/STEP%d/err_M%03d_no_interp_test.dat", step_no,model_no);


    FILE *fp_k = fopen(k_filename, "w");
    FILE *fp_pk = fopen(pk_filename, "w");
    FILE *fp_err = fopen(err_filename, "w");

    gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,nk_pt);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();


    gsl_spline_init(sp, k_pt, D, nk_pt);


    for (kk = 0; kk < nk; kk++) {
      //Correct PT by a factor of f_cb due to differing definitions of delta in Amol's PT
      //and the power spectrum output from HACC
      fprintf(fp_pk, "%lf ", Pk_pt[kk]*f_cb*f_cb);
      fprintf(fp_k, "%lf ", k_pt[kk]);

      // This part applies the growth factor correction to P_cb from the PM runs

      for (pm_no = 0; pm_no < 16; pm_no++) {
        double D_interp;
        if ((k_pm[kk] < k_pt[nk_pt-1]) && (k_pt[kk]!=0))
          D_interp = gsl_spline_eval(sp, k_pm[kk], acc);
        else
          D_interp = 1;

        fprintf(fp_k, "%lf ", k_pm[kk]);
        fprintf(fp_pk, "%lf ", Pk_pm[16*kk+pm_no]*D_interp*D_interp);
        fprintf(fp_err, "%lf ", Pk_pm_err[16*kk+pm_no]*D_interp*D_interp);
      }

      // This part applies the growth factor correction to P_cb from the HACC runs
      double D_interp;
      if ((k_hacc[kk] < k_pt[nk_pt-1])&&(k_pt[kk]!=0))
        D_interp = gsl_spline_eval(sp, k_hacc[kk], acc);
      else
        D_interp = 1;

      fprintf(fp_k, "%lf ", k_hacc[kk]);
      fprintf(fp_pk, "%lf ", Pk_hacc[kk]*D_interp*D_interp);
      fprintf(fp_err, "%lf ", Pk_hacc_err[kk]*D_interp*D_interp);
      fprintf(fp_k, "\n");  fprintf(fp_pk, "\n");
      fprintf(fp_err, "\n");
    }
    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);

    free(Pk_pt); free(Pk_pm); free(Pk_hacc); free(D); free(Pk_nu);
    free(k_pt); free(k_pm); free(k_hacc);
    fclose(fp_k); fclose(fp_pk); fclose(fp_err);
  }
  free(k);
  return(0);
}

double get_hubble(int model_no) {
  double h;
  double Om, Omb, Omnu, w0, wa, s8, n_s;
  char name[5];
  char filename[256];
  FILE *fp = fopen("/Users/jkwan/emulators/design.dat","r");

  if (fp == NULL) {
    fprintf(stderr, "Couldn't read: %s\n", filename);
    return(-1);
  }
  //      M000   0.1335   0.02258    0.8   0.71   0.963   -1.0   0.0   0.0
  for (int i=0; i < model_no+1; i++) {
    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf", name, &Om, &Omb, &s8, &h, &n_s,  &w0, &wa, &Omnu);
  }

  fclose(fp);
  return(h);
}

double get_f_cb(int model_no) {
  double f_cb;
  double h, Om, Omb, Omnu, w0, wa, s8, n_s;

  char filename[256], name[5];
  FILE *fp = fopen("/Users/jkwan/emulators/design.dat","r");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't read: %s\n", filename);
    return(-1);
  }
  //      M000   0.1335   0.02258    0.8   0.71   0.963   -1.0   0.0   0.0
  for (int i=0; i < model_no+1; i++) {
    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf", name, &Om, &Omb, &s8, &h, &n_s,  &w0, &wa, &Omnu);
  }

  f_cb = (Om-Omnu)/Om;
  fclose(fp);
  return(f_cb);
}


int process_HACC_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *Pk_err) {
  int ncol;
  char filename_base[256], filename[256];
  char * buf = (char *)malloc(256*sizeof(char));

  sprintf(filename, "/Users/jkwan/emulators/MiraU/matter_power/M%03d/L2100/HACC000/analysis/Pow/m%03d.pk.%3d.6400.pk", model_no, model_no, step_no);
  FILE *fp_test = fopen(filename,"r");
  if (fp_test == NULL) // TitanU
      sprintf(filename_base, "/Users/jkwan/emulators/TitanU/matter_power");
  else  // MiraU
      sprintf(filename_base, "/Users/jkwan/emulators/MiraU/matter_power");

  sprintf(filename, "%s/M%03d/L2100/HACC000/analysis/Pow/m%03d.pk.%03d.6400.pk", filename_base, model_no, model_no, step_no);

  if (model_no == 0)
    sprintf(filename, "/Users/jkwan/emulators/M000/m%03d.pk.%3d.6400.pk", model_no, step_no);

  fprintf(stderr, "%s\n", filename);

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Can't open file: %s\n", filename);
    return(-1);
  } else {
    fgets (buf, 256, fp);
  }

  char delims[] = "[";
  char *header = NULL;

  // Count the columns in the file.
  // strsep looks for columns within buf separated by the [ character
  // This is because sometimes files have an extra column of zeros on the end.

  ncol = 0;
  // sometimes there is no header, only count columns if there is a header
  if (buf[0]!='#') {
    ncol = 4;
    printf("%c\n", buf[0]);
  } else {
    while ((header = strsep(&buf, delims)) != NULL) {
      ncol++;
    }
  }

  free(buf);

  int nread = 1;
  int n = 0;
  double *any = (double *) malloc(ncol*sizeof(double));
  double *k_h  = (double *) malloc(8000*sizeof(double));
  double *Pk_h = (double *) malloc(8000*sizeof(double));
  double *err_h = (double *) malloc(8000*sizeof(double));

  fprintf(stderr, "ncol = %d\n", ncol);
  if (model_no == 0)
    ncol = 3;

  while (nread == 1) {
    for (int i = 0; i < ncol; i++)
      nread = fscanf(fp, "%lf", &any[i]);

    //      nread = fscanf(fp, "%lf %lf %lf %lf", &(k_h[n]), &(Pk_h[n]), &any1, &any2);
    k_h[n] = any[0]*h;
    Pk_h[n] = any[1]/pow(h,3.);

    //      err_h[n] = Pk_h[n]/sqrt(k_counts[n]);
    err_h[n] = Pk_h[n]/sqrt(any[2]);

    // This is any[4] for M048, M059,  M064, M080, M091-M111
    // any[2] for M061-M062, M064
    // and any[3] for everything else
    // Note also that some of the columns change between snapshots.
    n++;
    }
  n--;

  gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  gsl_spline *sp_err = gsl_spline_alloc(gsl_interp_cspline,n);
  gsl_interp_accel *acc_err = gsl_interp_accel_alloc();

  gsl_spline_init(sp, k_h, Pk_h, n);
  gsl_spline_init(sp_err, k_h, err_h, n);

  // for (kk = 0; kk < nk; kk++)
  //   {
  //     if (k[kk] > k_h[0] && k[kk] < k_h[n-1])
  // 	{
  // 	  Pk[kk] = gsl_spline_eval(sp, k[kk], acc);
  // 	  Pk_err[kk] = gsl_spline_eval(sp_err, k[kk], acc_err);
  // 	}
  //     else
  // 	{
  // 	  Pk[kk] = 0.;
  // 	  Pk_err[kk] = 0;
  // 	}
  //     //      fprintf(stderr, "%lf %lf \n", k[kk], log10(Pk[kk]));
  //   }

  gsl_spline_free(sp);
  gsl_interp_accel_free(acc);

  gsl_spline_free(sp_err);
  gsl_interp_accel_free(acc_err);

  for (int kk = 0; kk < nk; kk++) {
    if (kk < n) {
      k[kk] = k_h[kk];
      Pk[kk] = Pk_h[kk];
      Pk_err[kk] = err_h[kk];
    } else {
      k[kk] = 0.;
      Pk[kk] = 0.;
      Pk_err[kk] = 0;
    }
  // fprintf(stderr, "%lf %lf \n", k[kk], log10(Pk[kk]));
  }

  free(k_h); free(Pk_h); free(err_h); free(any);
  fclose(fp);   fclose(fp_test);
  return(n);
}

int process_PM_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *Pk_err) {
  // Convert to Mpc (no h)

  int ncol = 4;
  int nread = ncol;

  double junk1;
  double *k_h = (double *) malloc(3000*sizeof(double));
  double *Pk_h = (double *) malloc(3000*sizeof(double));
  double *err_h = (double *) malloc(3000*sizeof(double));

  char filename[256], filename_base[256], buf[256];
  sprintf(filename, "/Users/jkwan/emulators/MiraU/matter_power/M%03d/L1300/PM000/analysis/Pow/m%03d.pk.%3d", model_no, model_no, step_no);

  FILE *fp_test = fopen(filename,"r");

  if (fp_test == NULL)
    sprintf(filename_base, "/Users/jkwan/emulators/TitanU/matter_power");
  else
    sprintf(filename_base, "/Users/jkwan/emulators/MiraU/matter_power");

  int n = 0;
  for (int pm_no = 0; pm_no < 16; pm_no++) {
    n = 0;
    sprintf(filename, "%s/M%03d/L1300/PM%03d/analysis/Pow/m%03d.pk.%3d", filename_base, model_no, pm_no, model_no, step_no);
    if (model_no == 0) {
      sprintf(filename, "/Users/jkwan/emulators/M000/PM%03d/analysis/Pow/m%03d.pk.%3d", pm_no, model_no, step_no);
      ncol=3;
    }
    FILE *fp = fopen(filename, "r");
    // fprintf(stderr, "%s\n", filename);

    if (fp == NULL) {
      fprintf(stderr, "Couldn't read %s\n", filename);
      return(-1);
    } else {
      fgets (buf, 256, fp);
    }

    nread = ncol;
    while (nread == ncol) {
      if (ncol==4)
        nread = fscanf(fp, "%lf %lf %lf %lf", &(k_h[n]), &(Pk_h[n]), &junk1, &(err_h[n]));
      else
        nread = fscanf(fp, "%lf %lf %lf", &(k_h[n]), &(Pk_h[n]), &(err_h[n]));
      k_h[n] *= h;
      Pk_h[n] /= pow(h,3);

      err_h[n] = Pk_h[n]/sqrt(err_h[n]);
      // err_h[n] /= pow(h,3);
      //	  fprintf(stderr, "%lf %lf %lf\n", k_h[n], Pk_h[n], err_h[n]);
      n++;
    }
    n--;

    // The following code does the spline for P(k) but I took it out for Earl and Kelly
    // To get all the k's on the same value, uncomment the part below

    // gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,n);
    // gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // gsl_spline *sp_err = gsl_spline_alloc(gsl_interp_cspline,n);
    // gsl_interp_accel *acc_err = gsl_interp_accel_alloc();

    // gsl_spline_init(sp, k_h, Pk_h, n);
    // gsl_spline_init(sp_err, k_h, err_h, n);

    // int kk;
    // for (kk = 0; kk < nk; kk++)
    // 	{
    // 	  if (k[kk] > k_h[0] && k[kk] < k_h[n-1])
    // 	    {
    // 	      Pk[16*kk+pm_no] = gsl_spline_eval(sp, k[kk], acc);
    // 	      Pk_err[16*kk+pm_no] = gsl_spline_eval(sp_err, k[kk], acc_err);
    // 	    }
    // 	  else
    // 	    {
    // 	      Pk[16*kk+pm_no] = 0.;
    // 	      Pk_err[16*kk+pm_no] = 0;
    // 	    }

    // 	  //fprintf(stderr, "%lf %lf\n", k[kk], Pk[kk]);
    // 	}

    // gsl_spline_free(sp);       gsl_spline_free(sp_err);
    // gsl_interp_accel_free(acc);       gsl_interp_accel_free(acc_err);
    fclose(fp);

    for (int kk = 0; kk < nk; kk++) {
      if (kk < n) {
        k[kk] = k_h[kk];
        Pk[16*kk+pm_no] = Pk_h[kk];
        Pk_err[16*kk+pm_no] = err_h[kk];
      } else {
        k[kk] = 0.;
        Pk[16*kk+pm_no] = 0.;
        Pk_err[16*kk+pm_no] = 0;
      }
      // fprintf(stderr, "%lf %lf\n", k[kk], Pk[kk]);
    }
  }
  fclose(fp_test);
  free(k_h); free(Pk_h); free(err_h);
  return(n);
}

int process_PT_runs(double h, int model_no, int step_no, int nk, double *k, double *Pk, double *D, double *Pk_nu) {
  char sed_call[256];
  sprintf(sed_call, "sed '/^#/ d' /Users/jkwan/emulators/PT/redTime_M%03d.dat > junk.dat", model_no);
  /* sprintf(sed_call, "sed '/^#/ d' /Users/jkwan/emulators/PT/redTime_E%03d.dat > junk.dat", model_no); */

   if (model_no == 0)
     sprintf(sed_call, "sed '/^#/ d' /Users/jkwan/emulators/M%03d/redTime_M%03d.dat > junk.dat", model_no, model_no);

  system(sed_call);
  //  FILE *fp = fopen(filename, "r");
  FILE *fp = fopen("junk.dat", "r");

  int nread = 1, nk_pt=202, nz=26;

  /* comment this part for E0XX models */
  if (model_no > 36 || model_no == 0)
    nz = 27;

  double any[15];
  size_t nchar = 1000;
  char *line  = (char *)malloc((nchar+1)*sizeof(char));

  double *k_h = (double *) malloc(nk_pt*sizeof(double));
  double *Pk_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *Pk_nu_h = (double *) malloc(nz*nk_pt*sizeof(double));
  double *D_h = (double *) malloc(nz*nk_pt*sizeof(double));

  int z_no;
  int steps[8] = {163, 189, 247, 300, 347, 401, 453, 499};
  int output_z[8] = {2, 4, 7, 11, 17, 21, 24, 25};

  // need to comment out for extra models
  if (model_no > 36 || model_no == 0)  {
      output_z[4] = 18;
      output_z[5] = 22;
      output_z[6] = 25;
      output_z[7] = 26;
  }

  for (int i = 0; i < 8; i++)
    if (step_no == steps[i])
      z_no = i;

  fprintf(stderr, "%d\n", z_no);

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

  // comment this out for extra models E0XX, leave uncommented for MXXX models:
  // I took this out when Earl and Kelly asked for unsplined P(k)
  // To get all the k's on the same value, uncomment the part below
  // The naming convention changes between model numbers so it is more complicated than it needs to be

  // if (steps[z_no] == 300 && model_no < 37 && model_no!=0)
  //   {
  //     double *Pk_interp = malloc(nk*sizeof(double));
  //     double *D_interp  = malloc(nk*sizeof(double));
  //     double *nu_interp  = malloc(nk*sizeof(double));
  //     double dx, dy, dz, dD, dnu;
  //     dz = 0.0129;
  //     dx = 0.038;
  //     for (kk = 0; kk < nk_pt; kk++)
  //       {
  //         /* double dz = 0.6431-z[11]; */
  //         /* double dx  = z[12]-z[11]; */
  //         dy = Pk_h[nk_pt*12+kk]-Pk_h[nk_pt*11+kk];
  //         Pk_interp[kk] = dy/dx * dz + Pk_h[nk_pt*11+kk];

  // 	  dD = D_h[nk_pt*12+kk]-D_h[nk_pt*11+kk];
  // 	  D_interp[kk] = dD/dx * dz + D_h[nk_pt*11+kk];

  //         dy = Pk_nu_h[nk_pt*12+kk]-Pk_nu_h[nk_pt*11+kk];
  // 	  nu_interp[kk] = dy/dx * dz + Pk_nu_h[nk_pt*11+kk];
  // 	  //          nu_interp[kk] = Pk_nu_h[nk_pt*12+kk];

  //       }

  //     gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
  //     gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline, nk_pt);

  //     gsl_spline_init (sp, k_h, Pk_interp, nk_pt);

  //     gsl_interp_accel *accD  = gsl_interp_accel_alloc ();
  //     gsl_spline *spD = gsl_spline_alloc(gsl_interp_cspline, nk_pt);

  //     gsl_spline_init (spD, k_h, D_interp, nk_pt);

  //     double D0 = gsl_spline_eval(spD, k[nk_pt-1], accD);


  //     gsl_interp_accel *acc_nu  = gsl_interp_accel_alloc ();
  //     gsl_spline *sp_nu = gsl_spline_alloc(gsl_interp_cspline, nk_pt);

  //     gsl_spline_init (sp_nu, k_h, nu_interp, nk_pt);

  //     for (kk = 0; kk < nk; kk++)
  //       {

  // 	  if (k[kk] > k_h[0] && k[kk] < k_h[nk_pt-1])
  // 	    {
  // 	      Pk[kk] = gsl_spline_eval(sp, k[kk], acc);
  // 	      D[kk] = gsl_spline_eval(spD, k[kk], accD)/D0;
  // 	      Pk_nu[kk] = gsl_spline_eval(sp_nu, k[kk], acc_nu);
  // 	      //	      fprintf(stderr, "%d %lf\n", kk, k[kk]);
  // 	    }
  // 	  else
  // 	    {
  // 	      Pk[kk] = 1.;
  // 	      D[kk] = D[kk-1];
  // 	      Pk_nu[kk] = 1.;
  // 	      //	      fprintf(stderr, "%lf\n", k[kk]);
  // 	    }
  // 	}

  //     free(Pk_interp); free(D_interp); free(nu_interp);
  //     gsl_spline_free(sp); gsl_interp_accel_free(acc);
  //     gsl_spline_free(spD); gsl_interp_accel_free(accD);
  //     gsl_spline_free(sp_nu); gsl_interp_accel_free(acc_nu);

  //   }
  // else
    // {
    //   gsl_spline *sp = gsl_spline_alloc(gsl_interp_cspline,nk_pt);
    //   gsl_interp_accel *acc = gsl_interp_accel_alloc();

    //   gsl_spline_init(sp, k_h, &(Pk_h[nk_pt*output_z[z_no]]), nk_pt);

    //   gsl_spline *spD = gsl_spline_alloc(gsl_interp_cspline,nk_pt);
    //   gsl_interp_accel *accD = gsl_interp_accel_alloc();

    //   gsl_spline_init(spD, k_h, &(D_h[nk_pt*output_z[z_no]]), nk_pt);

    //   double D0 = gsl_spline_eval(spD, k_h[nk_pt-1], accD);

    //   //      fprintf(stderr, "knorm = %lf %lf\n", k_h[nk_pt-1], D0);

    //   gsl_interp_accel *acc_nu  = gsl_interp_accel_alloc ();
    //   gsl_spline *sp_nu = gsl_spline_alloc(gsl_interp_cspline, nk_pt);

    //   gsl_spline_init (sp_nu, k_h, &(Pk_nu_h[nk_pt*output_z[z_no]]), nk_pt);


    //   for (kk = 0; kk < nk; kk++)
    // 	if (k[kk] < k_h[nk_pt-1] && k[kk] > k_h[0])
    // 	  {
    // 	    Pk[kk] = gsl_spline_eval(sp, k[kk], acc);

    // 	    D[kk] = gsl_spline_eval(spD, k[kk], accD)/D0;
    // 	    //	    fprintf(stderr, "%lf %lf\n", k[kk], D[kk]);
    // 	    Pk_nu[kk] = gsl_spline_eval(sp_nu, k[kk], acc_nu);

    // 	  }
    // 	else
    // 	  {
    // 	    Pk[kk] = 1.;
    // 	    D[kk] = 1.;//D[kk-1];
    // 	    Pk_nu[kk] = 1.;
    // 	  }

    //   gsl_spline_free(sp); gsl_interp_accel_free(acc);
    //   gsl_spline_free(spD); gsl_interp_accel_free(accD);
    //   gsl_spline_free(sp_nu); gsl_interp_accel_free(acc_nu);
    // }

  double D0 = D_h[nk_pt*(output_z[z_no]+1)-1];
  //  fprintf(stderr, "%lf\n", D0);

  for (int kk=0; kk < nk; kk++) {
    if (kk < nk_pt) {
      k[kk] = k_h[kk];
      Pk[kk] = Pk_h[nk_pt*output_z[z_no]+kk];
      D[kk] = D_h[nk_pt*output_z[z_no]+kk]/D0;
      Pk_nu[kk] = Pk_nu_h[nk_pt*output_z[z_no]+kk];
      //	    fprintf(stderr, "%lf\n", D_h[nk_pt*output_z[z_no]+kk]/D0);
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
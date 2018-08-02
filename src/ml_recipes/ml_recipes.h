/*
 *  ml_recipes.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef ML_RECIPES_H
#define ML_RECIPES_H

//The stuff in this directory is a fairly light port of the godawful bloody CC code in ARTIST's modellib. I didn't dare change too much because of the proliferation of meaningless, ueber-short variable names. Way to go for incomprehensible code. :-/

#include "../ml_types.h"
#include "../constants.h"
#include "../ml_funcs.h"
#include "../ml_models.h"
#include "../py_utils.h"
#include "../error_codes.h"

#ifdef WITH_C99
#include <complex.h>
#endif

#define ML_MEAN_MOL_WT  2.3
//#define ML_oneOnMuMp    (1.0/ML_MEAN_MOL_WT/MPROTON)
#define ML_oneOnMuMp    (1.0/ML_MEAN_MOL_WT/AMU)
#define ML_gsize     1e-5    /* I have no clue. */

#define BoEb56_q       1.03    /* Relative stepsize in the integration (r_i+1 / r_i) */

#define CG97_npgrid    300

#define DDN01_nwav    100
#define DDN01_nrmod    100

#define Me09_ngrid    100
#define Me09_mmw      2.0

#define Ul76_mmw      2.0


int	Al03_onFinalizeConfiguration(void);
double	Al03_density(const double x, const double y, const double z);
double	Al03_temperature(const double x, const double y, const double z);
void	Al03_velocity(const double x, const double y, const double z, double *vel);
void	Al03_bmag(const double x, const double y, const double z, double *b);

int	BoEb56_onFinalizeConfiguration(void);
double	BoEb56_density(double x, double y, double z);
double	BoEb56_temperature(double x, double y, double z);
void	BoEb56_velocity(double x, double y, double z, double *vel);

int	CG97_onFinalizeConfiguration(void);
double	CG97_density(const double x, const double y, const double z);
double	CG97_t_dust(const double x, const double y, const double z);
void	CG97_velocity(const double x, const double y, const double z, double *vel);

int	DDN01_onFinalizeConfiguration(void);
double	DDN01_t_dust(const double x, const double y, const double z);
double	DDN01_density(const double x, const double y, const double z);
void	DDN01_velocity(const double x, const double y, const double z, double* v);

int	LiSh96_onFinalizeConfiguration(void);
double	LiSh96_density(const double x, const double y, const double z);
void	LiSh96_bmag(const double x, const double y, const double z, double* B);
void	LiSh96_velocity(const double x, const double y, const double z, double* v);
double	LiSh96_temperature(const double x, const double y, const double z);

int	Ma88_onFinalizeConfiguration(void);
double	Ma88_density(const double x, const double y, const double z);
double	Ma88_temperature(const double x, const double y, const double z);
double	Ma88_abundance(const double x, const double y, const double z);
void	Ma88_velocity(const double x, const double y, const double z, double* v);

int	Me09_onFinalizeConfiguration(void);
double	Me09_density(const double x, const double y, const double z);
void	Me09_velocity(const double x, const double y, const double z, double* v);

int 	Shu77_onFinalizeConfiguration(void);
double	Shu77_density(const double x, const double y, const double z);
double	Shu77_temperature(const double x, const double y, const double z);
double	Shu77_t_dust(const double x, const double y, const double z);
void	Shu77_velocity(const double x, const double y, const double z, double* v);

int 	Ul76_onFinalizeConfiguration(void);
double	Ul76_density(const double x, const double y, const double z);
void	Ul76_velocity(const double x, const double y, const double z, double* v);

#endif /* ML_RECIPES_H */


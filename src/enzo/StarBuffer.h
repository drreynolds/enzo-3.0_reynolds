/***********************************************************************
/
/  STAR PARTICLE STRUCTURE FOR COMMUNICATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef __STARBUFFER_H
#define __STARBUFFER_H

#define MAX_ACCR 100

struct StarBuffer {
  FLOAT	pos[MAX_DIMENSION];
  float		vel[MAX_DIMENSION];
  float		delta_vel[MAX_DIMENSION];
  int		naccretions;
  float		accretion_rate[MAX_ACCR];
  FLOAT	accretion_time[MAX_ACCR];
  float		Mass;
  float		FinalMass;
  float		DeltaMass;
  float		BirthTime;
  float		LifeTime;
  int		FeedbackFlag;
  int		Identifier;
  int		level;
  int		GridID;
  star_type	type;
};

#endif

/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file cl-timing.h
 * \brief Note that this file is intially written by Andreas Kloeckner and
 * further changed by Christian Eder.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GBLA_CL_TIMING_H
#define GBLA_CL_TIMING_H
#ifdef __APPLE__

#include <sys/time.h>

typedef struct timeval timestamp_type;

static inline void get_timestamp(timestamp_type *t)
{
  gettimeofday(t, NULL);
}

static inline double timestamp_diff_in_seconds(timestamp_type start,
timestamp_type end)
{
  /* Perform the carry for the later subtraction by updating start. */
  if (end.tv_usec < start.tv_usec) {
    int nsec = (start.tv_usec - end.tv_usec) / 1000000 + 1;
    start.tv_usec -= 1000000 * nsec;
    start.tv_sec += nsec;
  }
  if (end.tv_usec - start.tv_usec > 1000000) {
    int nsec = (end.tv_usec - start.tv_usec) / 1000000;
    start.tv_usec += 1000000 * nsec;
    start.tv_sec -= nsec;
  }

  return end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6;
}

#else

#include <time.h>

typedef struct timespec timestamp_type;

static inline void get_timestamp(timestamp_type *t)
{
  clock_gettime(CLOCK_REALTIME, t);
}

static inline double timestamp_diff_in_seconds(timestamp_type start, timestamp_type end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp.tv_sec + 1e-9*temp.tv_nsec;
}

#endif
#endif

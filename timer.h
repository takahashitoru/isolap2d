#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include "elapsed.h"

typedef struct {
  double time;
  double work;
} timerType;

#ifdef __cplusplus
extern "C" {
#endif
  void allocTimer(timerType **timer);
  void initTimer(timerType *timer);
  void startTimer(timerType *timer);
  void stopTimer(timerType *timer);
  void freeTimer(timerType **timer);
  void printTimer(FILE *fp, char *s, timerType *timer);
  double getTimer(timerType timer);
#ifdef __cplusplus
}
#endif

#endif /* TIMER_H */

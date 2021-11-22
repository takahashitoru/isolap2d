#include "timer.h"

void allocTimer(timerType **timer)
{
  *timer = (timerType *)malloc(sizeof(timerType));
  if (*timer == NULL) {
    fprintf(stderr, "%s: fail to allocate. Exit.\n", __FUNCTION__);
    exit(1);
  }
}

void initTimer(timerType *timer)
{
  timer->time = 0.0;
}  

void startTimer(timerType *timer)
{
  timer->work = elapsed();
}

void stopTimer(timerType *timer)
{
  timer->time += elapsed() - timer->work;
}

void freeTimer(timerType **timer)
{
  free(*timer);
}

void printTimer(FILE *fp, char *s, timerType *timer)
{
  fprintf(fp, "# %s: %s = %f\n", __FUNCTION__, s, timer->time);
}

double getTimer(timerType timer)
{
  return(timer.time);
}

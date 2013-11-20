#include <stdlib.h>
#include <string.h>
#include "timer.h"
#include "ts_util.h"

void timer_init(timer_ctx_t* timer)
{
  memset(timer, 0, sizeof(*timer));
}

void timer_start(timer_ctx_t* timer)
{
  timer->start = tsTOD();
}

void timer_stop(timer_ctx_t* timer)
{
  timer->end = tsTOD();
}

char *timer_diff_as_str(timer_ctx_t* timer)
{
  struct timespec diff = tsSubtract(timer->end, timer->start);

  return tsShow(diff, false, "min %M, sec %S");
}

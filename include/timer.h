#ifndef TIMER_H
#define TIMER_H 1

typedef struct timer_ctx_s
{
  struct timespec start;
  struct timespec end;
}
timer_ctx_t;

void timer_init(timer_ctx_t* timer);
void timer_start(timer_ctx_t* timer);
void timer_stop(timer_ctx_t* timer);
char *timer_diff_as_str(timer_ctx_t* timer);

#endif /* TIMER_H */

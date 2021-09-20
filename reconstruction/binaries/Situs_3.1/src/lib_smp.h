/* header file for lib_smp.c */

/*
 * workq.h
 *
 * This header file defines the interfaces for a "work queue"
 * manager. A "manager object" is created with several
 * parameters, including the required size of a work queue
 * entry, the maximum desired degree of parallelism (number of
 * threads to service the queue), and the address of an
 * execution engine routine.
 *
 * The application requests a work queue entry from the manager,
 * fills in the application-specific fields, and returns it to
 * the queue manager for processing. The manager will create a
 * new thread to service the queue if all current threads are
 * busy and the maximum level of parallelism has not yet been
 * reached.
 *
 * The manager will dequeue items and present them to the
 * processing engine until the queue is empty; at that point,
 * processing threads will begin to shut down. (They will be
 * restarted when work appears.)
 */

#ifndef __SITUS_LIB_SMP
#define __SITUS_LIB_SMP

#include <stdarg.h>
#include <errno.h>

#ifdef _SMP_
#include <pthread.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Powell_Args {
  unsigned long   iter;
  double          *pow_init;
#ifdef _SMP_
  pthread_mutex_t print_mutex;
#endif
} POWARG;

#ifdef _SMP_

/*
 * Structure to keep track of work queue requests.
 */
typedef struct workq_ele_tag {
  struct workq_ele_tag   *next;
  void (*engine)(void *arg);      /* user engine */
  void   *data;
} workq_ele_t;

/*
 * Structure describing a work queue.
 */
typedef struct workq_tag {
  pthread_mutex_t     mutex;
  pthread_cond_t      cv;             /* wait for work */
  pthread_mutex_t     barriermutex;
  pthread_cond_t      barrier;        /* barrier */
  pthread_attr_t      attr;           /* create detached threads */
  workq_ele_t         *first, *last;  /* work queue */
  int                 valid;          /* set when valid */
  int                 quit;           /* set when workq should quit */
  int                 parallelism;    /* number of threads required */
  int                 counter;        /* current number of threads */
  int                 idle;           /* number of idle threads */
} workq_t;

#define WORKQ_VALID     0xdec1992

/*
 * Define work queue functions
 */
extern workq_t *workq_init(int threads);
extern void *workq_server(void *args);
extern void workq_barrier(workq_t *wq);
extern void workq_destroy(workq_t *wq);
extern void workq_add(workq_t *wq, void (*engine)(void *arg), void *data);

#endif

/* Library internal utility functions */

inline void print(char *fmt, ...);
inline void error(char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif

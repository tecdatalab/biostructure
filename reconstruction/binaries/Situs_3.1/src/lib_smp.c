/*********************************************************************
*                           L I B _ S M P                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Jochen Heyd, Julio Kovacs and Willy Wriggers 2005-2015         *
**********************************************************************
*                                                                    *
* SMP thread functions for parallel colores                          *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_smp.h"

#ifdef _SMP_

/*
 * Initialize a work queue.
 */
workq_t *workq_init(int threads)
{
  workq_t *wq;
  int status;

  wq = (workq_t *) alloc_vect(sizeof(workq_t), 1);

  status = pthread_attr_init(&wq->attr);
  if (status != 0) {
    error("Pthread_attr_init failed\n");
  }
  status = pthread_attr_setdetachstate(
             &wq->attr, PTHREAD_CREATE_DETACHED);
  if (status != 0) {
    pthread_attr_destroy(&wq->attr);
    error("Pthread_attr_setdetachstate failed\n");
  }
  status = pthread_mutex_init(&wq->mutex, NULL);
  if (status != 0) {
    pthread_attr_destroy(&wq->attr);
    error("Pthread_mutex_init failed\n");
  }
  status = pthread_mutex_init(&wq->barriermutex, NULL);
  if (status != 0) {
    pthread_attr_destroy(&wq->attr);
    error("Pthread_mutex_init failed\n");
  }
  status = pthread_cond_init(&wq->cv, NULL);
  if (status != 0) {
    pthread_mutex_destroy(&wq->mutex);
    pthread_attr_destroy(&wq->attr);
    error("Pthread_cond_init failed\n");
  }
  status = pthread_cond_init(&wq->barrier, NULL);
  if (status != 0) {
    pthread_mutex_destroy(&wq->mutex);
    pthread_attr_destroy(&wq->attr);
    error("Pthread_cond_init failed\n");
  }
  wq->quit = 0;                       /* not time to quit */
  wq->first = wq->last = NULL;        /* no queue entries */
  wq->parallelism = threads;          /* max servers */
  wq->counter = 0;                    /* no server threads yet */
  wq->idle = 0;                       /* no idle servers */
  wq->valid = WORKQ_VALID;
  return wq;
}


/*
 * Thread start routine to serve the work queue.
 */
void *workq_server(void *arg)
{
  struct timeval timeout_sys;
  struct timespec timeout;
  workq_t *wq = (workq_t *)arg;
  workq_ele_t *we;
  int status, timedout;

  /*
   * We don't need to validate the workq_t here... we don't
   * create server threads until requests are queued (the
   * queue has been initialized by then!) and we wait for all
   * server threads to terminate before destroying a work
   * queue.
   */
  status = pthread_mutex_lock(&wq->mutex);
  if (status != 0)
    return NULL;

  while (1) {
    timedout = 0;
    gettimeofday(&timeout_sys, NULL);
    timeout.tv_sec  = timeout_sys.tv_sec + 2;
    timeout.tv_nsec = 0;

    while (wq->first == NULL && !wq->quit) {
      /*
       * Server threads time out after spending 2 seconds
       * waiting for new work, and exit.
       */
      wq->idle++;
      status = pthread_cond_timedwait(
                 &wq->cv, &wq->mutex, &timeout);
      wq->idle--;
      if (status == ETIMEDOUT) {
        timedout = 1;
        break;
      } else if (status != 0) {
        /*
         * This shouldn't happen, so the work queue
         * package should fail. Because the work queue
         * API is asynchronous, that would add
         * complication. Because the chances of failure
         * are slim, I choose to avoid that
         * complication. The server thread will return,
         * and allow another server thread to pick up
         * the work later. Note that, if this was the
         * only server thread, the queue won't be
         * serviced until a new work item is
         * queued. That could be fixed by creating a new
         * server here.
         */
        print(
          "WorkQ: Worker wait failed, %d (%s)\n",
          status, strerror(status));
        wq->counter--;
        pthread_mutex_unlock(&wq->mutex);
        return NULL;
      }
    }

    // printf ("Work queue: %#lx, counter: %d, idle: %d, quit: %d\n",
    //         wq->first, wq->counter, wq->idle, wq->quit);
    we = wq->first;

    if (we != NULL) {
      wq->first = we->next;
      if (wq->last == we)
        wq->last = NULL;
      status = pthread_mutex_unlock(&wq->mutex);
      if (status != 0)
        return NULL;
      we->engine(we->data);
      free_vect_and_zero_ptr(&we);
      status = pthread_mutex_lock(&wq->mutex);
      if (status != 0)
        return NULL;
    }

    /* Broadcast barrier if we hit it */
    if (wq->first == NULL) {
      pthread_cond_broadcast(&wq->barrier);
    }

    /*
     * If there are no more work requests, and the servers
     * have been asked to quit, then shut down.
     */
    if (wq->first == NULL && wq->quit) {
      wq->counter--;


      /*
       * NOTE: Just to prove that every rule has an
       * exception, I'm using the "cv" condition for two
       * separate predicates here.  That's OK, since the
       * case used here applies only once during the life
       * of a work queue -- during rundown. The overhead
       * is minimal and it's not worth creating a separate
       * condition variable that would be waited and
       * signaled exactly once!
       */
      if (wq->counter == 0)
        pthread_cond_broadcast(&wq->cv);
      pthread_mutex_unlock(&wq->mutex);
      return NULL;
    }

    /*
     * If there's no more work, and we wait for as long as
     * we're allowed, then terminate this server thread.
     */
    if (wq->first == NULL && timedout) {
      wq->counter--;
      break;
    }
  }

  pthread_mutex_unlock(&wq->mutex);
  return NULL;
}


/*
 * Add an item to a work queue.
 */
void workq_add(workq_t *wq, void (*engine)(void *arg), void *element)
{
  workq_ele_t *item;
  pthread_t id;
  int status;

  if (wq->valid != WORKQ_VALID)
    error("WorkQ Add: workq invalid\n");

  /*
   * Create and initialize a request structure.
   */
  item = (workq_ele_t *) alloc_vect(sizeof(workq_ele_t), 1);
  item->engine = engine;
  item->data   = element;
  item->next   = NULL;
  status = pthread_mutex_lock(&wq->mutex);
  if (status != 0) {
    free_vect_and_zero_ptr(&item);
    error("WorkQ Add: Mutex lock failed\n");
  }

  /*
   * Add the request to the end of the queue, updating the
   * first and last pointers.
   */
  if (wq->first == NULL)
    wq->first = item;
  else
    wq->last->next = item;
  wq->last = item;

  /*
   * if any threads are idling, wake one.
   */
  if (wq->idle > 0) {
    status = pthread_cond_signal(&wq->cv);
    if (status != 0) {
      pthread_mutex_unlock(&wq->mutex);
      error("WorkQ Add: Waking threads failed\n");
    }
  } else if (wq->counter < wq->parallelism) {
    /*
     * If there were no idling threads, and we're allowed to
     * create a new thread, do so.
     */
    status = pthread_create(
               &id, &wq->attr, workq_server, (void *)wq);
    if (status != 0) {
      pthread_mutex_unlock(&wq->mutex);
      error("WorkQ Add: Creating thread failed\n");
    }
    wq->counter++;
  }
  pthread_mutex_unlock(&wq->mutex);
  return;
}

/*
 * Barrier for work queue
 */
void workq_barrier(workq_t *wq)
{
  int status;

  if (wq->valid != WORKQ_VALID)
    error("WorkQ Barrier: workq invalid\n");
  status = pthread_mutex_lock(&wq->barriermutex);
  if (status != 0)
    error("WorkQ Barrier: Mutex lock failed\n");

  if (wq->first != NULL) {
    while (wq->first != NULL) {
      status = pthread_cond_wait(&wq->barrier, &wq->barriermutex);
      if (status != 0) {
        pthread_mutex_unlock(&wq->barriermutex);
        error("WorkQ Barrier: Conditional wait failed\n");
      }
    }
  }

  pthread_mutex_unlock(&wq->barriermutex);

  return;
}


/*
 * Destroy a work queue.
 */
void workq_destroy(workq_t *wq)
{
  int status;

  if (wq->valid != WORKQ_VALID)
    error("WorkQ Destroy: workq invalid\n");
  status = pthread_mutex_lock(&wq->mutex);
  if (status != 0)
    error("WorkQ Destroy: First mutex lock failed\n");
  wq->valid = 0;                 /* prevent any other operations */

  /*
   * Check whether any threads are active, and run them down:
   *
   * 1.       set the quit flag
   * 2.       broadcast to wake any servers that may be asleep
   * 4.       wait for all threads to quit (counter goes to 0)
   *          Because we don't use join, we don't need to worry
   *          about tracking thread IDs.
   */
  if (wq->counter > 0) {
    wq->quit = 1;
    /* if any threads are idling, wake them. */
    if (wq->idle > 0) {
      status = pthread_cond_broadcast(&wq->cv);
      if (status != 0) {
        pthread_mutex_unlock(&wq->mutex);
        error("WorkQ Destroy: Waking threads failed\n");
      }
    }

    /*
     * Just to prove that every rule has an exception, I'm
     * using the "cv" condition for two separate predicates
     * here. That's OK, since the case used here applies
     * only once during the life of a work queue -- during
     * rundown. The overhead is minimal and it's not worth
     * creating a separate condition variable that would be
     * waited and signalled exactly once!
     */
    while (wq->counter > 0) {
      status = pthread_cond_wait(&wq->cv, &wq->mutex);
      if (status != 0) {
        pthread_mutex_unlock(&wq->mutex);
        error("WorkQ Destroy: Waiting for threads to finish failed\n");
      }
    }
  }
  status = pthread_mutex_unlock(&wq->mutex);
  if (status != 0)
    error("WorkQ Destroy: Mutex unlock failed\n");
  status = pthread_mutex_destroy(&wq->mutex);
  if (status != 0)
    error("WorkQ Destroy: Mutex destroy failed\n");
  status = pthread_cond_destroy(&wq->cv);
  if (status != 0)
    error("WorkQ Destroy: Mutex cond destroy failed\n");
  status = pthread_cond_destroy(&wq->barrier);
  if (status != 0)
    error("WorkQ Destroy: Barrier cond destroy failed\n");
  status = pthread_attr_destroy(&wq->attr);
  if (status != 0)
    error("WorkQ Destroy: Attribute destroy failed\n");

  free_vect_and_zero_ptr(&wq);

  return;
}
#endif

/* error: print error message and exit */

extern inline void error(char *fmt, ...)
{
  va_list args;

  fflush(stdout);

  fprintf(stderr, "lib_smp> ");
  fprintf(stderr, "ERROR - ");

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  if (fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
    fprintf(stderr, " %s", strerror(errno));
  fprintf(stderr, "\n");
  exit(2); /* conventional value for failed execution */
}

/* print: printf which flushes buffer */

extern inline void print(char *fmt, ...)
{
  va_list args;

  printf("lib_smp> ");
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);

  fflush(stdout);

  return;
}

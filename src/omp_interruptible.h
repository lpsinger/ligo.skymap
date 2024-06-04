/*
 * Copyright (C) 2017-2024  Leo Singer
 *
 * These preprocessor macros help make long-running Python C extensions,
 * possibly that contain OpenMP parallel for loops, respond gracefully to
 * signals.
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


/*
 * In a normal C program that does not mess with signals, when the user
 * types Ctrl-C, the process is sent the SIGINT signal. The default SIGINT
 * signal handler terminates the process swiftly, even if the program is in
 * the middle of a CPU-intensive loop, even an OpenMP parallel for loop.
 *
 * It's different for a Python program. Python itself attaches handlers for
 * most (all?) signals. Python's SIGINT handler sets a flag to remind itself to
 * raise a KeyboardInterrupt exception on the main thread before interpreting
 * the next instruction.
 *
 * If Python is in the middle of executing a long-running method in a Python C
 * extension, then the interpreter will remain unresponsive until the method
 * has returned, when it can raise the KeyboardInterrupt. This delay can be
 * very annoying to the user.
 *
 * This header provides a few macros to temporarily change the SIGINT handler
 * and provide a flag that C functions can check periodically to terminate
 * early. Here's a skeleton code sample that includes an OpenMP loop to show
 * how to use the macros. We start with our basic function, foo, which contains
 * an OpenMP parallel for loop:
 *
 *      int foo(int n)
 *      {
 *          int retval = 0;
 *          #pragma omp parallel for
 *          for (int i = 0; i < n; i ++)
 *          {
 *              ... // The actual work occurs here.
 *          }
 *          return retval;
 *      }
 *
 * We add the macros:
 *
 *      #include "omp_interruptible.h"
 *
 *      // Native C function that does the work.
 *      int foo(int n)
 *      {
 *          int retval = 0;
 *          OMP_BEGIN_INTERRUPTIBLE  // Replace SIGINT handler.
 *          #pragma omp parallel for
 *          for (int i = 0; i < n; i ++)
 *          {
 *              // Exit loop early if SIGINT has fired.
 *              // Note: you can replace OMP_EXIT_LOOP_EARLY with a simple
 *              // `break;` statement or check it in the loop conditional,
 *              // if you are not using an OpenMP loop.
 *              if (OMP_WAS_INTERRUPTED)
 *                  OMP_EXIT_LOOP_EARLY
 *              ...  // The actual work occurs here.
 *          }
 *          if (OMP_WAS_INTERRUPTED)
 *              retval = -1;
 *          OMP_END_INTERRUPTIBLE  // Restore SIGINT handler.
 *          return retval;
 *      }
 *
 * Finally, here's the Python C extension:
 *
 *      #include <Python.h>
 *
 *      static PyObject *mod_foo(PyObject *module, PyObject *args)
 *      {
 *          int reval;
 *
 *          // Run the underlying C function, releasing the global interpreter
 *          // lock (GIL) in the mean time so that other Python threads (if
 *          // any) can run.
 *          Py_BEGIN_ALLOW_THREADS
 *          int retval = foo(1000);
 *          Py_END_ALLOW_THREADS
 *
 *          // Important: call PyErr_CheckSignals() to give Python a chance to
 *          // raise a KeyboardInterrupt exception, if needed.
 *
 *          // Indicate success or failure of the method to the interpreter.
 *          PyErr_CheckSignals();
 *          if (retval == 0)
 *              Py_RETURN_NONE;
 *          else
 *              return NULL;
 *      }
 *
 *      static PyMethodDef methods[] = {
 *          {"foo", (PyCFunction)mod_foo, METH_NOARGS, "doc string here"},
 *          {NULL, NULL, 0, NULL}
 *      };
 *
 *       static PyModuleDef moduledef = {
 *           PyModuleDef_HEAD_INIT,
 *           "mod", NULL, -1, methods,
 *           NULL, NULL, NULL, NULL
 *       };
 *
 *      PyMODINIT_FUNC PyInit_mod(void)
 *      {
 *          return PyModule_Create(&moduledef);
 *      }
 *
 * Note that only one section of a code at a time will have the signal handler
 * active. If multiple threads call OMP_BEGIN_INTERRUPTIBLE at the same time,
 * only one of them is guaranteed to be cancelled.
 */


#ifndef OMP_INTERRUPTIBLE_H
#define OMP_INTERRUPTIBLE_H

#include "branch_prediction.h"

#include <assert.h>
#include <signal.h>
#include <stdlib.h>
#include <pthread.h>


static int omp_was_interrupted = 0;
static const int omp_interruptible_signal = SIGINT;
static pthread_mutex_t omp_interruptible_lock = PTHREAD_MUTEX_INITIALIZER;
static struct sigaction omp_interruptible_old_action;


/* A utility to safely assert that a function that may have side-effects
 * returns a nonzero value.
 *
 * On a non-debug build, `assert(!foo())` will not emit any code, whereas
 * `must_succed(foo())` will compile to the same code as `foo()`. */
static void must_succeed(int result)
{
    assert(!result);
}


static void omp_interruptible_restore_handler()
{
    must_succeed(sigaction(
        omp_interruptible_signal,
        &omp_interruptible_old_action,
        NULL));
}


static void omp_interruptible_handler(int sig)
{
    omp_was_interrupted = 1;
    #pragma omp flush
    omp_interruptible_restore_handler();
    raise(sig);
}


static const struct sigaction omp_interruptible_action = {
    .sa_handler = omp_interruptible_handler
};


static int omp_interruptible_begin() {
    int result = pthread_mutex_trylock(&omp_interruptible_lock);
    if (LIKELY(!result)) {
        omp_was_interrupted = 0;
        must_succeed(sigaction(
            omp_interruptible_signal,
            &omp_interruptible_action,
            &omp_interruptible_old_action));
    }
    return result;
}


static void omp_interruptible_end(int omp_interruptible_begin_result) {
    if (LIKELY(!omp_interruptible_begin_result)) {
        omp_interruptible_restore_handler();
        must_succeed(pthread_mutex_unlock(&omp_interruptible_lock));
    }
}


static int omp_interruptible_get(int omp_interruptible_begin_result) {
    if (UNLIKELY(omp_interruptible_begin_result))
        return 0;
    else
        return omp_was_interrupted;
}


#define OMP_BEGIN_INTERRUPTIBLE int omp_interruptible_begin_result = omp_interruptible_begin();
#define OMP_END_INTERRUPTIBLE omp_interruptible_end(omp_interruptible_begin_result);
#define OMP_WAS_INTERRUPTED UNLIKELY(omp_interruptible_get(omp_interruptible_begin_result))


#if _OPENMP
#define OMP_EXIT_LOOP_EARLY continue;
#else
#define OMP_EXIT_LOOP_EARLY break;
#endif


#endif /* OMP_INTERRUPTIBLE_H */

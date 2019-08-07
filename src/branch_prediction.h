/* GCC branch prediction hints */

#ifndef BRANCH_PREDICTION_H
#define BRANCH_PREDICTION_H

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

#endif /* BRANCH_PREDICTION_H */

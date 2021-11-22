#ifndef DEBUGMACROS_H
#define DEBUGMACROS_H

#define INFO(fmt, ...) fprintf(stderr, "# %s: " fmt, __FUNCTION__, ## __VA_ARGS__)
#define MESG(message) fprintf(stderr, "# %s: %s", __FUNCTION__, message)
#if defined(MYDEBUG)
#define DBG(fmt, ...) fprintf(stderr, "### %s: " fmt, __FUNCTION__, ## __VA_ARGS__)
#define MSG(message) fprintf(stderr, "### %s: %s", __FUNCTION__, message)
#else
#define DBG(fmt, ...)
#define MSG(message)
#endif

#if defined(_DEBUG) || defined(MYDEBUG)
#include <assert.h>
#define ASSERT(s) assert(s)
#else
#define ASSERT(s)
#endif

#endif /* DEBUGMACROS_H */

#ifndef MBA_DIFF_H
#define MBA_DIFF_H

/* diff - compute a shortest edit script (SES) given two sequences
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef LIBMBA_API
#ifdef WIN32
# ifdef LIBMBA_EXPORTS
#  define LIBMBA_API  __declspec(dllexport)
# else /* LIBMBA_EXPORTS */
#  define LIBMBA_API  __declspec(dllimport)
# endif /* LIBMBA_EXPORTS */
#else /* WIN32 */
# define LIBMBA_API extern
#endif /* WIN32 */
#endif /* LIBMBA_API */

#include <mba/varray.h>
#include <mba/hashmap.h> /* cmp_fn */

typedef const void *(*idx_fn)(const void *s, int idx, void *context);

typedef enum {
	DIFF_MATCH = 1,
	DIFF_DELETE,
	DIFF_INSERT
} diff_op;

struct diff_edit {
	short op;
	int off; /* off into s1 if MATCH or DELETE but s2 if INSERT */
	int len;
};

/* consider alternate behavior for each NULL parameter
 */
LIBMBA_API int diff(const void *a, int aoff, int n,
		const void *b, int boff, int m,
		idx_fn idx, cmp_fn cmp, void *context, int dmax,
		struct varray *ses, int *sn,
		struct varray *buf);

#ifdef __cplusplus
}
#endif

#endif /* MBA_DIFF_H */

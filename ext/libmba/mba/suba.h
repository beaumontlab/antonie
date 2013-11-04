#ifndef MBA_SUBA_H
#define MBA_SUBA_H

/* suba - sub-allocate memory from larger chunk of memory
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

#include <stddef.h>
#include <mba/allocator.h>

#define SUBA_PTR_SIZE(ptr) ((ptr) ? (*((size_t *)(ptr) - 1)) : 0)

extern struct allocator *suba_init(void *mem, size_t size, int rst, size_t mincell);
extern void *suba_alloc(struct allocator *suba, size_t size, int zero);
extern void *suba_realloc(struct allocator *suba, void *ptr, size_t size);
extern int suba_free(void *suba, void *ptr);

extern void *suba_addr(const struct allocator *suba, const ref_t ref);
extern ref_t suba_ref(const struct allocator *suba, const void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* MBA_SUBA_H */


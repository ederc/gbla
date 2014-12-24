#ifndef __GB_macros_H
#define __GB_macros_H

#define Mjoin(pre,nam) my_join(pre , nam)
#define my_join(pre, nam) pre ## _ ## nam

#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define min(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b; })

#define SWAP(a,b)  \
	t = a ; \
	a = b ; \
	b = t



#endif /* __GB_macros_H */



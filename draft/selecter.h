#ifndef __GB_selecter_h
#define __GB_selecter_h

uint32_t select_int8_t()   { return 1 | (0<<1) ; }
uint32_t select_uint8_t()  { return 0 | (0<<1) ; }
uint32_t select_int16_t()  { return 1 | (1<<1) ; }
uint32_t select_uint16_t() { return 0 | (1<<1) ; }
uint32_t select_int32_t()  { return 1 | (2<<1) ; }
uint32_t select_uint32_t() { return 0 | (2<<1) ; }
uint32_t select_int64_t()  { return 1 | (3<<1) ; }
uint32_t select_uint64_t() { return 0 | (3<<1) ; }

#endif

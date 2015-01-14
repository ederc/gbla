#ifndef __GB_printer_h
#define __GB_printer_h

void print_double(double m)
{
	fprintf(stderr,"%f",m);
}
void print_float(float m)
{
	fprintf(stderr,"%f",m);
}
void print_int16_t(int16_t m)
{
	fprintf(stderr,"%d",m);
}
void print_uint16_t(uint16_t m)
{
	fprintf(stderr,"%u",m);
}
void print_int32_t(int32_t m)
{
	fprintf(stderr,"%d",m);
}
void print_uint32_t(uint32_t m)
{
	fprintf(stderr,"%u",m);
}
void print_int64_t(int64_t m)
{
	fprintf(stderr,"%ld",m);
}
void print_uint64_t(uint64_t m)
{
	fprintf(stderr,"%lu",m);
}


#endif /* __GB_printer_h */

/* vim: set ft=c: */


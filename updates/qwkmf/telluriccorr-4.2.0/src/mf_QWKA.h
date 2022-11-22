#ifndef QWKA_HEADER_H
#define QWKA_HEADER_H
typedef enum {VERBOSE, COMPARE, TERMINATE, TAPE5} QWKA_FLAG;
cpl_boolean mf_QWKAParseEnvVar(QWKA_FLAG FLAG);
#endif


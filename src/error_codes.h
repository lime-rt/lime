/*
 *  error_codes.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef ERROR_CODES_H
#define ERROR_CODES_H

#define PY_BUILDVALUE_NULL     1
#define PY_CALLOBJECT_NULL     2
#define PY_LISTCHECK_ERR       3
#define PY_EMPTY_LIST          4
#define PY_LIST_OVERFLOW       5
#define PY_NON_NUMERIC         6
#define PY_MACROS_FAIL         7
#define PY_CONVERT_FAIL        8
#define PY_DICT_SET_FAIL       9

#define ML_UNRECOG_MODEL      10
#define ML_UNRECOG_PARAM      11

#define GPT_NOT_CLASS          12
#define GPT_LIST_ATTR_NOT_FND  13
#define GPT_BAD_LIST           14
#define GPT_BAD_N_LIST_ELEM    15
#define GPT_MALLOC_FAIL        16
#define GPT_NOT_TUPLE          17
#define GPT_BAD_TUPLE_ITEM     18
#define GPT_ITEM_NOT_STR       19
#define GPT_FAIL_STR_COPY      20
#define GPT_ITEM_NOT_INT       21

#define PY_STRING_READ_FAIL    22
#define PY_IMPORT_FAIL         23

#define CA_ATTR_READ_FAIL            24
#define CA_LIST_ITEM_READ_FAIL       25
#define CA_LISTLIST_ITEM_READ_FAIL   26
#define CA_TYPES_DONT_MATCH          27
#define CA_ILLEGAL_PAR               28

#define RPI_ATTR_READ_FAIL           29
#define RPI_LIST_ITEM_READ_FAIL      30
#define RPI_MALLOC_FAIL              31


#endif /* ERROR_CODES_H */


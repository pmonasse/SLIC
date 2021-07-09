/*
 * SPDX-License-Identifier: GPL-3.0-or-later or BSD-2-Clause
 * @file io_jpg.h
 * @brief JPEG input/output
 *
 * Copyright (c) 2021, Pascal Monasse <pascal.monasse@enpc.fr>
 */

#ifndef IO_JPG_H
#define IO_JPG_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>

unsigned char *io_jpg_read_u8(const char* fname, size_t*nx,size_t*ny,size_t*nc);
unsigned char *io_jpg_read_u8_rgb(const char *fname, size_t *nxp, size_t *nyp);
unsigned char *io_jpg_read_u8_gray(const char *fname, size_t *nxp, size_t *nyp);
int io_jpg_write_u8(const char *fname, const unsigned char *data,
                    size_t nx, size_t ny, size_t nc);

#ifdef __cplusplus
}
#endif

#endif


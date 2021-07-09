/*
 * SPDX-License-Identifier: GPL-3.0-or-later or BSD-2-Clause
 * @file io_jpg.c
 * @brief JPEG input/output
 *
 * Copyright (c) 2021, Pascal Monasse <pascal.monasse@enpc.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include "io_jpg.h"
#include <stdlib.h>
#include <stdio.h>
#include <jpeglib.h>
#include <setjmp.h>
#include <string.h>

struct my_error_mgr {
  struct jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

METHODDEF(void) jpeg_error(j_common_ptr cinfo) {
  struct my_error_mgr *myerr = (struct my_error_mgr*) cinfo->err;
  (*cinfo->err->output_message)(cinfo);
  longjmp(myerr->setjmp_buffer, 1);
}

unsigned char* io_jpg_read_u8(const char *fname,
                              size_t* nxp, size_t* nyp, size_t* ncp) {
    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;
    unsigned char* data = NULL;
    FILE* fp = NULL;
    int stride;

    if (NULL == fname || NULL == nxp || NULL == nyp || NULL == ncp)
        return NULL;
    
    /* Open the JPG input file */
    if(strcmp(fname,"-")==0)
        fp = stdin;
    else if((fp=fopen(fname,"rb"))==NULL)
        return NULL;

    /* Error handling */
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = &jpeg_error;
    if(setjmp(jerr.setjmp_buffer)) {
        jpeg_destroy_decompress(&cinfo);
        free(data);
        if(stdin != fp)
            fclose(fp);
        return NULL;
    }

    /* Read header */
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, fp);
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);
    *nxp = (size_t)cinfo.output_width;
    *nyp = (size_t)cinfo.output_height;
    *ncp = (size_t)cinfo.output_components;

    /* Read contents */
    stride = cinfo.output_width * cinfo.output_components;
    data = malloc(*nxp * *nyp * *ncp);
    while(cinfo.output_scanline < cinfo.output_height) {
        JSAMPLE* row = data+cinfo.output_scanline*stride;
        jpeg_read_scanlines(&cinfo, &row, 1);
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    if(stdin != fp)
        fclose(fp);

    return data;
}

unsigned char *io_jpg_read_u8_rgb(const char *fname, size_t *nxp, size_t *nyp) {
    unsigned char *data, *datargb=NULL, *p, *prgb;
    size_t nc=3, i;
    data = io_jpg_read_u8(fname, nxp, nyp, &nc);
    if(data==NULL || nc==3)
        return data;
    if(nc==1) {
        prgb = datargb = malloc(*nxp * *nyp * 3);
        p=data;
        for(i=0; i<*nxp**nyp; i++, p++) {
            *prgb++ = *p;
            *prgb++ = *p;
            *prgb++ = *p;
        }
    }
    free(data);
    return datargb;
}

unsigned char rgb_gray(unsigned char r, unsigned char g, unsigned char b) {
    return (unsigned char) (6969 * r + 23434 * g + 2365 * b) / 32768;
}


unsigned char *io_jpg_read_u8_gray(const char *fname, size_t *nxp, size_t *nyp){
    unsigned char *data, *datagray=NULL, *prgb, *p;
    size_t nc=1, i;
    data = io_jpg_read_u8(fname, nxp, nyp, &nc);
    if(data==NULL || nc==1)
        return data;
    else {
        p = datagray = malloc(*nxp * *nyp);
        prgb=data;
        for(i=0; i<*nxp**nyp; i++, prgb+=3)
            *p++ = rgb_gray(prgb[0],prgb[1],prgb[2]);
    }
    free(data);
    return datagray;
}

int io_jpg_write_u8(const char *fname, const unsigned char *data,
                    size_t nx, size_t ny, size_t nc) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE* fp = NULL;
    int stride;
    JSAMPLE* row;

    if(nx<=0 || ny<=0 || (nc!=1 && nc!=3) || fname==NULL || data==NULL)
        return -1;
    
    /* Open the JPG output file */
    if(strcmp(fname,"-")==0)
        fp = stdout;
    else if((fp=fopen(fname,"wb"))==NULL)
        return -1;

    /* Prepare compressor */
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, fp);
    cinfo.image_width = (JDIMENSION)nx;
    cinfo.image_height = (JDIMENSION)ny;
    cinfo.input_components = (JDIMENSION)nc;
    cinfo.in_color_space = (nc==3)? JCS_RGB: JCS_GRAYSCALE;
    jpeg_set_defaults(&cinfo);

    /* Write contents */
    jpeg_start_compress(&cinfo, TRUE);
    stride = nx * nc;
    while(cinfo.next_scanline < cinfo.image_height) {
        row = (JSAMPLE*)(data + stride*cinfo.next_scanline);
        jpeg_write_scanlines(&cinfo, &row, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);
    if(stdout != fp)
        fclose(fp);
    return 0;
}

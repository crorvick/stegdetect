/*
 * Copyright 2001 Niels Provos <provos@citi.umich.edu>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *      This product includes software developed by Niels Provos.
 * 4. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <err.h>
#include <string.h>
#include <setjmp.h>

#include <jpeglib.h>

#include "config.h"
#include "jphide_table.h"
#include "common.h"

struct njvirt_barray_control {
  JBLOCKARRAY mem_buffer;       /* => the in-memory buffer */
  JDIMENSION rows_in_array;     /* total virtual array height */
  JDIMENSION blocksperrow;      /* width of array (and of memory buffer) */
  JDIMENSION maxaccess;         /* max rows accessed by access_virt_barray */
  JDIMENSION rows_in_mem;       /* height of memory buffer */
  JDIMENSION rowsperchunk;      /* allocation chunk size in mem_buffer */
  JDIMENSION cur_start_row;     /* first logical row # in the buffer */
  JDIMENSION first_undef_row;   /* row # of first uninitialized row */
  boolean pre_zero;             /* pre-zero mode requested? */
  boolean dirty;                /* do current buffer contents need written? */
  boolean b_s_open;             /* is backing-store data valid? */
  jvirt_barray_ptr next;        /* link to next virtual barray control block */
  void *b_s_info;  /* System-dependent control info */
};

#define MAX_COMMENTS	10
u_char *comments[MAX_COMMENTS+1];
size_t commentsize[MAX_COMMENTS+1];
int ncomments;

typedef struct njvirt_barray_control *njvirt_barray_ptr;
njvirt_barray_ptr *dctcoeff;

JBLOCKARRAY dctcompbuf[MAX_COMPS_IN_SCAN];
int hib[MAX_COMPS_IN_SCAN], wib[MAX_COMPS_IN_SCAN];
struct jpeg_decompress_struct jinfo;

/* Comment processing */
u_char
jpeg_getc(j_decompress_ptr cinfo)
{
	struct jpeg_source_mgr *datasrc = cinfo->src;

	if (datasrc->bytes_in_buffer == 0) {
		if (! (*datasrc->fill_input_buffer) (cinfo))
			err(1, __FUNCTION__": fill_input");
	}
	datasrc->bytes_in_buffer--;
	return GETJOCTET(*datasrc->next_input_byte++);
}

METHODDEF(boolean)
comment_handler(j_decompress_ptr cinfo)
{
	u_int32_t length;
	u_char *p;
	
	length = jpeg_getc(cinfo) << 8;
	length += jpeg_getc(cinfo);
	length -= 2;

	p = malloc(length);
	if (p == NULL)
		return (FALSE);

	commentsize[ncomments] = length;
	comments[ncomments++] = p;

	while (length-- > 0) {
		*p++ = jpeg_getc(cinfo);
	}

	return (TRUE);
	
}

void
comments_init(void)
{
	memset(comments, 0, sizeof(comments));
	memset(commentsize, 0, sizeof(commentsize));
	ncomments = 0;
}

void
comments_free(void)
{
	int i;

	for (i = 0; i < ncomments; i++)
		free(comments[i]);
	ncomments = 0;
}

void
stego_set_callback(void (*cb)(int, short), enum order order)
{
	extern void (*stego_mcu_order)(int, short);
	extern void (*stego_natural_order)(int, short);

	switch (order) {
	case ORDER_MCU:
		stego_mcu_order = cb;
		break;
	case ORDER_NATURAL:
		stego_natural_order = cb;
		break;
	}
}

char *
fgetl(char *s, int size, FILE *stream)
{
        char *res, *pos;
        int c;

        if ((res = fgets(s, size, stream))) {
                if (!*res) return res;

                pos = res + strlen(res) - 1;
                if (*pos == '\n') {
                        *pos = 0;
                        if (pos > res)
				if (*--pos == '\r') *pos = 0;
                } else
			if ((c = getc(stream)) == '\n') {
				if (*pos == '\r') *pos = 0;
			} else
				while (c != EOF && c != '\n')
					c = getc(stream);
        }

        return (res);
}

short **pjdcts, *cbdcts;
int *pjbits, ncbbits;

void
jsteg_cb(int where, short val)
{
	int count;
	
	if ((val & 0x01) == val)
		return;

	count = *pjbits;
	
	if (count >= ncbbits) {
		if (ncbbits != 0)
			ncbbits *= 2;
		else
			ncbbits = 256;
		cbdcts = realloc(cbdcts, ncbbits * sizeof(short));
		if (cbdcts == NULL)
			err(1, "realloc");
		*pjdcts = cbdcts;
	}

	cbdcts[count] = val;

	(*pjbits)++;
}

int
prepare_jsteg(short **pdcts, int *pbits)
{
	pjdcts = pdcts;
	pjbits = pbits;

	stego_set_callback(jsteg_cb, ORDER_MCU);
	*pdcts = cbdcts = NULL;
	ncbbits = *pjbits = 0;

	return (0);
}

short **podcts, *cbodcts;
int *pobits, ncbobits;

void
outguess_cb(int where, short val)
{
	int count;
	
	if ((val & 0x01) == val || where == 0)
		return;

	count = *pobits;
	
	if (count >= ncbobits) {
		if (ncbobits != 0)
			ncbobits *= 2;
		else
			ncbobits = 256;
		cbodcts = realloc(cbodcts, ncbobits * sizeof(short));
		if (cbodcts == NULL)
			err(1, "realloc");
		*podcts = cbodcts;
	}

	cbodcts[count] = val;

	(*pobits)++;
}

int
prepare_outguess(short **pdcts, int *pbits)
{
	podcts = pdcts;
	pobits = pbits;

	stego_set_callback(outguess_cb, ORDER_NATURAL);
	*pdcts = cbodcts = NULL;
	ncbobits = *pobits = 0;

	return (0);
}

int
prepare_all(short **pdcts, int *pbits)
{
	int comp, row, col, val, bits, i;
	short *dcts;

	bits = 0;
	for (comp = 0; comp < 3; comp++)
		bits += hib[comp] * wib[comp] * DCTSIZE2;

	dcts = malloc(bits * sizeof (short));
	if (dcts == NULL) {
		warn(__FUNCTION__": malloc");
		return (-1);
	}

	bits = 0;
	for (comp = 0; comp < 3; comp++) 
		for (row = 0 ; row < hib[comp]; row++)
			for (col = 0; col < wib[comp]; col++)
				for (i = 0; i < DCTSIZE2; i++) {
					val = dctcompbuf[comp][row][col][i];
					
					dcts[bits++] = val;
				}

	*pdcts = dcts;
	*pbits = bits;

	return (0);
}

int
jsteg_size(short *dcts, int bits, int *off)
{
	int count, width, len, i, mode;
	short val;

	width = mode = len = count = 0;
	for (i = 0; i < bits; i++) {
		val = dcts[i];
		
		if (mode == 0)
			width = (width << 1) | (val & 0x1);
		else
			len = (len << 1) | (val & 0x1);

		count++;

		if (mode == 0 && count >= 5) {
			mode = 1;
			count = 0;
		} else if (mode == 1 && count >= width)
			goto out;
	}

 out:
	/* Save where we left off reading the length */
	if (off != NULL)
		*off = (i + 1);
	
	return (len);
}

int
prepare_normal(short **pdcts, int *pbits)
{
	int comp, row, col, val, bits, i;
	short *dcts = NULL;

	bits = 0;
	for (comp = 0; comp < 3; comp++)
		bits += hib[comp] * wib[comp] * DCTSIZE2;

	if (pdcts != NULL) {
		dcts = malloc(bits * sizeof (short));
		if (dcts == NULL) {
			warn(__FUNCTION__": malloc");
			return (-1);
		}
	}
	
	bits = 0;
	for (comp = 0; comp < 3; comp++) 
		for (row = 0 ; row < hib[comp]; row++)
			for (col = 0; col < wib[comp]; col++)
				for (i = 0; i < DCTSIZE2; i++) {
					val = dctcompbuf[comp][row][col][i];
					
					/* Skip 0 and 1 coeffs */
					if ((val & 1) == val)
						continue;

					if (dcts != NULL)
						dcts[bits] = val;
					bits++;
				}

	if (pdcts != NULL)
		*pdcts = dcts;
	*pbits = bits;

	return (0);
}

int
prepare_jphide(short **pdcts, int *pbits)
{
	int comp, val, bits, i, mbits;
	int lwib[MAX_COMPS_IN_SCAN];
	int spos, nheight, nwidth, j, off;
	short *dcts = NULL;
	char *back[3];

	for (i = 0; i < 3; i++)
		lwib[i] = 64 * wib[i] - 1;

	mbits = 0;
	for (comp = 0; comp < 3; comp++) {
		int off = hib[comp] * wib[comp] * DCTSIZE2;
		mbits += off;

		/* XXX - wasteful */
		back[comp] = calloc(off, sizeof (char));
		if (back[comp] == NULL) {
			warn(__FUNCTION__": calloc");
			goto err;
		}
	}

	if (pdcts != NULL) {
		dcts = malloc(mbits * sizeof (short));
		if (dcts == NULL) {
			warn(__FUNCTION__": malloc");
			goto err;
		}
	}

	comp = ltab[0];
	spos = ltab[1];
	nheight = 0;
	nwidth = spos - 64;
	bits = j = 0;

	while (bits < mbits) {
		nwidth += DCTSIZE2;
		if (nwidth > lwib[comp]) {
			nwidth = spos;
			nheight++;
			if (nheight >= hib[comp]) {
				j += 3;
				if (ltab[j] < 0)
					goto out;

				comp = ltab[j];
				nwidth = spos = ltab[j + 1];
				nheight = 0;
			}
		}
		/* Protect IV */
		if (comp == 0 && nheight == 0 && nwidth <= 7)
			continue;

		val = dctcompbuf[comp][nheight][nwidth / DCTSIZE2][nwidth % DCTSIZE2];

		/* XXX - Overwrite so that we remember where we are */
		off = nheight * wib[comp] * DCTSIZE2 + nwidth;
		if (back[comp][off])
			break;
		back[comp][off] = 1;

		if (dcts != NULL)
			dcts[bits] = val;
		bits++;
	}
 out:
	for (i = 0; i < 3; i++)
		free(back[i]);

	if (pdcts != NULL)
		*pdcts = dcts;
	*pbits = bits;

	return (0);
 err:
	for (i = 0; i < 3; i++)
		if (back[i] != NULL)
			free(back[i]);
	return (-1);
}

struct my_error_mgr {
	struct jpeg_error_mgr pub;    /* "public" fields */

	jmp_buf setjmp_buffer;        /* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
	my_error_ptr myerr = (my_error_ptr) cinfo->err;

	/* Return control to the setjmp point */
	longjmp(myerr->setjmp_buffer, 1);
}

void
jpg_finish(void)
{
	jpeg_finish_decompress(&jinfo);
	comments_free();
}

void
jpg_destroy(void)
{
	jpeg_destroy_decompress(&jinfo);
	comments_free();
}

int
jpg_open(char *filename)
{
	char outbuf[1024];
	int i;
	struct my_error_mgr jerr;
	jpeg_component_info *compptr;
	FILE *fin;

	comments_init();
	
	if ((fin = fopen(filename, "r")) == NULL) {
		int error = errno;

		fprintf(stderr, "%s : error: %s\n",
			filename, strerror(error));
		return (-1);
	}

	jinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	/* Establish the setjmp return context for my_error_exit to use. */
	if (setjmp(jerr.setjmp_buffer)) {
		/* Always display the message. */
		(*jinfo.err->format_message) ((j_common_ptr)&jinfo, outbuf);

		fprintf(stderr, "%s : error: %s\n", filename, outbuf);

		/* If we get here, the JPEG code has signaled an error.
		 * We need to clean up the JPEG object, close the input file,
		 * and return.
		 */
		jpeg_destroy_decompress(&jinfo);

		fclose(fin);
		return (-1);
	}
	jpeg_create_decompress(&jinfo);
	jpeg_set_marker_processor(&jinfo, JPEG_COM, comment_handler);
	jpeg_stdio_src(&jinfo, fin);
	jpeg_read_header(&jinfo, TRUE);

	jinfo.quantize_colors = TRUE;
	dctcoeff = (njvirt_barray_ptr *)jpeg_read_coefficients(&jinfo);

	fclose(fin);

	if (jinfo.out_color_space != JCS_RGB) {
		fprintf(stderr, "%s : error: is not a RGB image\n", filename);
		goto out;
	}
 
	i = jinfo.num_components;
	if (i != 3) {
		fprintf(stderr,
			"%s : error: wrong number of color components: %d\n",
			filename, i);
		goto out;
	}
	
	for(i = 0; i < 3; i++) {
		compptr = jinfo.cur_comp_info[i];
		/*
		fprintf(stderr,
			"input_iMCU_row: %d, v_samp_factor: %d\n",
			jinfo.input_iMCU_row * compptr->v_samp_factor,
			(JDIMENSION) compptr->v_samp_factor);

		fprintf(stderr, "hib: %d, wib: %d\n",
			jinfo.comp_info[i].height_in_blocks,
			jinfo.comp_info[i].width_in_blocks);
		*/

		wib[i] = jinfo.comp_info[i].width_in_blocks;
		hib[i] = jinfo.comp_info[i].height_in_blocks;
		dctcompbuf[i] = dctcoeff[i]->mem_buffer;
	}

	return (0);
out:
	jpg_destroy();

	return (-1);
}

int
file_hasextension(char *name, char *ext)
{
	int nlen = strlen(name);
	int elen = strlen(ext);

	if (nlen < elen)
		return (0);

	return (strcasecmp(name + nlen - elen, ext) == 0);
}

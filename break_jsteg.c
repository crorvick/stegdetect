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
#include <fcntl.h>
#include <string.h>
#include <err.h>
#include <md5.h>
#include <errno.h>
#include <unistd.h>

#include <jpeglib.h>

#include "config.h"
#include "common.h"
#include "arc4.h"
#include "break_jsteg.h"

#ifndef MIN
#define		MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define JSTEGBYTES	8

struct jstegobj {
	int skip;
	u_int8_t coeff[JSTEGBYTES];
};

int break_jsteg(struct jstegobj *, struct arc4_stream *);

void *
break_jsteg_read(char *filename)
{
	struct jstegobj *jstegob;
	int fd;

	fd = open(filename, O_RDONLY, 0);
	if (fd == -1) {
		fprintf(stderr, "%s : error: %s\n",
			filename, strerror(errno));
		return (NULL);
	}

	jstegob = malloc(sizeof(struct jstegobj));
	if (jstegob == NULL)
		err(1, "malloc");

	if (read(fd, jstegob, sizeof(*jstegob)) != sizeof(*jstegob)) {
		close(fd);
		free(jstegob);
		return (NULL);
	}

	jstegob->skip = ntohl(jstegob->skip);

	close(fd);

	return (jstegob);
}

int
break_jsteg_write(char *filename, void *arg)
{
	struct jstegobj *jstegob = arg;
	int fd;

	fd = open(filename, O_WRONLY|O_CREAT|O_TRUNC, 0644);
	if (fd == -1)
		return (-1);

	jstegob->skip = htonl(jstegob->skip);

	if (write(fd, jstegob, sizeof(*jstegob)) != sizeof(*jstegob)) {
		close(fd);
		return (-1);
	}

	close(fd);

	return (0);
}

void
break_jsteg_destroy(void *obj)
{
	free(obj);
}

void *
break_jsteg_prepare(short *dcts, int bits)
{
	struct jstegobj *jstegob;
	int i, j, max, bytes, jsbits;
	u_char *p;
	short val;
	
	jstegob = malloc(sizeof(struct jstegobj));
	if (jstegob == NULL)
		err(1, "malloc");

	max = sizeof(jstegob->coeff) * 8;
	bytes = jsteg_size(dcts, bits, &i);
	jsbits = bytes * 8;
	jstegob->skip = bytes - sizeof(jstegob->coeff);

	if (jsbits < max || i + jsbits > bits) {
		warnx(__FUNCTION__": bad size in bits, %d", bits);
		return (NULL);
	}

	p = jstegob->coeff;
	i += jsbits - max;
	for (j = 0; j < max; j++) {
		val = dcts[i++];

		if (j != 0 && (j % 8) == 0)
			p++;
		
		*p = (*p << 1) | (val & 0x01);
	}
	
	return (jstegob);
}

int
crack_jsteg(char *filename, char *word, void *obj)
{
	static u_char oword[57];	/* version 3 marker */
	static int init;
	static struct arc4_stream as;
	
	struct arc4_stream tas;
	struct jstegobj *jstegob = obj;
	int changed = 0;

	if (strcmp(word, oword)) {
		strlcpy(oword, word, sizeof(oword));
		changed = 1;
		init = 0;
	}

	if (!init || changed) {
		arc4_fixedkey(&as, word, strlen(word));
		init = 1;
	}

	tas = as;
	if (break_jsteg(jstegob, &tas)) {
		fprintf(stdout, "%s : jsteg(%s)\n", filename, word);
		return (1);
	}

	return (0);
}

int
break_jsteg(struct jstegobj *js, struct arc4_stream *as)
{
	u_char plain[JSTEGBYTES];
	u_char *p;
	int i;

	arc4_skipbytes(as, js->skip);

	p = js->coeff;
	for (i = 0; i < sizeof(js->coeff); i++)
		plain[i] = p[i] ^ arc4_getbyte(as);

	p = plain + sizeof(plain);
	
	if (memcmp(p - 7, "korejwa", 7) == 0)
		return (1);
	else if (memcmp(p - 4, "cMk", 3) == 0) {
		int n = *(p - 1);
		if (n == 4 || n == 5)
			return (1);
	}
	
	return (0);
}

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
#include <netinet/in.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <err.h>
#include <string.h>

#include <jpeglib.h>

#include "config.h"
#include "common.h"

#define VERSION "0.2"

#define DBG_PRINTHIST	0x0001
#define DBG_CHIDIFF	0x0002
#define DBG_CHICALC	0x0004
#define DBG_CHIEND	0x0008
#define DBG_PRINTONES	0x0010
#define DBG_CHI		0x0020
#define DBG_ENDVAL	0x0040
#define DBG_BINSRCH	0x0080

#define FLAG_DOOUTGUESS	0x0001
#define FLAG_DOJPHIDE	0x0002
#define FLAG_DOJSTEG	0x0004
#define FLAG_DOINVIS	0x0008

float chi2cdf(float chi, int dgf);

char *progname;

float DCThist[257];
float scale = 1;		/* Sensitivity scaling */

int debug = 0;

void
buildDCThist(short *data, int x, int y)
{
	int i, min, max;
	int off, count, sum;
	static short *olddata;
	static int oldx, oldy;

	if (olddata != data || x < oldx || y < oldy ||
	    x - oldx + y - oldy >= y - x) {
		olddata = data;
		oldx = x;
		oldy = y;

		memset(DCThist, 0, sizeof(DCThist));
	} else {
		for (i = oldx; i < x; i++) {
			off = data[i];

			/* Don't know what to do about DC! */
			if (off < -128)
				continue;
			else if (off > 127)
				continue;

			DCThist[off + 128]--;
		}

		olddata = data;
		oldx = x;

		x = oldy;

		oldy = y;
	}

	min = 2048;
	max = -2048;

	/* Calculate coefficent frequencies */
	sum = count = 0;
	for (i = x; i < y; i++) {
		if ((i & ~63) == i) {
			if (debug & DBG_PRINTONES)
				fprintf(stdout, "%d] %d\n", i, count);
			sum += count;
			count = 0;
		}

		off = data[i];
		if (off == 1)
			count++;

		if (off < min)
			min = off;
		if (off > max)
			max = off;

		/* Don't know what to do about DC! */
		if (off < -128)
			continue;
		else if (off > 127)
			continue;
		
		DCThist[off + 128]++;
	}

	if (debug & DBG_PRINTHIST) {
		for (i = 0; i < 256; i++) {
			fprintf(stdout, "%4d: %8.1f\n", i - 128, DCThist[i]);
		}

		fprintf(stdout, "Min: %d, Max: %d, Sum-1: %d\n",
			min, max, sum);
	}
}

/*
 * Self calibration on bad test example.
 */

int
unify_false_jsteg(float *hist, float *theo, float *obs, float *discard)
{
	int i, size = 0;

	/* Build theoretical histogram */
	for (i = 0; i < 128; i++) {
		if (i == 64 || i == 65 || i == 0)
			continue;

		theo[size] = (float)(hist[2*i - 1] + hist[2*i])/2;
		obs[size++] = hist[2*i];
	}

	return (size);
}

int
unify_false_outguess(float *hist, float *theo, float *obs, float *discard)
{
	int i, size = 0;
	int one, two;

	/* Build theoretical histogram */
	for (i = 0; i < 128; i++) {
		if (i == 64 || i == 65 || i == 0)
			continue;

		one = hist[2*i - 1];
		two = hist[2*i];

		theo[size] = (float)(one + two)/2;
		obs[size++] = two;
	}

	return (size);
}


int
unify_false_jphide(float *hist, float *theo, float *obs, float *discard)
{
	int i, size = 0;

	/* Build theoretical histogram */
	for (i = 0; i < 128; i++) {
		if (i == 63 || i == 64 || i == 65)
			continue;

		if (i < 64) {
			theo[size] = (float)(hist[2*i] + hist[2*i + 1])/2;
			obs[size++] = hist[2*i + 1];
		} else {
			theo[size] = (float)(hist[2*i - 1] + hist[2*i])/2;
			obs[size++] = hist[2*i];
		}
	}

	return (size);
}

int
unify_normal(float *hist, float *theo, float *obs, float *discard)
{
	int i, size = 0;

	/* Build theoretical histogram */
	for (i = 0; i < 128; i++) {
		if (i == 64)
			continue;

		theo[size] = (float)(hist[2*i] + hist[2*i + 1])/2;
		obs[size++] = hist[2*i + 1];
	}

	return (size);
}

int
unify_outguess(float *hist, float *theo, float *obs, float *pdiscard)
{
	int i, size = 0;
	int one, two, sum, discard;
	float f, fbar;

	sum = 0;
	for (i = 0; i < 256; i++) {
		if (i == 64 || i == 65)
			continue;
		sum += hist[i];
	}

	discard = 0;
	/* Build theoretical histogram */
	for (i = 0; i < 128; i++) {
		if (i == 64)
			continue;

		one = hist[2*i];
		two = hist[2*i + 1];

		if (one > two) {
			f = one;
			fbar= two;
		} else {
			f = two;
			fbar = one;
		}

		/* Try to check if outguess could have been used here
		 * If the smaller coefficient is less than a quarter
		 * of the larger one, then outguess has probably been
		 * not used, and we included the coefficient.
		 * Otherwise check if outguess modifications could
		 * have reduced the difference significantly.
		 */
		if ((fbar > f/4) &&
		    ((f - f/3) - (fbar + f/3) > 10)) {
			if ((debug & DBG_CHIDIFF) && (one || two))
				fprintf(stdout,
					"%4d: %8.3f - %8.3f skipped (%f)\n",
					i*2 - 128,
					(float)two,
					(float)(one + two)/2,
					(float)(one + two)/sum);
				
			discard += one + two;
			continue;
		}

		theo[size] = (float)(one + two)/2;
		obs[size++] = two;
	}

	*pdiscard = (float)discard/sum;
	return (size);
}


int
unify_jphide(float *hist, float *theo, float *obs, float *discard)
{
	int i, size;

	size = 0;
	/* First pass */
	for (i = 0; i < 256; i++) {
		/* Exclude special cases */
		if (i >= (-1 + 128) && i <= (1 + 128))
			continue;

		/* Lower bit = 0 */
		if (i < 128 && !(i & 1))
			continue;
		else if (i & 1)
			continue;

		theo[size] = (hist[i] + hist[i + 1])/2;
		obs[size++] = hist[i];
	}

	/* Special case for 1 and -1 *
	   theo[size] = (hist[-1 + 128] + hist[1 + 128] + hist[0 + 128])/2;
	   obs[size++] = hist[-1 + 128] + hist[1 + 128];
	*/
	return (size);
}


float
chi2test(short *data, int bits,
	 int (*unify)(float *, float *, float *, float *),
	 int a, int b)
{
	float DCTtheo[128];
	float DCTobs[128];
	float chi, sumchi, ymt, ytt, f, discard;
	int i, dgf, size;

	if (a < 0)
		a = 0;
	if (b > bits)
		b = bits;

	if (a >= b)
		return (-1);

	buildDCThist(data, a, b);

	discard = 0;
	size = (*unify)(DCThist, DCTtheo, DCTobs, &discard);

	ymt = ytt = 0;
	sumchi = 0;
	dgf = 0;
	for (i = 0; i < size; i++) {
		ymt += DCTobs[i];
		ytt += DCTtheo[i];

		if (debug & DBG_CHIDIFF) {
			if (DCTobs[i] || DCTtheo[i])
				fprintf(stdout, "%4d: %8.3f - %8.3f\n", i,
					DCTobs[i],
					DCTtheo[i]);
		}

		if (ytt >= 5) {
			/* Calculate chi^2 */
			chi = ymt - ytt;


			if (debug & DBG_CHICALC) {
				fprintf(stdout,
					"     (%8.3f - %8.3f)^2 = %8.3f / %8.3f = %8.3f | %8.3f\n",
					ymt, ytt,
					chi*chi, ytt, chi*chi/ytt, sumchi);
			}


			chi = chi*chi;
			chi /= ytt;

			sumchi += chi;

			dgf++;
			ymt = ytt = 0;
		}
	}

	f = 1 - chi2cdf(sumchi, dgf - 1);

	if (debug & DBG_CHIEND) {
		fprintf(stdout,
			"Categories: %d, Chi: %f, Q: %f, dis: %f -> %f\n",
			dgf, sumchi, f, discard, f * (1 - discard));
	}

	return (f * (1 - discard));
}

#define BINSEARCHVAR \
	float _max, _min, _good; \
	int _iteration

#define BINSEARCH(imin, imax, imaxiter) \
	percent = (imax); \
	_good = (imax); \
	_iteration = 0; \
	_min = (imin); \
	_max = (imax); \
	while (_iteration < (imaxiter))

#define BINSEARCH_NEXT(thresh) \
	if (debug & DBG_BINSRCH) \
		fprintf(stdout, "sum: %f, percent: %f,  good: %f\n", \
			sum, percent, _good); \
	if (_iteration == 0) { \
		if (sum >= (thresh)) \
			break; \
		_good = percent; \
		percent = _min; \
	} else \
	if (sum < (thresh)) { \
		_good = percent; \
		if (_good == _min) /* XXX */\
			break; /* XXX */\
		_max = percent; \
		percent = (_max - _min)/2 + _min; \
	} else { \
		_min = percent; \
		percent = (_max - _min)/2 + _min; \
	} \
	_iteration++

#define BINSEARCH_IFABORT(thresh) \
	percent = _good; \
	if (_good >= (thresh))

int
histogram_chi_jsteg(short *data, int bits)
{
	int off, length, minlen, maxlen;
	float f, sum, percent, i, count, where;
	float max, aftercount, scale, fs, start;
	BINSEARCHVAR;

	if (bits == 0)
		goto abort;
	start = 100*400/bits;
	if (start >= 7)
		goto abort;

	if (start < 0.5)
		start = 0.5;

	BINSEARCH(start, 7, 6) {
		sum = 0;
		for (i = 1; i <= 100; i += percent) {
			off = i*bits/100;
			f = chi2test(data, bits, unify_false_jsteg, 0, off);
			if (f == 0)
				break;
			if (f > 0.4)
				sum += f * (i > 1 ? percent : 1);
			if ((debug & DBG_CHI) && f != 0)
				fprintf(stdout, "%04f[:] %8.5f%%\n",
					i, f * 100);
		}

		BINSEARCH_NEXT(1);
	}

	BINSEARCH_IFABORT(7) {
	abort:
		if (debug & DBG_ENDVAL)
			fprintf(stdout,
				"Accumulation: no detection possible\n");
		return (-1);
	}

	where = count = 0;
	aftercount = max = 0;
	scale = 0.95;
	sum = 0;
	for (i = 1; i <= 100; i += percent) {
		off = i*bits/100;
		f = chi2test(data, bits, unify_normal, 0, off);
		if (f == 0)
			break;
		if (f > 0.4) {
			sum += f;
			count++;
		}
		if (f >= (max * scale)) {
			if (f > max) {
				max = f;
				/* More latitude for high values */
				fs = (max - 0.4) / 0.6;
				if (fs > 0)
					scale = 1 - (0.15*fs + (1 - fs)*0.05);
			}
			aftercount = -3;
			where = i;
		} else if (f > 0.05 * max) {
			if (aftercount >= 0)
				aftercount += f;
			else {
				aftercount++;
				where = i;
			}
		}

		if ((debug & DBG_CHI) && f != 0)
			fprintf(stdout, "%04f: %8.5f%%\n",
				i, f * 100);
	}

	length = jsteg_size(data, bits, NULL);
	minlen = where*bits/100/8;
	maxlen = (where+1)*bits/100/8;
	if (debug & DBG_ENDVAL) {
		fprintf(stdout,
		    "Accumulation (%4.1f%%): %f%% - %f (%f) (%d:%d - %d)\n",
		    percent, sum * 100, aftercount*100, count,
		    length, minlen, maxlen);
	}

	if (aftercount > 0)
		sum -= aftercount;
	/* Require a positive sum and at least two working samples */
	if (sum < 0 || count < 3)
		sum = 0;
	if (length < minlen/2 || length > maxlen*2)
		sum = 0;

	return (scale * sum / (2 * percent));
}

float norm_outguess[21] = {
	0.5,
	0.7071067811865475244,
	1,
	1.2247448713915890491,
	1.4142135623730950488,
	1.581138830084189666,
	1.73205080756887729353,
	1.87082869338697069279,
	2,
	2.1213203435596425732,
	2.23606797749978969641,
	2.34520787991171477728,
	2.4494897427831780982,
	2.54950975679639241501,
	2.6457513110645905905,
	2.73861278752583056728,
	2.8284271247461900976,
	2.91547594742265023544,
	3,
	3.08220700148448822513,
	3.162277660168379332
};

int
histogram_chi_outguess(short *data, int bits)
{
	int i, off, range;
	float percent, count;
	float f, sum, norm;
	BINSEARCHVAR;

	BINSEARCH(1, 10, 7) {
		range = percent*bits/100;
		sum = 0;
		for (i = 0; i <= 100; i ++) {
			off = i*bits/100;
			f = chi2test(data, bits, unify_false_outguess,
				     off - range, off + range);
			sum += f;
			if ((debug & DBG_CHI) && f != 0)
				fprintf(stdout, "%04d[:] %8.5f%%\n",
					i, f * 100);
		}

		BINSEARCH_NEXT(0.3);
	}

	/* XXX */
	BINSEARCH_IFABORT(10)
		return (0);

	range = percent*bits/100;
	count = 0;
	sum = 0;
	for (i = 0; i <= 100; i ++) {
		off = i*bits/100;
		f = chi2test(data, bits, unify_outguess,
			     off - range, off + range);
		if (f > 0.25)
			sum += f;
		if (f > 0.001)
			count++;
		if ((debug & DBG_CHI) && f != 0)
			fprintf(stdout, "%04d: %8.5f%%\n", i, f * 100);
	}

	count /= percent;

	off = percent * 2;
	if (off >= sizeof(norm_outguess)/sizeof(float))
		off = sizeof(norm_outguess)/sizeof(float) - 1;

	norm = sum / norm_outguess[off];

	if (debug & DBG_ENDVAL)
		fprintf(stdout,
			"Accumulation (%4.1f%%): %8.3f%% (%8.3f%%) (%4.1f)\n",
			percent,
			sum * 100,
			norm * 100,
			count);

	/* XXX - some wild adjustment */
	if (count < 15) {
		sum -= (15 - count) * 0.5;
		if (sum < 0)
			sum = 0;
	}

	return (scale * sum / 0.5);
}

int
histogram_chi_jphide(short *data, int bits)
{
	int i, off, range;
	float f, sum, percent;
	BINSEARCHVAR;

	BINSEARCH(1, 10, 7) {
		range = percent*bits/100;
		sum = 0;
		for (i = 1; i <= 100; i ++) {
			off = i*bits/100;
			f = chi2test(data, bits, unify_false_jphide,
				     0, off);
			sum += f;
			if ((debug & DBG_CHI) && f != 0)
				fprintf(stdout, "%04d[:] %8.5f%%\n",
					i, f * 100);
		}

		BINSEARCH_NEXT(3);
	}

	BINSEARCH_IFABORT(10)
		return (0);

	range = percent * bits/100;
	sum = 0;
	for (i = 1; i <= 100; i ++) {
		off = i*bits/100;
		f = chi2test(data, bits, unify_jphide,
			     0, off);
		if (f > 0.3)
			sum += f;
		else if (f < 0.2)
			break;

		if ((debug & DBG_CHI) && f != 0)
			fprintf(stdout, "%04d: %8.5f%%\n", i, f * 100);
	}

	if (debug & DBG_ENDVAL)
		fprintf(stdout, "Accumulation (%4.1f%%): %f%%\n",
			percent,
			sum * 100);

	return (scale * sum / 7);
}

#define RANGE	500
#define ADAPT	0.2

void
histogram_test(short *data, int bits)
{
	int i;
	float ratio;
	int runsum, which = -1, currun;
	int lastsum;
	int one;

	one = 0;
	ratio = 0.5;
	runsum = 0;
	lastsum = 0;
	currun = 0;

	for (i = 0; i < bits; i++) {
		if (which != (abs(data[i]) & 1)) {
			which = abs(data[i]) & 1;
			runsum += currun * currun;
			currun = 0;
		} else
			currun ++;

		if (abs(data[i]) & 1)
			one++;
		if (i && (i % RANGE == 0)) {
			ratio = (1-ADAPT)*ratio + ADAPT*(float)one/RANGE;
			one = 0;

			runsum += currun * currun;
			if (lastsum == 0)
				lastsum = runsum;
			else
				lastsum = (1-ADAPT)*lastsum + ADAPT*runsum;
			
			fprintf(stdout, "%7.3f: %8.3f %5.3f\n",
				(float)i/bits*100,
				(float)lastsum / RANGE,
				ratio);

			which = -1;
			currun = 0;
			runsum = 0;
		}
	}
}

#define RLE_MAX 1024

void
histogram_rle(short *data, int bits)
{
	int run, which, i;
	int hist[3][RLE_MAX], score[3], sum[3];

	which = 3;
	run = 0;

	memset(hist, 0, sizeof(hist));

	for (i = 0; i < bits; i++) {
		if (data[i] < -1 || data[i] > 1)
			continue;
		if (data[i] != which) {
			if (run > RLE_MAX)
				run = RLE_MAX - 1;
			if (which != 3) {
				hist[which + 1][run]++;
			}
			which = data[i];
			run = 0;
		}

		run++;
	}

	score[0] = sum[0] = 0;
	score[1] = sum[1] = 0;
	score[2] = sum[2] = 0;
	for (i = 0; i < RLE_MAX; i++) {
		score[0] += hist[0][i] * i * i;
		sum[0] += hist[0][i] * i;
		score[1] += hist[1][i] * i * i;
		sum[1] += hist[1][i] * i;
		score[2] += hist[2][i] * i * i;
		sum[2] += hist[2][i] * i;

		if (hist[0][i] + hist[1][i] + hist[2][i] == 0)
			continue;

		fprintf(stdout, "%4d: %5d %5d %5d\n",
			i, hist[0][i], hist[1][i], hist[2][i]);
	}

	fprintf(stderr, "-1: %10.3f 0: %10.3f 1: %10.3f\n",
		score[0]/(float)sum[0],
		score[1]/(float)sum[1],
		score[2]/(float)sum[2]
		);
}

void
usage(void)
{
	fprintf(stderr,
		"Usage: %s [-V] [-s <float>] [-d <num>] [-t <tests>] [file.jpg ...]\n",
		progname);
}

char *
quality(char *prepend, int q)
{
	static char quality[1024];
	char stars[4];
	int i;

	for (i = 0; i < q && i < 3; i++)
		stars[i] = '*';
	stars[i] = 0;

	snprintf(quality, sizeof(quality), "%s(%s)", prepend, stars);

	return (quality);
}

void
detect(char *filename, int scans)
{
	char outbuf[1024];
	int bits, jbits;
	int res, flag;
	short *dcts, *jdcts;

	if (scans & FLAG_DOJSTEG) {
		prepare_jsteg(&jdcts, &jbits);
	}
	
	if (jpg_open(filename) == -1) {
		if (jdcts != NULL)
			free(jdcts);
		return;
	}

	if (scans & FLAG_DOJSTEG) {
		stego_set_callback(NULL, ORDER_MCU);
	}
	
	flag = 0;
	sprintf(outbuf, "%s :", filename);

	if (scans & FLAG_DOINVIS) {
		extern u_char *comments[];
		extern size_t commentsize[];
		extern int ncomments;
		u_char *p;
		u_int32_t length;

		if (ncomments < 2)
			goto no_invisiblesecrets;

		/* Check first file */
		p = comments[1];
		length = p[3] << 24;
		length |= p[2] << 16;
		length |= p[1] << 8;
		length |= p[0];

		if (commentsize[1] == length + 4) {
			char tmp[128];

			flag = 1;
			snprintf(tmp, sizeof(tmp), " invisible(%d)", length);
			strlcat(outbuf, tmp, sizeof(outbuf));
		}
		
	no_invisiblesecrets:
	}
	
	if (scans & FLAG_DOJSTEG) {
		/* Set via the callback */
		dcts = jdcts;
		bits = jbits;

		if (dcts == NULL)
			goto jsteg_error;
		
		res = histogram_chi_jsteg(dcts, bits);
		if (res > 0) {
			strlcat(outbuf, quality(" jsteg", res),
				sizeof(outbuf));
			flag = 1;

			/* If this detects positivly so will outguess|jphide */
			scans &= ~(FLAG_DOOUTGUESS|FLAG_DOJPHIDE);
		}

		/* Special case to disable other methods for images, that
		 * will likelty to be false positive
		 */

		if (res == -1) {
			strlcat(outbuf, " skipped (false positive likely)",
				sizeof(outbuf));
			flag = 1;
			scans &= ~(FLAG_DOOUTGUESS|FLAG_DOJPHIDE);
		}

		free(dcts);
	jsteg_error:
	}

	if ((scans & FLAG_DOOUTGUESS) && prepare_normal(&dcts, &bits) != -1) {
		res = histogram_chi_outguess(dcts, bits);
		if (res) {
			strlcat(outbuf, quality(" outguess(old)", res),
				sizeof(outbuf));
			flag = 1;
		}
		free(dcts);
	}

	if ((scans & FLAG_DOJPHIDE) && prepare_jphide(&dcts, &bits) != -1) {
		res = histogram_chi_jphide(dcts, bits);
		if (res) {
			strlcat(outbuf, quality(" jphide", res),
				sizeof(outbuf));
			flag = 1;
		}
		free(dcts);
	}

	if (!flag)
		strlcat(outbuf, " negative", sizeof(outbuf));

	fprintf(stdout, "%s\n", outbuf);

	jpg_finish();
	jpg_destroy();
}

int
main(int argc, char *argv[])
{
	int i, scans;
	extern char *optarg;
	extern int optind;
	char ch;

	progname = argv[0];

	scans = FLAG_DOOUTGUESS | FLAG_DOJPHIDE | FLAG_DOJSTEG | FLAG_DOINVIS;

	/* read command line arguments */
	while ((ch = getopt(argc, argv, "s:Vd:t:")) != -1)
		switch((char)ch) {
		case 's':
			if ((scale = atof(optarg)) == 0) {
				usage();
				exit(1);
			}
			break;
		case 'V':
			fprintf(stdout, "Stegdetect Version %s\n", VERSION);
			exit(1);
		case 'd':
			debug = atoi(optarg);
			break;
		case 't':
			scans = 0;
			for (i = 0; i < strlen(optarg); i++)
				switch(optarg[i]) {
				case 'o':
					scans |= FLAG_DOOUTGUESS;
					break;
				case 'j':
					scans |= FLAG_DOJSTEG;
					break;
				case 'p':
					scans |= FLAG_DOJPHIDE;
					break;
				case 'i':
					scans |= FLAG_DOINVIS;
					break;
				default:
					usage();
					exit(1);
				}
			break;
		default:
			usage();
			exit(1);
		}

	argc -= optind;
	argv += optind;

	setvbuf(stdout, NULL, _IOLBF, 0);

	if (argc > 0) {
		while (argc) {
			detect(argv[0], scans);
			
			argc--;
			argv++;
		}
	} else {
		char line[1024];

		while (fgetl(line, sizeof(line), stdin) != NULL)
			detect(line, scans);
	}

	exit(0);
}

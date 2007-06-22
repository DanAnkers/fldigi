// ----------------------------------------------------------------------------
// rtty.h  --  RTTY modem
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code 
// distribution.
//  gmfsk Copyright (C) 2001, 2002, 2003
//  Tomi Manninen (oh2bns@sral.fi)
//
// fldigi is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ----------------------------------------------------------------------------

#ifndef _RTTY_H
#define _RTTY_H

#include "complex.h"
#include "modem.h"
#include "trx.h"
#include "misc.h"
#include "filters.h"
#include "fftfilt.h"
#include "digiscope.h"
#include "id.h"

#define	RTTY_SampleRate	8000
#define	RTTYMaxSymLen	1024

#define	LETTERS	0x100
#define	FIGURES	0x200

#define dispwidth 100

enum RTTY_RX_STATE {
	RTTY_RX_STATE_IDLE = 0,
	RTTY_RX_STATE_START,
	RTTY_RX_STATE_DATA,
	RTTY_RX_STATE_PARITY,
	RTTY_RX_STATE_STOP,
	RTTY_RX_STATE_STOP2
};

enum RTTY_PARITY {
	PARITY_NONE = 0,
	PARITY_EVEN,
	PARITY_ODD,
	PARITY_ZERO,
	PARITY_ONE
};

extern double _SHIFT[];
extern double _BAUD[];
extern int _BITS[];

//enum TTY_MODE { LETTERS, FIGURES };

class rtty : public modem {
private:
// waterfall ID
	id				*wfid;

	double shift;
	int symbollen;
	int nbits;
//	RTTY_PARITY parity;
	int stoplen;
	int msb;
	bool useFSK;

	double		phaseacc;
	double		rtty_squelch;
	double		rtty_shift;
	double		rtty_baud;
	int 		rtty_bits;
	RTTY_PARITY	rtty_parity;
	int			rtty_stop;
	bool 		rtty_reverse;
	bool		rtty_msbfirst;

	C_FIR_filter	*hilbert;
	C_FIR_filter	*lpfilt;
	Cmovavg *bitfilt;
	fftfilt *bpfilt;

	double pipe[RTTYMaxSymLen];
	int pipeptr;

	double bbfilter[RTTYMaxSymLen];
	unsigned int filterptr;

	RTTY_RX_STATE rxstate;

	int counter;
	int bitcntr;
	int rxdata;
	double posfreq, negfreq;
	double poserr, negerr;
	double avgratio;
	int poscnt, negcnt;

	double prevsymbol;
	complex prevsmpl;

	int rxmode;
	int txmode;
	int preamble;
	
	double sum, sum2, mean, stddev;
	
	void clear_syncscope();
	void update_syncscope();
	inline complex mixer(complex in);
	unsigned char bitreverse(unsigned char in, int n);
	int decode_char();
	int rttyparity(unsigned int);
	int rx(bool bit);
	double nco(double freq);
	void send_symbol(int symbol);
	void send_stop();
	void send_char(int c);
	void send_idle();
//	void keyline(int);
	int rttyxprocess();
	int baudot_enc(unsigned char data);
	char baudot_dec(unsigned char data);
	void Metric();
public:
	rtty(trx_mode mode);
	~rtty();
	void init();
	void rx_init();
	void tx_init(cSound *sc);
	void restart();
	int rx_process(double *buf, int len);
	int tx_process();

	void searchDown();
	void searchUp();

};

#endif
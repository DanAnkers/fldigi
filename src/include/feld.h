//
//	feld.h  --  FELDHELL modem
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code 
// distribution.
//  gmfsk Copyright (C) 2001, 2002, 2003
//  Tomi Manninen (oh2bns@sral.fi)
//  Copyright (C) 2004
//  Lawrence Glaister (ve7it@shaw.ca)
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


#ifndef _FELD_H
#define _FELD_H

#include "modem.h"
#include "id.h"
#include "filters.h"
#include "fftfilt.h"

#define	FeldSampleRate	8000
#define FeldMaxSymLen	1024
#define	FeldColumnRate	17.5
#define FeldBandWidth	245.0

#define	RxColumnLen	30
#define	TxColumnLen	14

#define	RxPixRate	((RxColumnLen)*(FeldColumnRate))
#define	TxPixRate	((TxColumnLen)*(FeldColumnRate))

#define	DownSampleInc	((double)(RxPixRate)/(samplerate))
#define	UpSampleInc	((double)(TxPixRate)/(samplerate))

#define	PIXMAP_W	14
#define	PIXMAP_H	(TxColumnLen)

class feld : public modem {
enum FELD_STATE {PREAMBLE, POSTAMBLE, DATA};
protected:
// waterfall ID
	id				*wfid;
//rx
	double rxphacc;
	double rxcounter;
	double agc;
	double peakhold;
	double minhold;

	C_FIR_filter	*hilbert;
	fftfilt			*bpfilt;
	Cmovavg			*bbfilt;
	Cmovavg			*minmaxfilt;
//	double bbfilter[MaxSymLen];
//	unsigned int filterptr;

//tx
	FELD_STATE	tx_state;
	double txphacc;
	double txcounter;
	double hell_bandwidth;
	
	int depth;
	int dxmode;
	int halfwidth;
	bool blackboard;
	bool hardkeying;

	int preamble;
	int postamble;
	int prevsymb;
	complex prev;
	
	double OnShape[80];
	double OffShape[80];
	
	int col_data[2*RxColumnLen];
	int col_pointer;
	int fntnbr;
	
	complex mixer(complex);
	double nco(double);
	void	rx(complex);
	void	FSKHELL_rx(complex);
	void	send_symbol(int currsymbol, int nextsymbol);
	void	send_null_column();
	void	tx_char(char);
	void	initKeyWaveform();
public:
	feld(trx_mode);
	~feld();
	void	init();
	void	rx_init();
	void	tx_init(cSound *sc);
	void 	restart() {};
	int		rx_process(double *buf, int len);
	int		tx_process();
	int		get_font_data(unsigned char c, int col);
};


#endif
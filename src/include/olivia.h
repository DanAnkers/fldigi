//
// Olivia.h
//
// Copyright (C) 2006
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.  Adapted from code contained in gmfsk source code 
// distribution.
//	Copyright (C) 2005
//	Tomi Manninen (oh2bns@sral.fi)
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
//
//

#ifndef _OLIVIA_H
#define _OLIVIA_H

#include "trx.h"
#include "id.h"
#include "jalocha/pj_mfsk.h"

class olivia : public modem {
private:
// waterfall ID
	id				*wfid;

	MFSK_Transmitter < float >*Tx;
	MFSK_Receiver < float >*Rx;
//	MFSK_Transmitter < double >*Tx;
//	MFSK_Receiver < double >*Rx;

	short int	*txbuffer;
	double		*txfbuffer;
	int 		txbufferlen;

	short int	*rxbuffer;
	int			rxbufferlen;

	double		phaseacc;
	complex		prevsymbol;
	int			preamble;
	unsigned int	shreg;

	int			escape;
	int			smargin;
	int			sinteg;
	bool		olivia_esc;
	int			tones;
	int			bw;

	int			preamblesent;
	int			postamblesent;
	double		preamblephase;

	double		txbasefreq;
	double		lastfreq;
	
	double		nco(double freq);
	void		send_preamble();
	void		send_postamble();
	
public:
	olivia();
	~olivia();
	void init();
	void rx_init();
	void tx_init(cSound *sc);
	void restart();
	int rx_process(double *buf, int len);
	int tx_process();
	int unescape(int c);
};


#endif
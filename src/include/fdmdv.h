// ----------------------------------------------------------------------------
// fdmdv.h  --  fdmdv modem
//
// Copyright (C) 2010
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.
//
// Fldigi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#ifndef _FDMDV_H
#define _FDMDV_H

#include "modem.h"
#include "openlpc.h"
#include "complex.h"
#include "psk.h"

#define	fdmdv_SampleRate	8000
#define FDMDV_FRAME_0 0
#define FDMDV_FRAME_1 1

#define CARRIER_BW 75
#define CARRIERS	15
#define BUFFER_FRAMES	1024 // Frames of date to buffer. Must be enough for rx_process not to overflow the buffer.

class fdmdv : public psk {
private:
	char varicoded_message[800]; // Each char stores 1 bit. Message can be 80 characters, varicode can be up to 8 bits plus 2 inter-character bits
	int message_length;
	int message_pointer;
	openlpc_encoder_state* voice_encoder;
	openlpc_decoder_state* voice_decoder;
	unsigned char bpsk_bit;
	unsigned char data_fifo[CARRIERS*BUFFER_FRAMES];
	int fifo_read_ptr;
	int fifo_write_ptr;

	void fdmdv_write_bpsk(int sym, int carrier, double* buffer, int len);
	void fdmdv_write_qpsk(int sym, int carrier, double* buffer, int len);
	void psk_init(void);
	void fifo_process(void);

// PSK stuff
        int                             symbollen;
        bool                    _qpsk;
        bool                    _pskr;
        double                  phaseacc[CARRIERS];
        complex                 prevsymbol[CARRIERS];
        unsigned int            shreg;
// rx variables & functions

        C_FIR_filter            *fir1;
        C_FIR_filter            *fir2;
//      C_FIR_filter            *fir3;
        double                  *fir1c;
        double                  *fir2c;
        Cmovavg                 *snfilt;
        Cmovavg                 *imdfilt;

        double                  I1[NUM_FILTERS];
        double                  I2[NUM_FILTERS];
        double                  Q1[NUM_FILTERS];
        double                  Q2[NUM_FILTERS];
        double                  COEF[NUM_FILTERS];
        double                  m_Energy[NUM_FILTERS];
        int                             m_NCount;
        bool                    imdValid;

        encoder                 *enc;
        viterbi                 *dec;
        unsigned int    bitshreg[CARRIERS];
        int                     rxbitstate;

        double                  phase[CARRIERS];
        double                  freqerr;
        int                             bits;
        double                  bitclk[CARRIERS];
        double                  syncbuf[16];
        double                  scope_pipe[2*PipeLen];//[PipeLen];
        unsigned int    pipeptr;
        unsigned int    dcdshreg;

        int                     dcd;
        int                             dcdbits;
        complex                 quality;
        int                             acquire;

        viewpsk*                pskviewer;
        pskeval*                evalpsk;

        unsigned char           rx_symbol(complex symbol, char type, int carrier);
        void                    rx_bit(int bit);
        void                    rx_bit2(int bit);
        void                    rx_qpsk(int bits);
        void                    rx_pskr(unsigned char symbol);
        double                  scopedata[16];
// IMD & s/n variables
        double                  k0, k1, k2;
        double                  I11, I12, I21, I22, I31, I32;
        double                  snratio, s2n, imdratio, imd;
        double                  E1, E2, E3;
        double                  afcmetric;
        
// tx variables & functions
        double                  *tx_shape;
        int                     preamble;
        void                    tx_symbol(int sym);
        void                    tx_bit(int bit);
        void                    tx_char(unsigned char c);
        void                    tx_flush();
        void                    update_syncscope();
        void                    signalquality();
        void                    findsignal();
        void                    phaseafc();
        void                    afc();
        void                    coreafc();

        void                    initSN_IMD();
        void                    resetSN_IMD();
        void                    calcSN_IMD(complex z);


protected:
	int varicode_encode(char* decoded, char* encoded);
	int varicode_decode(char* decoded, char* encoded);

public:
	fdmdv(trx_mode pskmode);
	~fdmdv();
	void init();
	void rx_init();
	void voicerx_init(SoundBase *vsc);
	void tx_init(SoundBase *sc);
	void voicetx_init();
	void restart();
	int rx_process(const double *buf, int len);
	int tx_process();
	char message[80];

};

#endif

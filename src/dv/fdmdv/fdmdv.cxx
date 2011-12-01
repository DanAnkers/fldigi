// ----------------------------------------------------------------------------
// fdmdv.cxx  --  fdmdv modem
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

#include <config.h>

#include "fdmdv.h"
#include "openlpc.h"
#include "modem.h"
#include "digiscope.h"
#include "fl_digi.h"
#include "configuration.h"

#include "debug.h"


#define FDMDV_VOICE_MS		40	//40ms voice samples are passed to the codec
#define FDMDV_BW					1125
#define FDMDV_DATA_BITS		2		//2 data bits per frame
#define FDMDV_CODEC_BITS	56	//56 bits come back from the 40ms voice
#define FDMDV_CODEC_BYTES ceil(FDMDV_CODEC_BITS/8)
#define FDMDV_SPACING			75	//75Hz spacing between carriers

void fdmdv::tx_init(SoundBase *sc)
{
	scard = sc;
	
	// message_length = varicode_encode(message, varicoded_message);
	message_length = 10;
	message_pointer = 0;
}

void fdmdv::voicetx_init()
{
}

void fdmdv::rx_init()
{
	put_MODEstatus(mode);
}

void fdmdv::voicerx_init(SoundBase *vsc)
{
	vscard = vsc;
}

void fdmdv::init()
{
	modem::init();
	rx_init();
	set_scope_mode(Digiscope::BLANK);
}

fdmdv::~fdmdv()
{
}

void fdmdv::restart()
{
	set_bandwidth(FDMDV_BW);
}

fdmdv::fdmdv(trx_mode pskmode) : psk(pskmode)
{
	mode = MODE_FDMDV;
	samplerate = 8000;
	vsamplerate = 8000;
	strncpy(message,"MD1CLV Testing",80);
	bpsk_bit = 0;

	// PSK stuff

        cap |= CAP_AFC | CAP_AFC_SR;

        // QPSK50
        symbollen = 160;
        dcdbits = 32; //?
        cap |= CAP_REV;

        // create impulse response for experimental FIR filters
        double fir1c[64];
        double fir2c[64];

        fir1 = new C_FIR_filter();
        fir2 = new C_FIR_filter();

        switch (progdefaults.PSK_filter) {
                case 1:
                // use the original gmfsk matched filters
                        for (int i = 0; i < 64; i++) {
                                fir1c[i] = gmfir1c[i];
                                fir2c[i] = gmfir2c[i];
                        }
                        fir1->init(FIRLEN, symbollen / 16, fir1c, fir1c);
                        fir2->init(FIRLEN, 1, fir2c, fir2c);
                        break;
                case 2:
                // creates fir1c matched sin(x)/x filter w hamming
                        wsincfilt(fir1c, 1.0 / symbollen, false);
                        fir1->init(FIRLEN, symbollen / 16, fir1c, fir1c);
                // creates fir2c matched sin(x)/x filter w hamming
                        wsincfilt(fir2c, 1.0 / 16.0, false);
                        fir2->init(FIRLEN, 1, fir2c, fir2c);
                        break;
                case 3:
                // creates fir1c matched sin(x)/x filter w hamming
                        wsincfilt(fir1c, 1.0 / symbollen, false);
                        fir1->init(FIRLEN, symbollen / 16, fir1c, fir1c);
                // 1/22 with Hamming window nearly identical to gmfir2c
                        wsincfilt(fir2c, 1.0 / 22.0, false);
                        fir2->init(FIRLEN, 1, fir2c, fir2c);
                        break;
                case 4:
                        fir1->init_lowpass (FIRLEN, 16, 1.5 / symbollen);
                        wsincfilt(fir2c, 1.5 / 16.0, true);
                        fir2->init(FIRLEN, 1, fir2c, fir2c);
                case 0:
                default :
                // creates fir1c matched sin(x)/x filter w blackman
                        wsincfilt(fir1c, 1.0 / symbollen, true);
                        fir1->init(FIRLEN, symbollen / 16, fir1c, fir1c);
                // creates fir2c matched sin(x)/x filter w blackman
                        wsincfilt(fir2c, 1.0 / 16.0, true);
                        fir2->init(FIRLEN, 1, fir2c, fir2c);
        }

        snfilt = new Cmovavg(16);
        imdfilt = new Cmovavg(16);

        tx_shape = new double[symbollen];

        // raised cosine shape for the transmitter
        for ( int i = 0; i < symbollen; i++)
                tx_shape[i] = 0.5 * cos(i * M_PI / symbollen) + 0.5;

        samplerate = PskSampleRate;
        fragmentsize = symbollen;
        bandwidth = samplerate / symbollen;
        snratio = s2n = imdratio = imd = 0;

        sigsearch = 0;
        for (int i = 0; i < 16; i++)
                syncbuf[i] = 0.0;
        E1 = E2 = E3 = 0.0;
        acquire = 0;

	restart();
}

// Demodulate -> Decode -> Play/Display
int fdmdv::rx_process(const double *buf, int len)
{
	double wbuf = buf[0];
	double* bufptr = &wbuf;
	vscard->Write(bufptr, len);
	return 0;
}

//=====================================================================
// fdmdv transmit
//=====================================================================

int fdmdv::tx_process()
{
	// Voice frame is 40mS
	// By a marvellous coincidence, this is 320 samples...
	// ... exactly the same as OPENLPC_FRAMESIZE_1_4
	int len = OPENLPC_FRAMESIZE_1_4;
	float voice_buffer[ len ];
	short sbuffer[ len ];
	unsigned char encodedbuffer[ FDMDV_CODEC_BITS + FDMDV_DATA_BITS ]; // One bit is stored in each char
	int encodedlen;
	// The modulated buffer is sent to the rig soundcard
	double modulatedbuffer[ len ];

	openlpc_encoder_state* voice_coder;
	voice_coder = create_openlpc_encoder_state();
	init_openlpc_encoder_state(voice_coder, OPENLPC_FRAMESIZE_1_4);
	// Get a 40ms voice sample and convert from floats to shorts
	vscard->Read(voice_buffer, len);
	for(int i=0 ; i < len ; i++) {
		sbuffer[i] = (short) voice_buffer[i];
	}

	// Frames are created in pairs.
	// First frame is 28 bits of voice
	// Second frame is 26 bits of voice + 2 bits data
	// First frame is indicated by no phase change on the BPSK carrier
	// Second frame is indicated by phase change on BPSK carrier
	
	// Encode the voice
	encodedlen = openlpc_encode(sbuffer, encodedbuffer, voice_coder);
	if(encodedlen != 55) {
	  // We've got a problem!
  }
	// This fills encodedbuffer[0] through encodedbuffer[55] with 1 bit of voice data
	// per element

	// Add the next 2 data bits into encodedbuffer[56] and [57]
  encodedbuffer[56] = varicoded_message[message_pointer++];
  encodedbuffer[57] = varicoded_message[message_pointer++];
	if (message_pointer >= message_length) message_pointer = 0;

	// Create carriers for 2 frames.
	// BPSK carrier first, then QPSK two at a time
	// I'm assuming encodedbuffer[0] is the LSB, and that the data bits
	// go at MSB and MSB-1 of the second frame
	// The first frame goes into the first half of modulatedbuffer
	// The second frame goes into the first half
	for (int i = 0; i < 2; i++) {
		fdmdv_write_bpsk(i == 0 ? FDMDV_FRAME_0:FDMDV_FRAME_1, FDMDV_SPACING*7, (double*) (modulatedbuffer + i*len/2), len/2);
		for(int carrier=0 ; carrier < 7 ; carrier++) {
			int freq_offset = FDMDV_SPACING*carrier; // Carrier spacing is 75Hz
			char frame_offset = 28*i;
			char symbol = (encodedbuffer[ 2*carrier + frame_offset ] << 1) +
			               encodedbuffer[ 2*carrier + frame_offset + 1 ];
			fdmdv_write_qpsk(symbol, freq_offset, (double*) (modulatedbuffer + i*len/2), len/2);
			symbol = (encodedbuffer[ 2*(carrier+7) + frame_offset ] << 1) +
			          encodedbuffer[ 2*(carrier+7) + frame_offset + 1 ];
			freq_offset += FDMDV_SPACING * 8;
    	fdmdv_write_qpsk(symbol, freq_offset, (double*) (modulatedbuffer + i*len/2), len/2);
  	}
	}

	ModulateXmtr(modulatedbuffer, len);
	return 0;
}

int varicode_encode(char* decoded, char* encoded)
{
	// Dummy function
  return 0;
}

int varicode_decode(char* decoded, char* encoded)
{
	// Dummy function
  return 0;
}

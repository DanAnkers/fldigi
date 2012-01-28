// ----------------------------------------------------------------------------
// fdmdv.cxx  --  fdmdv modem
//
// Copyright (C) 2012
//		Daniel Ankers, MD1CLV
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
#include "fdmdvvaricode.h"

#include "debug.h"


#define FDMDV_VOICE_MS		40	//40ms voice samples are passed to the codec
#define FDMDV_BW					1125
#define FDMDV_DATA_BITS		2		//2 data bits per frame
#define FDMDV_CODEC_BITS	56	//56 bits come back from the 40ms voice
#define FDMDV_CODEC_BYTES ceil(FDMDV_CODEC_BITS/8)
#define FDMDV_SPACING			75	//75Hz spacing between carriers
#define BPSK_CARRIER			7
#define BPSK							0
#define QPSK							1

using namespace std;

void fdmdv::tx_init(SoundBase *sc)
{
	scard = sc;
	
	// message_length = fdmdv_varicode_encode(message, varicoded_message);
	message_length = 10;
	message_pointer = 0;
}

void fdmdv::voicetx_init()
{
}

void fdmdv::rx_init()
{
	fill(phaseacc, phaseacc+CARRIERS, 0);
	fill(prevsymbol, prevsymbol+CARRIERS, complex (1.0, 0.0));
	quality         = complex (0.0, 0.0);
	//fill(shreg, shreg+CARRIERS, 0);
	dcdshreg = 0;
	dcd = 0;
	fill(bitclk, bitclk+CARRIERS, 0);
	freqerr = 0.0;
	sigsearch = 0;
	put_MODEstatus(mode);
	//resetSN_IMD();
	imdValid = false;
	afcmetric = 0.0;
	bpsk_bit = 0;
	fifo_read_ptr = 0;
	fifo_write_ptr = 0;
	fill(data_fifo, data_fifo+(CARRIERS-1)*5, 0);
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

	// PSK stuff

	cap |= CAP_AFC | CAP_AFC_SR;

	// QPSK50
	symbollen = 160;
	dcdbits = 50; //?
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
	double delta;
	complex z, z2;

	voice_decoder = create_openlpc_decoder_state();
	init_openlpc_decoder_state(voice_decoder, OPENLPC_FRAMESIZE_1_4);

	while (len-- > 0) {
		for(int carrier = 0; carrier < 15; carrier++)
		{
			delta = TWOPI * (frequency + carrier * FDMDV_SPACING) / samplerate;
			// Mix with the internal NCO
			z = complex ( *buf * cos(phaseacc[carrier]), *buf * sin(phaseacc[carrier]) );

			phaseacc[carrier] += delta;
			if (phaseacc[carrier] > M_PI)
				phaseacc[carrier] -= TWOPI;

			// Filter and downsample
			if (fir1->run( z, z )) { // fir1 returns true every Nth sample
				// final filter
				fir2->run( z, z2 ); // fir2 returns value on every sample
				//calcSN_IMD(z);

				// Sync correction
				int idx = (int) bitclk[carrier];
				double sum = 0.0;
				double ampsum = 0.0;
				syncbuf[idx] = 0.8 * syncbuf[idx] + 0.2 * z2.mag();

				for (int i = 0; i < 8; i++) {
					sum += (syncbuf[i] - syncbuf[i+8]);
					ampsum += (syncbuf[i] + syncbuf[i+8]);
				}
				// added correction as per PocketDigi
				sum = (ampsum == 0 ? 0 : sum / ampsum);

				bitclk[carrier] -= sum / 5.0;
				// Bit Clock - this should be synchronised over all carriers (?)
				// but we track per carrier just in case
				bitclk[carrier] += 1;

				if (bitclk[carrier] < 0) bitclk[carrier] += 16.0;
				if (bitclk[carrier] >= 16.0) {
					bitclk[carrier] -= 16.0;
					data_fifo[fifo_write_ptr++] = rx_symbol(z2, carrier==BPSK_CARRIER?BPSK:QPSK, carrier);
					if(fifo_write_ptr == fifo_read_ptr)
					{ //Buffer overflow
						printf("Overflow\n");
					}
					if(fifo_write_ptr == CARRIERS*BUFFER_FRAMES)
						fifo_write_ptr = 0;
				}
			}
		}
		buf++;
	}
	fifo_process();

	return 0;
}

void fdmdv::fifo_process(void)
{
	/* 
	Run through the FIFO which contains raw symbols and assemble them
	into FDMDV frame pairs (made up of FDMDV_CODEC_BITS voice bits and
	FDMDV_DATA_BITS data bits.)

	The first time this is called it is possible that we are missing the
	first frame of data - in that case the frame will be discarded.
	Following that if we receive half a frame pair then we leave fifo_read_ptr
	pointing to the start of the frame and deal with it the next time
	fifo_process is called.
	*/

	bitset<FDMDV_CODEC_BITS + FDMDV_DATA_BITS> encodedbits;
	unsigned char encodedbuffer[8];
	short sbuffer[OPENLPC_FRAMESIZE_1_4];
	double voice_buffer[OPENLPC_FRAMESIZE_1_4];
	int decodedlen;

	// Work out how many symbols we've got, taking into account buffer wraparound
	int symbols_to_process = 
		(fifo_write_ptr<fifo_read_ptr?CARRIERS*BUFFER_FRAMES:0)
		+ fifo_write_ptr-fifo_read_ptr;

	if((symbols_to_process%CARRIERS) != 0)
	{
		// This shouldn't happen - it means that either the write buffer or
		// the read buffer is in the wrong place
		// Log an error and return
		fprintf("Buffer appears to be in the wrong place - we have %d symbols\n", symbols_to_process);
		return;
	}
	if(symbols_to_process < 2*CARRIERS)
	{
		// We've been called without a full frame of data to process!
		printf("We don't have a full frame of data\n");
		return;
	}

	while(fifo_read_ptr != fifo_write_ptr)
	{
		// Find the first full frame
		// Indicated by a change in value on BPSK_CARRIER
		int first_bpsk_pos = fifo_read_ptr+BPSK_CARRIER;
		int second_bpsk_pos = (fifo_read_ptr >= CARRIERS*(BUFFER_FRAMES-1)?BPSK_CARRIER:fifo_read_ptr+CARRIERS+BPSK_CARRIER);
		if(data_fifo[first_bpsk_pos] != data_fifo[second_bpsk_pos])
		{
			// This looks like a full frame!
			int bitset_pos = 0;
			int tmp_read_ptr = fifo_read_ptr;

			// Put the data into a bitset
			for(int i = 0; i < CARRIERS*2; i++)
			{
				if(i%CARRIERS != BPSK_CARRIER)
				{
					encodedbits[bitset_pos++] = (data_fifo[tmp_read_ptr] & 1);
					encodedbits[bitset_pos++] = (data_fifo[tmp_read_ptr] & 2) >> 1;
				}
				// Check for FIFO wraparound
				if(++tmp_read_ptr > CARRIERS*BUFFER_FRAMES) tmp_read_ptr = 0;
			}

			// Convert the bitset into an array
			for(unsigned int i = 0; i < encodedbits.size(); i++)
			{
				int byte = floor(i/8); // Byte 0 is the most significant
				int bit = 7-(i%8);     // Bit 7 is the most significant
				encodedbuffer[byte] |= encodedbits[i] << bit; 
			}
			char datasymbol = encodedbuffer[7] & 3;
			if(datasymbol == 0)
			{ // We have received an entire character
				fdmdv_varicode_decode(shreg);
			} else {
				shreg = (shreg << 2) | datasymbol;
			}
			encodedbuffer[7] &= ~3;
			decodedlen = openlpc_decode(encodedbuffer, sbuffer, voice_decoder);
	
			for (int i = 0; i < OPENLPC_FRAMESIZE_1_4; i++)
			{
				voice_buffer[i] = (double)sbuffer[i];
			}
			vscard->Write(voice_buffer, OPENLPC_FRAMESIZE_1_4);
			fifo_read_ptr += CARRIERS;
		}
		fifo_read_ptr += CARRIERS;
		if(fifo_read_ptr >= CARRIERS*BUFFER_FRAMES) fifo_read_ptr = 0;
	}
	return;
}

unsigned char fdmdv::rx_symbol(complex symbol, char type, int carrier)
{
	unsigned char bits = 0;

	phase[carrier] = (prevsymbol[carrier] % symbol).arg();
	prevsymbol[carrier] = symbol;

	if (phase[carrier] < 0)
		phase[carrier] += TWOPI;

	if (type==QPSK) {
		bits = ((int) (phase[carrier] / M_PI_2 + 0.5)) & 3;
	} else { // bpsk
		bits = ((int) (phase[carrier] / M_PI + 0.5)) & 1;
	}

	return bits;
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
	unsigned char encodedbuffer[ OPENLPC_ENCODED_FRAME_SIZE ];
	int encodedlen;
	// The modulated buffer is sent to the rig soundcard
	double modulatedbuffer[ len ];

	voice_encoder = create_openlpc_encoder_state();
	init_openlpc_encoder_state(voice_encoder, OPENLPC_FRAMESIZE_1_4);
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
	encodedlen = openlpc_encode(sbuffer, encodedbuffer, voice_encoder);
	if(encodedlen != OPENLPC_ENCODED_FRAME_SIZE) {
	  // We've got a problem!
  }
	// This fills encodedbuffer[0] through encodedbuffer[6] with voice data
	// Convert these into a bitset:
	bitset<FDMDV_CODEC_BITS + FDMDV_DATA_BITS> encodedbits;
	for(unsigned int i = 0; i < encodedbits.size(); i++)
	{
		int byte = floor(i/8); // Byte 0 is the most significant
		int bit = 7-(i%8);     // Bit 7 is the most significant
		encodedbits[i] = (encodedbuffer[byte]>>bit) & 1;
	}

	// Add the 2 data bits into the LSBs of encodedbits
  encodedbits[FDMDV_CODEC_BITS] = varicoded_message[message_pointer++];
  encodedbits[FDMDV_CODEC_BITS + 1] = varicoded_message[message_pointer++];
	if (message_pointer >= message_length) message_pointer = 0;

	// Create carriers for 2 frames.
	// BPSK carrier first, then QPSK two at a time
	// The first frame goes into the first half of modulatedbuffer
	// The second frame goes into the second half
	for (int i = 0; i < 2; i++) {
		fdmdv_write_bpsk(i == 0 ? FDMDV_FRAME_0:FDMDV_FRAME_1, FDMDV_SPACING*7, (double*) (modulatedbuffer + i*len/2), len/2);
		for(int carrier=0 ; carrier < 7 ; carrier++) {
			int freq_offset = FDMDV_SPACING*carrier; // Carrier spacing is 75Hz
			char frame_offset = 28*i;

			// Each symbol is 2 bits
			// Carrier 0 uses bits 0 and 1
			// Carrier 1 uses bits 2 and 3
			// etc.
			char symbol = (char)((encodedbits[2*carrier+frame_offset]<<1)
			                     +encodedbits[2*carrier+frame_offset+1]);
			fdmdv_write_qpsk(symbol, freq_offset, (double*) (modulatedbuffer + i*len/2), len/2);
			symbol = (char)((encodedbits[2*(carrier+7)+frame_offset]<<1)
			                     +encodedbits[2*(carrier+7)+frame_offset+1]);
			freq_offset += FDMDV_SPACING * 8;
    	fdmdv_write_qpsk(symbol, freq_offset, (double*) (modulatedbuffer + i*len/2), len/2);
  	}
	}

	ModulateXmtr(modulatedbuffer, len);
	return 0;
}

int varicode_encode(char* decoded, char* encoded)
{
	char encoded_char;
	int len=0;
	int encoded_char_bitcount;
	for (unsigned int i=0; i < sizeof(decoded); i++)
	{
		encoded_char=fdmdv_varicode_encode(decoded[i]);
		if(encoded_char<3)
			encoded_char_bitcount = 2;
		else if(encoded_char<15)
			encoded_char_bitcount = 4;
		else
			encoded_char_bitcount = 6;
		// Right shift "encoded" array.
		// (How?)
		// Add 2 extra bits for end of character signal
		int correct_byte = 1; // FIXME
		encoded[correct_byte] &= !((2^(encoded_char_bitcount+2))-1);
		encoded[correct_byte] |= (encoded_char<<encoded_char_bitcount); 
		len+=encoded_char_bitcount+2;
	}
	return len;
}


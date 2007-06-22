#ifndef _CONFIGURATION_H
#define _CONFIGURATION_H

#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "main.h"
#include "rtty.h"

using namespace std;

struct configuration {
	bool	changed;
	double 	squelch;
	double	wfRefLevel;
	double	wfAmpSpan;
	double	CWsweetspot;
	double	RTTYsweetspot;
	double	PSKsweetspot;
	bool	StartAtSweetSpot;
// RTTY
	double		rtty_squelch;
	int			rtty_shift;
	int			rtty_baud;
	int 		rtty_bits;
	int			rtty_parity;
	int			rtty_stop;
	bool 		rtty_reverse;
	bool		rtty_msbfirst;
	bool		rtty_crcrlf;
	bool		rtty_autocrlf;
	int			rtty_autocount;
	int			rtty_afcspeed;
	bool		useFSKkeyline;		// use RTS for FSK
	bool		useFSKkeylineDTR;	// use DTR for FSK
	bool		FSKisLSB;
	bool		RTTY_USB;
	bool		useUART;
	bool		PreferXhairScope;
// CW
	bool		useCWkeylineRTS;	// use RTS for CW
	bool		useCWkeylineDTR;	// use DTR for CW
	int			CWweight;
	int			CWspeed;
	int			CWbandwidth;
	int			CWtrack;
	int			CWrange;
// FELD-HELL
	bool		FELD_IDLE;
// OLIVIA
	int			oliviatones;
	int			oliviabw;
// DOMINOEX
	double		DOMINOEX_BW;
// User interface data
	int		Font;
	int		FontSize;
	int		FontColor;
	uchar	red;
	uchar	green;
	uchar	blue;
	bool	MultiColorWF;
	int		wfPreFilter;
	bool	UseCursorLines;
	bool	UseCursorCenterLine;
	bool	UseBWTracks;
	int		feldfontnbr;
	bool	viewXmtSignal;
	bool	sendid;
	bool	macroid;
	int		QRZ;
// Rig Interface data
	bool	btnusb;
	int		btnPTTis;
	int		btnRTSDTRis; // obsolete
	int		btnPTTREVis; // obsolete
	bool	RTSptt;
	bool	DTRptt;
	bool	RTSplus;
	bool	DTRplus;
	int		choiceHAMLIBis;
	int		chkUSEMEMMAPis;
	int		chkUSEHAMLIBis;
	int		chkUSERIGCATis;
	string  HamRigName;
	string  HamRigDevice;
	int		HamRigBaudrate;
	string	CWFSKport;
// Operator data
	string	myCall;
	string	myQth;
	string	myName;
	string	myLocator;
	string  PTTdev;
	string	secText;
// Sound card
	string	SCdevice;
	int		RX_corr;
	int		TX_corr;
// Contest stuff
	bool	UseLeadingZeros;
	int		ContestStart;
	int		ContestDigits;
// Macro timer constants and controls
	bool	useTimer;
	int		macronumber;
	int		timeout;
	
// Mixer configuration
	double	RcvMixer;
	double	XmtMixer;
	bool	MicIn;
	bool	LineIn;
	bool	EnableMixer;
	double	PCMvolume;
	bool	MicMute;
	bool	LineMute;
	
// waterfall palette
	RGBint	cfgpal[9];
	
public:
	void readDefaults(ifstream &f);
	void writeDefaults(ofstream &f);
	void storeDefaults();
	void loadDefaults();
	void saveDefaults();
	int  openDefaults();
	void initOperator();
	void initInterface();
	void initMixerDevices();
	void writeDefaultsXML();
	
	void getRigs();
	string strBaudRate();
	
	friend std::istream &operator>>(std::istream &stream, configuration &c);
	friend std::ostream &operator<<(std::ostream &ostream, configuration c);
};

extern configuration progdefaults;

extern void mixerInputs();
extern void enableMixer(bool);

#endif
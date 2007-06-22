#ifndef MIXER_H
#define MIXER_H

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>

#include <sys/ioctl.h>
#include <fcntl.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sys/soundcard.h>
#include <math.h>

#include <iostream>
#include <string>

class MixerException {
public:
	char	szError[80];
	int		error;
	MixerException() { *szError = 0; error = 0; }
	MixerException(int e) {
		sprintf(szError,"Error: %d, %s", e, strerror(e));
		error = e;
	}
	MixerException(char *s) {
		sprintf(szError,"Error: %s", s);
		error = 1;
	}
};

class cMixer {
private:
	std::string mixer;
	int		mixer_fd;
	int		recmask;
	int		devmask;
	int		outmask;
	int		recsrc;
	int		num_out;
	int		outs[SOUND_MIXER_NRDEVICES];
	int		num_rec;
	int		recs[SOUND_MIXER_NRDEVICES];

// values on init for restoration	
	int		inpsrc0;
	double	inplevel0;
	double	pcmlevel0;
	double	vollevel0;

	int		NumMixers;
	int		NumDevice;
	int		Devices[10];
	char	szDevice[12];

	int		initMask();
	void	findNumMixers();
	double	ChannelVolume(int);

public:
	cMixer();
	~cMixer();
	void			openMixer(const char *dev = "/dev/mixer");
	void			closeMixer();
	
	void			setXmtLevel(double v);
	void			setRcvGain(double v);
	
	int				numMixers() { return NumMixers;}
	int				MixerNum(int i) { return Devices[i];}
	const char *	MixerName( int index );
	double			OutVolume();
	void			OutVolume(double vol);
	double			PCMVolume();
	void			PCMVolume(double volume );
	int				NumOutputVolumes();
	double			OutputVolume( int i );
	void			OutputVolume( int i, double volume );
	const char *	OutputVolumeName( int i );
	int				GetNumInputSources();
	const char *	GetInputSourceName( int i);
	int				InputSourceNbr(char *source);
	double			InputVolume();
	void			InputVolume( double volume );
	int				GetCurrentInputSource();
	void			SetCurrentInputSource( int i );
	double			GetPlaythrough();
	void			SetPlaythrough( double volume );
	
};

#endif
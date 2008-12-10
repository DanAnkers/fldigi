#ifndef DXCC_H_
#define DXCC_H_

struct dxcc {
	const char* country;
	int cq_zone;
	int itu_zone;
	char continent[3];
	float latitude;
	float longitude;
	float gmt_offset;
	dxcc(const char* cn = "", int cq = 0, int itu = 0, const char* ct = "",
	     float lat = 0.0f, float lon = 0.0f, float tz = 0.0f);
};

bool dxcc_open(const char* filename);
void dxcc_close(void);
const dxcc* dxcc_lookup(const char* callsign);

#endif // DXCC_H_

#pragma once
#include <stdint.h>
#define GMUL65536(x) ( (( (x&32768) == 0) ? x<<1 : (x<<1 ^ 0x2b ) ) & 0xffff )
static uint16_t multiplyInGF65536(
	uint16_t a,
	uint16_t b
) {
	uint16_t res = 0;
	for (; b; b >>= 1) {
		if (b & 1) {
			res ^= a;
		}
		a = GMUL65536(a);
	}
	return res;
}
#define GMUL256(x) ( (( (x&128) == 0) ? x<<1 : (x<<1 ^ 0x1b ) ) & 0xff )
static uint8_t multiplyInGF256_0x1b(
	uint8_t a,
	uint8_t b
) {
	uint8_t res = 0;
	for (; b; b >>= 1) {
		if (b & 1) {
			res ^= a;
		}
		a = GMUL256(a);
	}
	return res;
}


#define GMUL16(x) ( ( (( x & 8 ) == 0) ? x<<1 : ( x<<1 ^ 0x3 ) ) & 0xf )
static uint8_t multiplyInGF16(
	uint8_t a,
	uint8_t b
) {
	uint8_t res = 0;
	for (; b; b >>= 1) {
		if (b & 1) {
			res ^= a;
		}
		a = GMUL16(a);
	}
	return res;

}
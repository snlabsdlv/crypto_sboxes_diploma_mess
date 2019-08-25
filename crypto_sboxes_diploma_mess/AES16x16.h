typedef unsigned int uint32_t;
typedef unsigned short uint16_t;
typedef unsigned char uint8_t;
/**
* Modular inverse in GF(2ˆ16) using the EEA algorithm.
*/
#define PX_16 0x002B
#define FPX_16 0x1002B
#define MSB_16 0x8000
#define HMSB_16 0x10000
#define MSB_8 0x80
#define HMSB_8 0x100
#define LSB 0x1

// Quotient and remainder struct
typedef struct
{
	uint16_t q;
	uint16_t r;
	uint8_t error;
} QR;
/**
* Polynomial division in GF(2ˆ16).
*/
/**
* Polynomial division in GF(2ˆ16).
*/
QR g16_div(uint32_t ai, uint16_t b)
{
	uint16_t a = (uint16_t)ai;
	int msb = MSB_16;
	int d = 0;
	QR result = { 0, 0 };
	// Align the denominator with the numerator
	while (b > 0 && !(b & MSB_16)) {
		++d;
		b <<= 1;
	}
	// If the polynomial MSB is set (17th bit), increment
	// the quotient and reduce the numerator.
	if (ai & HMSB_16) {
		result.q ^= 1 << (d + 1);
		a ^= b << 1;
	}
	for (; d > -1; d--) {
		if ((a & msb) && (b & msb)) {
			result.q ^= 1 << d;
			a ^= b;
		}
		msb >>= 1;
		b >>= 1;
	}
	result.r = a;
	return result;
}
uint16_t g16_add(uint16_t x, uint16_t y)
{
	return x ^ y;
}
uint16_t g16_mul(uint16_t x, uint16_t y)
{
	uint16_t accum = 0;
	uint16_t msb = 0;
	uint16_t i;
	for (i = 0; i < 16; i++)
	{
		if (y & LSB) accum ^= x;
		msb = (x & MSB_16); // fetch the MSB
		x <<= 1;
		if (msb) x ^= PX_16;
		y >>= 1;
	}
	return accum;
}

uint16_t g16_inv(uint16_t x)
{
	// Trivial special cases.
	if (x == 0) return 0;
	if (x == 1) return 1;
	uint16_t r0 = PX_16; // rem[i - 2]
	uint16_t r1 = x; // rem[i - 1]
	uint16_t a0 = 0; // aux[i - 2]
	uint16_t a1 = 1; // aux[i - 1]
	uint16_t tmp;
	QR qr;
	int firstRun = 0;
	while (r1 > 0)
	{
		if (firstRun != 0) qr = g16_div(r0, r1);
		else
		{
			qr = g16_div(FPX_16, r1);
			firstRun++;
		}
		r0 = r1; r1 = qr.r;
		tmp = a0; a0 = a1;
		a1 = g16_add(tmp, g16_mul(qr.q, a1));
	}
	return a0;
}
uint16_t g16_change_basis(uint16_t x, uint16_t* M)
{
	uint16_t y = 0;
	for (int i = 15; i >= 0; i--)
	{
		if (x & 1) y ^= M[i];
		x >>= 1;
	}
	return y;
}

static uint16_t A[16] =
{
	0x797b, 0x7c85, 0x9378, 0x151, 0x2312, 0x82f, 0x3f35, 0xe57e,
	0x29d, 0x7e12, 0xdc62, 0xadbb, 0xced3, 0x87a0, 0xe900, 0x2d9c
};
static uint16_t AINV[16] =
{
	0x6a5e, 0xc863, 0x3b62, 0xec10, 0x3931, 0xb56e, 0xd1e7, 0xa06c,
	0x585f, 0x230c, 0xf6e0, 0x5557, 0x577e, 0x4d26, 0x17be, 0xc637
};
// Constant from our new S-box
#define CC 0x45B7
uint16_t forward(uint16_t x)
{
	uint16_t inv = g16_inv(x);
	inv = g16_change_basis(inv, A);
	return inv ^ CC;
}
uint16_t inverse(uint16_t x)
{
	uint16_t inv = x ^ CC;
	inv = g16_change_basis(inv, AINV);
	return g16_inv(inv);
}
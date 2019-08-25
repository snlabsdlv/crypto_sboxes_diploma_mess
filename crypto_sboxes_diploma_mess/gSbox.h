#pragma once
#define ROR(x, r) ((x >> r) | (x << (8 - r)))
#define ROL(x, r) ((x << r) | (x >> (8 - r)))

#define g1(x, y) (x = ROR(x, 2), x += y, y = ROL(y, 1)+185,y ^= x)
#define g2(x, y) (x = ROL(x, 1), x += y, y = ROR(y, 2)+100,y ^= x)
#define g3(x, y) (x = ROL(x, 1), x += y, y = ROR(y, 2),y ^= x)
#define g4(x, y) (x = ROR(x, 2), x += y, y = ROL(y, 1),y ^= x)

unsigned short speck2(unsigned short plain)
{
	unsigned char left = plain >> 8, right = plain;
	g1(left, right);
	g1(left, right);
	g1(left, right);
	g1(left, right);
	g1(left, right);
	g1(left, right);
	return (left << 8) | right;
}
void speck2ptr(unsigned char plain[2])
{
	for (int i = 0; i < 6; i++)
		g1(plain[0], plain[1]);
}
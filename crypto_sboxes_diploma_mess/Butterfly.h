#pragma once
#include <math.h>
#include <stdint.h>
#include "Precomp_Pi16.h"
#include "Precomp_Pi256.h"

uint8_t butt8Jimenez(uint8_t in, uint8_t pi1[16], uint8_t pi2[16]) {
	uint8_t low4 = in & 0xf;
	uint8_t high4 = (in & 0xf0) >> 4;
	low4 = high4 ? g14[multiplyInGF16(low4, high4)] : pi1[low4];
	high4 = low4 ? multiplyInGF16(low4, g14[high4]) : pi2[high4];
	return (high4<<4) | low4;
}
uint8_t butt8FominB1(uint8_t in) {
	uint8_t low4 = in & 0xf;
	uint8_t high4 = (in & 0xf0) >> 4;
	low4 = high4 ? multiplyInGF16(low4, high4) : g14[low4];
	high4 = low4 ? multiplyInGF16(high4, g13[low4]) : g14[high4];
	return (high4 << 4) | low4;
}
uint8_t butt8FominB2(uint8_t in) {
	uint8_t low4 = in & 0xf;
	uint8_t high4 = (in & 0xf0) >> 4;
	low4 = high4 ? multiplyInGF16(low4, g2[high4]) : g14[low4];
	high4 = low4 ? multiplyInGF16(high4, g14[low4]) : g14[high4];
	return (high4 << 4) | low4;
}

uint16_t butt16Fomin(uint16_t in, uint8_t pi1[256], uint8_t pi2[256]) {
	uint8_t low8 = in&0xff;
	uint8_t high8 = in>>8;
	low8 = high8 ? multiplyInGF256_0x1b(low8, pi1[high8]) : g254_256[low8];
	high8 = low8 ? multiplyInGF256_0x1b(high8, pi2[low8]) : g254_256[high8];

	return (high8 << 8) | low8;
}

uint16_t butt16Jimenez(uint16_t in, uint8_t pi1[256], uint8_t pi2[256]) {
	uint8_t low8 = in & 0xff;
	uint8_t high8 = in >> 8;
	low8 = high8 ? g254_256[multiplyInGF256_0x1b(low8, high8)] : pi1[low8];
	high8 = low8 ? multiplyInGF256_0x1b(low8, g254_256[high8]) : pi2[high8];

	return (high8 << 8) | low8;
}

uint16_t butt16test(uint16_t in, int k, uint8_t pi1[256], uint8_t pi2[256]) {
	uint8_t low8 = in & 0xff;
	uint8_t high8 = in >> 8;
	low8 = multiplyInGF256_0x1b(k, high8);
	high8 = low8+pi2[high8];
	return (high8 << 8) | low8;
}

uint16_t butt16from4(uint16_t in, uint8_t pi1[16], uint8_t pi2[16]) {
	uint8_t b0 = in & 0xf;
	uint8_t b1 = (in & 0xf0) >> 4;
	uint8_t b2 = (in & 0xf00) >> 8;
	uint8_t b3 = (in & 0xf000) >> 12;

	b0 = b1 ? pi1[multiplyInGF16(b0, b1)] : pi1[b0];
	b1 = b2 ? multiplyInGF16(b1, pi2[b2]) : pi2[b1];
	b2 = b3 ? pi1[multiplyInGF16(b2, b3)] : pi1[b2];
	b3 = b0 ? multiplyInGF16(b3, pi2[b0]) : pi2[b3];

	return (uint16_t)((b3 << 12) | (b2 << 8) | (b1 << 4) | b0);
}

uint16_t butt16from4_2(uint16_t in, uint8_t pi1[16], uint8_t pi2[16]) {
	uint8_t b0 = in & 0xf;
	uint8_t b1 = (in & 0xf0) >> 4;
	uint8_t b2 = (in & 0xf00) >> 8;
	uint8_t b3 = (in & 0xf000) >> 12;

	b0 = b1 ? pi1[multiplyInGF16(b0, b1)] : pi1[b0];
	b1 = b0 ? multiplyInGF16(b1, pi2[b0]) : pi2[b1];
	b2 = b3 ? pi1[multiplyInGF16(b2, b3)] : pi1[b2];
	b3 = b2 ? multiplyInGF16(b3, pi2[b2]) : pi2[b3];
	b0 = b3 ? pi1[multiplyInGF16(b0, b3)] : pi1[b0];
	b3 = b0 ? multiplyInGF16(b3, pi2[b0]) : pi2[b3];

	return (uint16_t)((b3 << 12) | (b2 << 8) | (b1 << 4) | b0);
}
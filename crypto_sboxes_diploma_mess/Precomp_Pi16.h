#pragma once
#include <stdint.h>
// modulo 0x3

uint8_t g2[16] = { 0, 1, 4, 5, 3, 2, 7, 6, 12, 13, 8, 9, 15, 14, 11, 10 };
uint8_t g3[16] = { 0, 1, 8, 15, 12, 10, 1, 1, 10, 15, 15, 12, 8, 10, 8, 12 };
uint8_t g4[16] = { 0, 1, 3, 2, 5, 4, 6, 7, 15, 14, 12, 13, 10, 11, 9, 8 };
uint8_t g5[16] = { 0, 1, 6, 6, 7, 7, 7, 6, 1, 7, 1, 6, 1, 6, 7, 1 };
uint8_t g6[16] = { 0, 1, 12, 10, 15, 8, 1, 1, 8, 10, 10, 15, 12, 8, 12, 15 };
uint8_t g7[16] = { 0, 1, 11, 13, 9, 14, 6, 7, 12, 5, 8, 3, 15, 2, 4, 10 };
uint8_t g8[16] = { 0, 1, 5, 4, 2, 3, 7, 6, 10, 11, 15, 14, 8, 9, 13, 12 };
uint8_t g9[16] = { 0, 1, 10, 12, 8, 15, 1, 1, 15, 12, 12, 8, 10, 15, 10, 8 };
uint8_t g10[16] = { 0, 1, 7, 7, 6, 6, 6, 7, 1, 6, 1, 7, 1, 7, 6, 1 };
uint8_t g11[16] = { 0, 1, 14, 9, 11, 13, 7, 6, 8, 3, 10, 4, 12, 5, 2, 15 };
uint8_t g12[16] = { 0, 1, 15, 8, 10, 12, 1, 1, 12, 8, 8, 10, 15, 12, 15, 10 };
uint8_t g13[16] = { 0, 1, 13, 11, 14, 9, 6, 7, 10, 4, 15, 2, 8, 3, 5, 12 };
uint8_t g14[16] = { 0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8 };

// modulo 0x9
/*
uint8_t g2[16] = { 0, 1, 4, 5, 9, 8, 13, 12, 15, 14, 11, 10, 6, 7, 2, 3, };
uint8_t g3[16] = { 0, 1, 8, 15, 15, 3, 5, 15, 5, 3, 1, 1, 3, 8, 5, 8, };
uint8_t g4[16] = { 0, 1, 9, 8, 14, 15, 7, 6, 3, 2, 10, 11, 13, 12, 4, 5, };
uint8_t g5[16] = { 0, 1, 11, 1, 10, 1, 11, 11, 1, 11, 11, 10, 10, 10, 10, 1, };
uint8_t g6[16] = { 0, 1, 15, 3, 3, 5, 8, 3, 8, 5, 1, 1, 5, 15, 8, 15, };
uint8_t g7[16] = { 0, 1, 7, 5, 12, 8, 2, 9, 15, 6, 10, 11, 14, 4, 13, 3, };
uint8_t g8[16] = { 0, 1, 14, 15, 2, 3, 12, 13, 5, 4, 11, 10, 7, 6, 9, 8, };
uint8_t g9[16] = { 0, 1, 5, 8, 8, 15, 3, 8, 3, 15, 1, 1, 15, 5, 3, 5, };
uint8_t g10[16] = { 0, 1, 10, 1, 11, 1, 10, 10, 1, 10, 10, 11, 11, 11, 11, 1, };
uint8_t g11[16] = { 0, 1, 13, 3, 7, 5, 14, 4, 8, 12, 11, 10, 9, 2, 6, 15, };
uint8_t g12[16] = { 0, 1, 3, 5, 5, 8, 15, 5, 15, 8, 1, 1, 8, 3, 15, 3, };
uint8_t g13[16] = { 0, 1, 6, 15, 13, 3, 9, 2, 5, 7, 10, 11, 4, 14, 12, 8, };
uint8_t g14[16] = { 0, 1, 12, 8, 6, 15, 4, 14, 3, 13, 11, 10, 2, 9, 7, 5, };
*/
// modulo 0xf
/*
uint8_t g2[16] = { 0, 1, 4, 5, 15, 14, 11, 10, 2, 3, 6, 7, 13, 12, 9, 8, };
uint8_t g3[16] = { 0, 1, 8, 15, 2, 8, 4, 8, 15, 4, 2, 15, 1, 1, 2, 4, };
uint8_t g4[16] = { 0, 1, 15, 14, 8, 9, 7, 6, 4, 5, 11, 10, 12, 13, 3, 2, };
uint8_t g5[16] = { 0, 1, 1, 13, 1, 12, 13, 13, 1, 12, 12, 12, 13, 12, 13, 1, };
uint8_t g6[16] = { 0, 1, 2, 8, 4, 2, 15, 2, 8, 15, 4, 8, 1, 1, 4, 15, };
uint8_t g7[16] = { 0, 1, 4, 7, 15, 10, 3, 14, 2, 11, 9, 5, 12, 13, 6, 8, };
uint8_t g8[16] = { 0, 1, 8, 9, 2, 3, 10, 11, 15, 14, 7, 6, 13, 12, 5, 4, };
uint8_t g9[16] = { 0, 1, 15, 4, 8, 15, 2, 15, 4, 2, 8, 4, 1, 1, 8, 2, };
uint8_t g10[16] = { 0, 1, 1, 12, 1, 13, 12, 12, 1, 13, 13, 13, 12, 13, 12, 1, };
uint8_t g11[16] = { 0, 1, 2, 11, 4, 7, 9, 5, 8, 6, 14, 3, 13, 12, 10, 15, };
uint8_t g12[16] = { 0, 1, 4, 2, 15, 4, 8, 4, 2, 8, 15, 2, 1, 1, 15, 8, };
uint8_t g13[16] = { 0, 1, 8, 6, 2, 11, 14, 3, 15, 10, 5, 9, 12, 13, 7, 4, };
uint8_t g14[16] = { 0, 1, 15, 10, 8, 6, 5, 9, 4, 7, 3, 14, 13, 12, 11, 2, };
*/

uint8_t* gs[13] = { g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14 };